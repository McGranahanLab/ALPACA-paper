#!/usr/bin/env Rscript

# To run this script for 'primary' cohort, and both LUAD and LUSC,
# run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/DRIVERORDERING_supplementaryanalysis.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUAD_all_topfreqevents.csv \
--all_events_path output/primary/primary_default/cohort_outputs/eventordering/LUAD_all_allevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin \
--save_directory output/primary/primary_default/cohort_outputs/eventordering
"

cat("\nRunning driver ordering supplementary analysis\n")

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gridExtra))

### Command line options 
option_list <- list(
  make_option(c("--top_events_path"), type="character", action="store", help="", metavar="character"),
  make_option(c("--all_events_path"), type="character", action="store", help="", metavar="character"),
  make_option(c("--cohort_tree_directory"), type="character", help="Directory where tree data is stored"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--save_directory"), type="character", help="save_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

# TODO remove
# opt$top_events_path <- '/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/publication/output/primary/primary_default/cohort_outputs/eventordering/LUSC_all_topfreqevents.csv'
# opt$all_events_path <- '/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/publication/output/primary/primary_default/cohort_outputs/eventordering/LUSC_all_allevents.csv'
# opt$cohort_tree_directory <- '/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/publication/_assets/tree_data/primary'
# opt$bin_directory <- '/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/publication/bin'

top_events_path <- opt$top_events_path
all_events_path <- opt$all_events_path
cohort_tree_directory <- opt$cohort_tree_directory
bin_dir <- opt$bin_directory
save_directory <- opt$save_directory

assets_dir <- gsub('bin', '_assets', bin_dir)
histo <- strsplit(basename(top_events_path), "_")[[1]][1]
clonality <- strsplit(basename(top_events_path), "_")[[1]][2]
histo_clonality <- paste(histo, clonality, sep = "_")
print(histo_clonality)

### Functions
source(file.path(bin_dir, "ANALYSIS/analysis_functions.R"))
source(file.path(bin_dir, "ANALYSIS/DRIVERORDERING_functions.R"))

### Colour palette
tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))

###################################################################################################
### Preprocess driver events - except now include all driver events, not just top most frequent ###
###################################################################################################

### Preprocess data
top_events <- fread(top_events_path)
all_events <- fread(all_events_path)

all_events <- merge.data.table(all_events, unique(top_events[, .(tumour_id, clone, parent_clone)]), by = c("tumour_id", "clone"), all = TRUE)

all_events <- all_events[!is.na(parent_clone)]

all_events[, clone := as.character(clone)]
all_events[, Clonality := ifelse(parent_clone == "diploid", "Truncal", "Subclonal")]
all_events[, n_tums_with_event := length(unique(tumour_id)), by = event_id]
all_events[, n_clones_with_event := .N, by = event_id]
all_events[, n_clones_with_event_clonality := .N, by = .(event_id, Clonality)]
all_events[is.na(event_id), n_clones_with_event := NA]
all_events[is.na(event_id), n_clones_with_event_clonality := NA]
all_events[is.na(event_id), n_tums_with_event := NA]

unique_all_events <- sort(unique(all_events$event_id))

### Get ordering of number of events across the cohort
all_events_counts <- unique(all_events[!is.na(event_id), .(event_id, n_clones_with_event, n_tums_with_event)])
setorder(all_events_counts, -n_clones_with_event)
all_events_prevalence_order <- unique(all_events_counts$event_id)

### Get ordering of number of events across the cohort, split by clonal/subclonal
all_events_counts_clonality <- unique(all_events[!is.na(event_id), .(event_id, n_clones_with_event_clonality, n_clones_with_event, Clonality, n_tums_with_event)])

### Get mean/max CCF of each clone in the tumour
all_events_clones <- unique(all_events[, .(tumour_id, clone)])
all_events_clones[, c("mean_ccf_tumour", "max_ccf_tumour", "max_ccf_expandedsubclone") := 
  rbindlist(lapply(seq(nrow(all_events_clones)), function(i) {
    # Load tree data
    tum_id <- all_events_clones[i, tumour_id]
    c <- all_events_clones[i, clone]
    tree_object <- readRDS(file = paste0(cohort_tree_directory, '/' , tum_id, '_trees.rds'))
    ccf_cluster_table <- tree_object$nested_pyclone$ccf_cluster_table
    clonality_table <- tree_object$clonality_out$clonality_table_corrected
    # Extract CCF metrics
    clone_ccf <- ccf_cluster_table[rownames(ccf_cluster_table) == c, ]
    clone_clonality <- clonality_table[rownames(clonality_table) == c, ]
    max_ccf_region <- names(clone_ccf)[which(clone_ccf == max(clone_ccf))]
    max_ccf_expandedsubclone <- ifelse("clonal" %in% as.character(clone_clonality[max_ccf_region]), "Expanded subclone", "Non-expanded subclone")
    out_df <- data.table(mean_ccf_tumour = mean(clone_ccf),
                        max_ccf_tumour = max(clone_ccf),
                        max_ccf_expandedsubclone = max_ccf_expandedsubclone)
    return(out_df)
}))]
# Merge with events:
all_events <- merge.data.table(all_events, all_events_clones, by = c("tumour_id", "clone"), all.x = TRUE)

### Label the driver gene events
all_events[, driver_alterations_clone := paste(event_id, collapse = ";"), by = .(tumour_id, clone)]
all_events[grepl("NA", driver_alterations_clone), driver_alterations_clone := NA]

### Remove tumours with no events:
all_events[, n_clones_w_event_tum := sum(!is.na(driver_alterations_clone)), by = tumour_id]
all_events <- all_events[n_clones_w_event_tum > 0]

# Assign event type
all_events[grepl('chr arm LOH', event_name), `Event type` := "Chr arm LOH"]
all_events[grepl('chr arm amp', event_name), `Event type` := "Chr arm amp"]
all_events[grepl('gene amp', event_name), `Event type` := "Gene amp"]
all_events[grepl('snv', event_name), `Event type` := "SNV"]


##############################################################################
### Extra analysis comparing distribution of cluster CCFs of driver events ###
##############################################################################

sc_expansion_levels <- c("Non-expanded subclone", "Expanded subclone")

# Extract all subclonal events
sc_all_events <- all_events[Clonality == 'Subclonal' & !is.na(event_id)]
sc_all_events[, n_clones_pereventtype := .N, by = .(`Event type`)]
sc_all_events[, n_clones_sc_expansion_pereventtype := .N, by = .(`Event type`, max_ccf_expandedsubclone)]
sc_all_events[, n_clones_sc_expansion_pereventid := .N, by = .(event_id, max_ccf_expandedsubclone)]

sc_all_events[, frac_clones_sc_expansion_pereventtype := n_clones_sc_expansion_pereventtype / .N, by = .(`Event type`)]
sc_all_events[, frac_clones_sc_expansion_pereventid := n_clones_sc_expansion_pereventid / .N, by = .(event_id)]

# Extract top sc events
sc_top_events <- sc_all_events[top_event == TRUE]

top_events_clonalillu_fishout <- rbindlist(lapply(unique(sc_top_events$event_id), function(cur_event) {
	print(cur_event)
	event_df <- unique(sc_all_events[event_id == cur_event, .(max_ccf_expandedsubclone,
																														n_clones_sc_expansion_pereventtype,
																														n_clones_sc_expansion_pereventid,
																														`Event type`)])
	cur_event_type <- event_df[, unique(`Event type`)]
	if (!all(sc_expansion_levels %in% event_df$max_ccf_expandedsubclone)) {
		missing_level <- sc_expansion_levels[!sc_expansion_levels %in% event_df$max_ccf_expandedsubclone]
		event_other <- data.table(
			max_ccf_expandedsubclone = missing_level,
			n_clones_sc_expansion_pereventtype = unique(sc_all_events[max_ccf_expandedsubclone == missing_level & `Event type` == unique(event_df$`Event type`), n_clones_sc_expansion_pereventtype]),
			n_clones_sc_expansion_pereventid = 0
		)
		event_df[, `Event type` := NULL]
		event_df <- rbind(event_df, event_other)
	}
	event_df[, max_ccf_expandedsubclone := factor(max_ccf_expandedsubclone, levels = sc_expansion_levels)]
	setorder(event_df, max_ccf_expandedsubclone)
	
	# Run fisher test:
	ct <- as.matrix(event_df)
	rownames(ct) <- ct[,1]
	ct <- ct[, 2:3]
	ct <- apply(ct, c(1, 2), as.numeric)
	ft <- fisher.test(ct)
	ft_out <- data.table(`Event type` = cur_event_type,
													event_id = cur_event,
													odds_ratio = ft$estimate,
													ci_lower = ft$conf.int[1],
													ci_upper = ft$conf.int[2],
													conf_level = attributes(ft$conf.int)[['conf.level']],
													pval = ft$p.value)
	return(ft_out)
}))
top_events_clonalillu_fishout[pval < 0.05, pval_signif := '*']
top_events_clonalillu_fishout[pval < 0.01, pval_signif := '**']
top_events_clonalillu_fishout[pval < 0.001, pval_signif := '***']

topevents_df <- unique(sc_top_events[, .(event_id,
									n_clones_with_event_clonality,
									max_ccf_expandedsubclone,
									`Event type`,
									n_clones_sc_expansion_pereventid,
									frac_clones_sc_expansion_pereventid)])

bgevents_df <- unique(sc_top_events[, .(event_id = 'background',
												n_clones_pereventtype,
												max_ccf_expandedsubclone,
												`Event type`,
												n_clones_sc_expansion_pereventtype,
												frac_clones_sc_expansion_pereventtype)])
setnames(bgevents_df,
					c("n_clones_sc_expansion_pereventtype", "frac_clones_sc_expansion_pereventtype", "n_clones_pereventtype"),
					c("n_clones_sc_expansion_pereventid", "frac_clones_sc_expansion_pereventid", 'n_clones_with_event_clonality'))

plot_df <- rbind(topevents_df, bgevents_df)
plot_df <- merge.data.table(plot_df, unique(top_events_clonalillu_fishout[, .(event_id, pval, pval_signif)]), by = "event_id", all = TRUE)
plot_df[, n_clones_event_plot := n_clones_with_event_clonality]
plot_df[, n_clones_sc_expansion_pereventid_plot := n_clones_sc_expansion_pereventid]

plot_df[event_id == "background", n_clones_event_plot := NA]
plot_df[event_id == "background", n_clones_sc_expansion_pereventid_plot := NA]

plot_df[, sc_expansion_factor := factor(max_ccf_expandedsubclone, levels = sc_expansion_levels)]

# plot_df[!is.na(pval_signif), event_name_sig := paste(event_id, pval_signif, sep = " ")]
# plot_df[is.na(pval_signif), event_name_sig := event_id]
plot_df[, event_name_sig := event_id]


plot_df[, frac_expanded := sum((max_ccf_expandedsubclone == 'Expanded subclone') * frac_clones_sc_expansion_pereventid), by = .(event_id, `Event type`)]
setorder(plot_df, -frac_expanded)
event_order <- unique(plot_df$event_name_sig)
event_order <- event_order[event_order != "background"]
event_order <- c(event_order, "background")
event_face <- rep('plain', length(event_order))
event_face[grepl("[ *]", event_order)] <- 'bold'

plot_df[, event_ordered := factor(event_name_sig, levels = event_order)]

print(plot_df[pval < 0.05])

# print just p values:
print(plot_df[pval < 0.05,c('event_name_sig','pval', 'Event type')])

### Plot filled barplot, ordered by fraction and add n numbers
colour_pal <- unlist(tx_palette$categorical)[1:length(sc_expansion_levels)]
names(colour_pal) <- sc_expansion_levels

g1 <- ggplot(plot_df, aes(event_ordered, n_clones_sc_expansion_pereventid, fill = sc_expansion_factor)) +
	geom_bar(stat = "identity", position = position_fill()) +
	theme_cowplot() +
	scale_fill_manual(values = colour_pal, name = 'Subclone expansion status') +
	geom_text(aes(label = n_clones_sc_expansion_pereventid), position = position_fill(vjust = 0.5)) +
	geom_text(aes(label = pval_signif), position = position_fill(vjust = 0.75), data = plot_df[max_ccf_expandedsubclone == 'Expanded subclone'], size = 8) +
	facet_grid(cols = vars(`Event type`), scales = 'free', space = 'free') +
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
			strip.background = element_rect(fill = "transparent", colour = "black")) +
	labs(x = 'Event ID', y = 'N. subclones with event', title = paste0(histo, " (subclonal events only)"))
# Save into figures/ directory
pdf(paste0('figures/suppfig4a_', histo, '_sc_topevents_prop_sc_expansion_vsbg.pdf'), width = 16, height = 4)
print(g1)
dev.off()
