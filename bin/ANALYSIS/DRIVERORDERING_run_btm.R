#!/usr/bin/env Rscript

# To run this script for both cohorts 'primary' and 'mets', and LUAD_all_topfreqevents.csv and LUSC_all_topfreqevents.csv
# run the following command from the command line (make sure you are in the root directory of this repo): 
"
# LUAD, primary cohort
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUAD_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin

# LUSC, primary cohort
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUSC_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin

# LUAD, mets cohort
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/mets/mets_default/cohort_outputs/eventordering/LUAD_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/mets \
--bin_directory bin

# LUSC, mets cohort
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/mets/mets_default/cohort_outputs/eventordering/LUSC_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/mets \
--bin_directory bin
"

cat("\nRunning driver ordering analysis\n")

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
  make_option(c("--top_events_path"), type="character", action="store", help="File path to combined cohort edge copy number change file with chrarm added info", metavar="character"),
  make_option(c("--cohort_tree_directory"), type="character", help="Directory where tree data is stored"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

top_events_path <- opt$top_events_path
cohort_tree_directory <- opt$cohort_tree_directory
bin_dir <- opt$bin_directory

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

####################################
### Run driver ordering analysis ###
####################################
cohort <- basename(cohort_tree_directory)

### Preprocess data
top_events <- fread(top_events_path)

if (grepl("mets", basename(cohort_tree_directory))) {
  # Mets clone class groups:
  mets_clusters <- fread(file.path(assets_dir, "cloneInfo_df.csv"))
  setnames(mets_clusters, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
  mets_clusters[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]
  seeding_clones <- mets_clusters[seedingClones == TRUE, .(tumour_id, clone, event_id = 'Seeding')]
  
  seeding_events <- unique(top_events[, .(tumour_id, clone, parent_clone)])
  seeding_events <- merge.data.table(seeding_events, seeding_clones, by = c("tumour_id", "clone"), all.x = FALSE, all.y = FALSE)
  top_events <- rbind(top_events, seeding_events)
  setorder(top_events, tumour_id)
}

top_events[, clone := as.character(clone)]
top_events[, Clonality := ifelse(parent_clone == "diploid", "Truncal", "Subclonal")]
top_events[, n_tums_with_event := length(unique(tumour_id)), by = event_id]
top_events[, n_clones_with_event := .N, by = event_id]
top_events[, n_clones_with_event_clonality := .N, by = .(event_id, Clonality)]
top_events[is.na(event_id), n_clones_with_event := NA]
top_events[is.na(event_id), n_clones_with_event_clonality := NA]
top_events[is.na(event_id), n_tums_with_event := NA]

unique_top_events <- sort(unique(top_events$event_id))

### Get ordering of number of events across the cohort
top_events_counts <- unique(top_events[!is.na(event_id), .(event_id, n_clones_with_event, n_tums_with_event)])
setorder(top_events_counts, -n_clones_with_event)
top_events_prevalence_order <- unique(top_events_counts$event_id)

### Get ordering of number of events across the cohort, split by clonal/subclonal
top_events_counts_clonality <- unique(top_events[!is.na(event_id), .(event_id, n_clones_with_event_clonality, n_clones_with_event, Clonality, n_tums_with_event)])

### Get mean/max CCF of each clone in the tumour
top_events_clones <- unique(top_events[, .(tumour_id, clone)])
top_events_clones[, c("mean_ccf_tumour", "max_ccf_tumour", "max_ccf_clonalillusion") := 
  rbindlist(lapply(seq(nrow(top_events_clones)), function(i) {
    # Load tree data
    tum_id <- top_events_clones[i, tumour_id]
    c <- top_events_clones[i, clone]
    tree_object <- readRDS(file = paste0(cohort_tree_directory, '/' , tum_id, '_trees.rds'))
    ccf_cluster_table <- tree_object$nested_pyclone$ccf_cluster_table
    clonality_table <- tree_object$clonality_out$clonality_table_corrected
    # Extract CCF metrics
    clone_ccf <- ccf_cluster_table[rownames(ccf_cluster_table) == c, ]
    clone_clonality <- clonality_table[rownames(clonality_table) == c, ]
    max_ccf_region <- names(clone_ccf)[which(clone_ccf == max(clone_ccf))]
    max_ccf_clonalillusion <- ifelse("clonal" %in% as.character(clone_clonality[max_ccf_region]), 'Clonal illusion', 'Non clonal illusion')
    out_df <- data.table(mean_ccf_tumour = mean(clone_ccf),
                        max_ccf_tumour = max(clone_ccf),
                        max_ccf_clonalillusion = max_ccf_clonalillusion)
    return(out_df)
}))]
# Merge with events:
top_events <- merge.data.table(top_events, top_events_clones, by = c("tumour_id", "clone"), all.x = TRUE)

### Label the driver gene events
top_events[, driver_alterations_clone := paste(event_id, collapse = ";"), by = .(tumour_id, clone)]
top_events[grepl("NA", driver_alterations_clone), driver_alterations_clone := NA]

### Remove tumours with no events:
top_events[, n_clones_w_event_tum := sum(!is.na(driver_alterations_clone)), by = tumour_id]
top_events <- top_events[n_clones_w_event_tum > 0]

### Make an edge matrix in long format for all drivers, every tumour
top_events_tumours <- sort(unique(top_events$tumour_id))
top_events_driver_edge_df <- rbindlist(lapply(top_events_tumours, function(tum_id) {
  tum_events_df <- unique(top_events[tumour_id == tum_id, .(tumour_id, clone, parent_clone, driver_alterations_clone)])
  driv_edges_tum <- get_driver_ancestor_descendant_pairs_tumour(tum_events_df = tum_events_df)
  return(driv_edges_tum)
}))

### Get counts of ancestor -> descendant pairs, removing duplicate parent nodes giving rise to same descendant:
print('Get counts of ancestor -> descendant pairs, removing duplicate parent nodes giving rise to same descendant')
top_events_driver_binary_df <- unique(top_events_driver_edge_df[, .(tumour_id, ancestor = ancestor_event, descendant = descendant_event)])

### Run Bradley Terry model 
top_events_bt_ordering <- run_btm(driver_binary_df = top_events_driver_binary_df, all_unique_alterations = unique_top_events)
# Don't include those genes with a std error > 10
top_events_bt_ordering <- top_events_bt_ordering[`Std. Error` < 10]

####################################
### Plot driver ordering figures ###
####################################
cat('\nPlot driver ordering figures\n')
event_colour_mapping_df <- data.frame(colour = as.character(tx_palette$driver_event_class))
event_colour_mapping_df$`Event type` <- names(tx_palette$driver_event_class)
event_colour_mapping_df <- as.data.table(event_colour_mapping_df)
top_events_bt_ordering[, event_id := gene]
top_events_bt_ordering[grepl("Seeding", event_id), `Event type` := "Seeding"]
top_events_bt_ordering[grepl("del", event_id), `Event type` := "del"]
top_events_bt_ordering[grepl("LOH", event_id), `Event type` := "LOH"]
top_events_bt_ordering[grepl("[+]", event_id), `Event type` := "Amp"]
top_events_bt_ordering[grepl("WGD", event_id), `Event type` := "WGD"]
top_events_bt_ordering[is.na(`Event type`), `Event type` := "SNV"]
top_events_bt_ordering <- merge.data.table(top_events_bt_ordering, event_colour_mapping_df, by = "Event type", all.x = TRUE, all.y = FALSE)
top_events_colour_pal <- top_events_bt_ordering$colour
setorder(top_events_bt_ordering, Estimate)

# make gene text italic
top_events_fontface <- ifelse(top_events_bt_ordering$`Event type` == "SNV", "italic", "plain")

# merge event prevalence
top_events_density <- merge.data.table(top_events_bt_ordering, top_events_counts_clonality, by = 'event_id', all.x = TRUE, all.y = FALSE)
top_events_estimate <- merge.data.table(top_events, top_events_bt_ordering, by = 'event_id', all = FALSE)

### Plot: mean/max CCFs, BTm estimate, event prevalences
top_events_meanccf <- ggplot(top_events_estimate, aes(y = reorder(event_id, Estimate), x = mean_ccf_tumour))+
  geom_density_ridges(aes(fill = `Event type`), alpha = .85) +
  scale_fill_manual(values = tx_palette$driver_event_class, name = "Event type") +
  theme_cowplot() +
  labs(y = paste0(histo, " most frequent events (", clonality, ")"), x = "Mean CCF") +
  xlim(0, 100) +
  theme(axis.text.y = element_text(colour = top_events_colour_pal, face = top_events_fontface),
        plot.margin = unit(c(1, 0, 0, 0.5), "cm"),
        legend.position = "none"
  )

top_events_maxccf <- ggplot(top_events_estimate, aes(y = reorder(event_id, Estimate), x = max_ccf_tumour))+
  geom_density_ridges(aes(fill = `Event type`), alpha = .85) +
  scale_fill_manual(values = tx_palette$driver_event_class, name = "Event type") +
  theme_cowplot() +
  labs(y = NULL, x = "Max CCF") +
  xlim(0, 100) +
  theme(axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = unit(c(1, 0, 0, 0.5), "cm"),
      legend.position = "none"
  )

top_events_est <- ggplot(top_events_bt_ordering,
  aes(x = reorder(event_id, Estimate),
      y = Estimate,
      ymin = (Estimate - `Std. Error`),
      ymax = (Estimate + `Std. Error`),
      colour = `Event type`)
  ) +
  geom_pointrange(aes(x = reorder(event_id, Estimate),
                      ymin = (Estimate - `Std. Error`),
                      ymax = (Estimate + `Std. Error`),
                      y = Estimate)
  ) +
  scale_colour_manual(values = tx_palette$driver_event_class, name = "Event type") +
  theme_cowplot() +
  labs(y = "Relative rank", x = NULL) +
  coord_flip() +
  scale_y_reverse() +
  theme(strip.background = element_rect(fill = "white")) +
  theme(axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = unit(c(1, 0, 0, 0.5), "cm"),
      legend.position = "none"
  )

top_events_prev <- ggplot(top_events_density, aes(reorder(event_id, Estimate), n_clones_with_event_clonality)) +
  geom_col(aes(fill = Clonality)) +
  labs(x = NULL, y = "Event prevalence") +
  theme_cowplot() +
  coord_flip() +
  scale_fill_manual(values = tx_palette$clonality) +
  scale_y_reverse() +
  theme(axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), 
      axis.line.y = element_blank(),
      plot.margin = unit(c(1, 0, 0, 0.5), "cm"))

g <- ggarrange(
  top_events_meanccf,
  top_events_maxccf,
  top_events_est,
  top_events_prev,
  nrow = 1,
  widths = c(2, .8, 1, 2)
)
# Combine plots and save
bt_outpath <- gsub('.csv', '_btm.pdf', top_events_path)
ggsave(filename = bt_outpath, plot = g, width = 10, height = 9, limitsize = FALSE)

# Additionally save into figures/ directory
bt_outpath_fig <- basename(bt_outpath)
fig_prefix <- ifelse(cohort == "primary", "Fig3", "Suppfig4")
ggsave(filename = paste0('figures/', fig_prefix, '_', bt_outpath_fig), plot = g, width = 10, height = 9, limitsize = FALSE)

# Plot legend:
events_legend <- ggplot(top_events_estimate, aes(y = reorder(event_id, Estimate), x = max_ccf_tumour, fill = `Event type`))+
  scale_fill_manual(values = tx_palette$driver_event_class, name = "Event type") +
  geom_density_ridges(aes(fill = `Event type`)) +
  theme_cowplot() 
legend <- cowplot::get_legend(events_legend)
ggsave(filename = gsub(".pdf", ".legend.pdf", bt_outpath), plot = legend, width = 3, height = 3, limitsize = FALSE)




### Plot pairwise event matrix

### Show count + fraction of co-occurrences that we are able to time using ancestral-descendant nodes
driver_binary_df <- copy(top_events_driver_binary_df)
top_events_pairs_df <- get_pairwise_event_freq(driver_binary_df)

# get only lower triangle:
lowertri_pairs <- get_lower_tri_pairwise(top_events_pairs_df, top_events_prevalence_order)

## Plot:
top_events_pairs_df[, `Ancestor event` := factor(as.character(ancestor), levels = top_events_prevalence_order)]
top_events_pairs_df[, `Descendant event` := factor(as.character(descendant), levels = top_events_prevalence_order)]
top_events_pairs_df[, anc_dec := paste(c(as.character(ancestor), as.character(descendant)), collapse = ';'), by = seq(nrow(top_events_pairs_df))]
top_events_pairs_df[, n_cooccurrences_plot := n_cooccurrences]
top_events_pairs_df[n_cooccurrences == 0 | !(anc_dec %in% lowertri_pairs$pair), n_cooccurrences_plot := NA]
top_events_pairs_df[, same_pair := ancestor == descendant]

g <- ggplot(top_events_pairs_df, aes(`Descendant event`, forcats::fct_rev(`Ancestor event`))) +
  geom_point(aes(size = n_cooccurrences_plot, color = frac_anc_dec), shape = 15, stat = "identity") +
  geom_point(aes(size = n_cooccurrences_plot), shape = 0, color = "#374E55FF") +
  # scale_size_continuous(range = c(0, 6), name = "N. tumours with\ntemporal ordering\npossible") +
  scale_size_continuous(range = c(0, 6), name = "N. tumours with temporal ordering possible") +
  theme_cowplot() +
  scale_colour_gradient2(low = "blue",
                          mid = 'white',
                          high = "red2",
                          midpoint = .5,
                          breaks = c(0, 0.5, 1),
                          labels = c(0, 0.5, 1),
                          name = "Frac. tumours with ancestor -> descendant"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 14, hjust = 0),
        panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "#f0f0f0"),
        legend.position = "bottom",
        legend.box = "vertical", 
        legend.margin = margin(),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18),
        axis.text.x.top = element_text(vjust = 0.5)
        ) +
  scale_x_discrete(position = "top") +
  # labs(title = paste0(histo, " most frequent events (", clonality, ")"), y = "Ancestor event", x = "Descendant event")
  labs(title = paste0(histo, " most frequent events (", clonality, ")"), y = "Ancestor event", x = "Descendant event")

pairwise_outpath <- gsub('.csv', '_pairwise.pdf', top_events_path)
ggsave(filename = pairwise_outpath, plot = g, width = 9.5, height = 10, limitsize = FALSE)

# Additionally save into figures/ directory
pairwise_outpath_fig <- basename(pairwise_outpath)
ggsave(filename = paste0('figures/', fig_prefix, '_', pairwise_outpath_fig), plot = g, width = 9.5, height = 10, limitsize = FALSE)

###########
### END ###
###########