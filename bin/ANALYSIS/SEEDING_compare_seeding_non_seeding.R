#!/usr/bin/env Rscript
# To run this script run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/SEEDING_compare_seeding_non_seeding.R \
--clone_distances_path output/mets/mets_default/cohort_outputs/edge_events.csv \
--clone_metrics_path output/mets/mets_default/cohort_outputs/clone_ploidy_floh.csv \
--bin_directory bin \
--save_directory output/mets/mets_default/cohort_outputs/n_events
"
cat("\n Comparing #SCNAs between different clone groups in mets cohort\n")

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(scales))


### Command line options
option_list <- list(
  make_option(c("--clone_distances_path"), type="character", action="store", help="File path to combined cohort clone genomic distances file", metavar="character"),
  make_option(c("--clone_metrics_path"), type="character", action="store", help="File path to combined cohort clone genomic distances file", metavar="character"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--save_directory"), type="character", help="save_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
clone_distances_path <- opt$clone_distances_path
clone_metrics_path <- opt$clone_metrics_path
bin_dir <- opt$bin_directory
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

assets_dir <- gsub('bin', '_assets', bin_dir)

### Functions
source(file.path(bin_dir, "/ANALYSIS/analysis_functions.R"))

### Load and process data
# Clinical data
clin_data_path <- file.path(assets_dir, "clinical_data_tumour.tsv")
clin_small <- load_reduced_clinicaldata(clin_data_path)

# Clone genomic distances data
distances_df <- fread(clone_distances_path)

# Clone ploidy + floh 
clone_ploidy_floh <- fread(clone_metrics_path)

# Mets clone class groups:
mets_cloneInfo <- fread(file.path(assets_dir, "cloneInfo_df.csv"))
setnames(mets_cloneInfo, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
mets_cloneInfo[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]

# Relapse category
patient_df <- fread(file.path(assets_dir, "clinical_data.tsv"))
relapse_df <- unique(patient_df[, .(patient_id = `Patient ID`, `Relapse category`)])

# Seeding info
seeding_df <- fread(file.path(assets_dir, "seedingTableFull.csv"))

# Colour palette
tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))

# evo metrics
evo <- fread(file.path(assets_dir, "20221110_TRACERx421_evolutionary_metrics.tsv"))
scna_ith <- evo[tumour_id %in% unique(distances_df$tumour_id), .(tumour_id, frac_abberant_genom_subcl)]


boxplot_theme <- theme(strip.background = element_rect(fill = "transparent", colour = "black"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))

##############################
### Process and merge data ###
##############################

edge_counts_df <- process_edge_snv_scna_dist(edge_distance_df = distances_df)
setnames(
  edge_counts_df,
  c('total_interval_gains_bin', 'total_interval_losses_bin', 'total_interval_events_bin'),
  c('total_no_gains', 'total_no_losses', 'total_no_events')
)

# Merge with clinical info:
edge_counts_df[tumour_id == "CRUK0084-Cluster2", tumour_id := "CRUK0084"] 
edge_df <- merge.data.table(edge_counts_df, clin_small, by = "tumour_id", all.x = TRUE, all.y = FALSE)
edge_df[tumour_id == "CRUK0084", tumour_id := "CRUK0084-Cluster2"] 

# Merge with cluster annotation:
mets_cloneInfo[, clone := as.character(clone)]
edge_df[, clone := as.character(clone)]
mets_dist <- merge.data.table(edge_df, mets_cloneInfo, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)

mets_dist <- assign_clone_classes(mets_dist)

mets_dist[, mrca_seeding_tum := any(parent_clone == "diploid" & seedingClones == TRUE), by = tumour_id]
mets_dist[, mrca_seeding_tum_fac := ifelse(mrca_seeding_tum, "Truncal seeding", "Subclonal seeding")]
mets_dist[, mrca_seeding_class := ifelse(clone_class == 'Non-seeding MRCA', 'Non-seeding MRCA',
                                    ifelse(parent_clone == "diploid" & seedingClones == TRUE, 'Seeding MRCA',
                                      ifelse(is_mrca == "subclone", "subclone", NA)))]

# Merge with seeding clonality:
mets_dist <- merge.data.table(mets_dist, seeding_df, by = "tumour_id", all.x = TRUE, all.y = FALSE)

# save
mets_dist[seeding_class == "Non-seeding", seeding_vs_nonseedingprim := "Non-seeding primary"]
mets_dist[clone_class == "Seeding" & mrca_seeding_tum == FALSE, seeding_vs_nonseedingprim := clone_class]
seeding_class_df <- unique(mets_dist[!is.na(seeding_vs_nonseedingprim), .(tumour_id, clone, seeding_vs_nonseedingprim)])
# write.csv(seeding_class_df, '_assets/seeding_vs_nonseeding_prim.csv', row.names = FALSE, quote = FALSE)
# mets_dist[is.na(seeding_vs_nonseedingprim) & is_mrca == 'subclone' & clone_class == 'Seeding']

# Split categories: mrca, subclonal, metastatic
mets_dist[, met_categories := ifelse(parent_clone == "diploid", 'MRCA',
                                ifelse(clone_class == "Metastatic", "Metastatic\nclones", "Subclonal primary\nclones"))]
mets_dist[, met_categories := factor(met_categories, levels = c('MRCA', "Subclonal primary\nclones", "Metastatic\nclones"))]


# Merge with relapse:
mets_dist[, patient_id := gsub('-Cluster.*', '', tumour_id)]
mets_dist <- merge.data.table(mets_dist, relapse_df, by = "patient_id", all.x = TRUE, all.y = FALSE)

# Add tree level
mets_dist[, tree_level := sapply(cluster_ancestors, function(x){
    ancscs <- strsplit(x, ";")[[1]]
    if (length(ancscs) == 1) return(0)
    else {
        return(length(ancscs) - 1)
    }
})]
mets_dist[, tree_level_fac := as.character(tree_level)]
mets_dist[tree_level >= 6, tree_level_fac := '>5']
mets_dist[, tree_level_fac := factor(tree_level_fac, levels = c(as.character(1:5), '>5'))]

# n. scnas in primary tumour
mets_dist[, sum_cn_events_primary_tumour := sum(total_no_events * (metClones == FALSE)), by = .(tumour_id)]
mets_dist[, total_n_cn_events_per_clone := sum_cn_events_tumour / n_clones]

mets_dist[, sum_cn_events_primary_tumour_sc := sum(total_no_events * ((metClones == FALSE) & (parent_clone != "diploid"))), by = .(tumour_id)]
mets_dist[, n_prim_subclones := sum((metClones == FALSE) & (parent_clone != "diploid")), by = .(tumour_id)]
mets_dist[, total_n_sc_cn_events_per_prim_sc := sum_cn_events_primary_tumour_sc / n_prim_subclones]
mets_dist[, total_n_cn_events_per_prim_clone := sum_cn_events_primary_tumour / (n_prim_subclones + 1)]

# # Get cumulative # events
cumu_df <- get_cumulative_snv_scna_events(mets_dist)


############
### MAIN ###
############

####################################################################################################
### Are there more SCNAs in the seeding clone compared to non-seeding clones in the same tumour? ###
####################################################################################################
cat('\nAre there more SCNAs in the seeding clone compared to non-seeding clones in the same tumour?\n')
metscluster_n_events <- melt.data.table(mets_dist, measure.vars = c('total_no_losses', 'total_no_gains', 'total_no_events'))
metscluster_n_events[, metric := gsub("total_no_", "N. ", variable)]

colour_pal <- unlist(tx_palette$clone_class)[c("Primary", "Seeding")]
names(colour_pal) <- c("Non-seeding primary", 'Seeding')

### Figure 4e. 
# Comparing #SCNAs in non-seeding clones vs seeding clones 
g <- ggplot(metscluster_n_events[mrca_seeding_tum == FALSE & !is.na(seeding_vs_nonseedingprim)],
  aes(seeding_vs_nonseedingprim, value, fill = seeding_vs_nonseedingprim)) +
  facet_wrap(~metric)+
  scale_fill_manual(values = colour_pal, name = "Clone class") +
  theme_cowplot() +
  geom_jitter(size = .5, width = .25, height = 0, colour = "grey50")+
  geom_boxplot(alpha = .9, outlier.shape = NA) +
  boxplot_theme +
  stat_compare_means(size = 4)+
  labs(x = "Clone class", y = "#SCNAs")#, title = "Tumours with subclonal seeding only")
ggsave(
  filename = paste0(save_directory, "/seeding_vs_nonseeding_scnaevents_mrcassedintumFALSE.pdf"),
  plot = g,
  width = 6,
  height = 4
)
# Also save in figures directory
ggsave(
  filename = "figures/Fig4e_seeding_vs_nonseeding_scnaevents_mrcassedintumFALSE.pdf",
  plot = g,
  width = 6,
  height = 4
)

### Supplementary Figure 5b
# Comparing distribution of seeding clone branch length vs distribution of non-seeding clone branch lengths

comps <- as.list(data.frame(combn(as.character(unique(metscluster_n_events$clone_class)), 2)))
g <- ggplot(metscluster_n_events[mrca_seeding_tum == FALSE],
  aes(clone_class, value, fill = clone_class)) +
  facet_wrap(~metric + mrca_seeding_tum_fac, scales = "free")+
  geom_jitter(size = .5, height = 0, colour = "grey50")+
  geom_boxplot(alpha = .9, outlier.shape = NA) +
  scale_fill_manual(values = tx_palette$clone_class, name = "Clone class") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "transparent", colour = "black"),
        # axis.text.x = element_text(angle = 45, hjust = 1)
        axis.text.x = element_blank()
        ) +
  stat_compare_means(comparisons = comps, label = "p.signif")+
  labs(x = "Clone class", y = "No. copy number events")
ggsave(
  filename = paste0(save_directory, "/cloneclass_scnaevents_mrcassedintumFALSE.pdf"),
  plot = g,
  width = 8.7,
  height = 4.5
)
# Also save in figures directory
ggsave(
  filename = "figures/Suppfig5b_cloneclass_scnaevents_mrcassedintumFALSE.pdf",
  plot = g,
  width = 8.7,
  height = 4.5
)

### Supplementary Figure 5d 
# Comparing SCNA/SNV slopes in seeding vs non-seeding tumours:
g <- ggplot(mets_dist[mrca_seeding_tum == FALSE & !is.na(seeding_vs_nonseedingprim)],
  aes(seeding_vs_nonseedingprim, slope_raw)) +
  geom_jitter(size = .5, width = .25, height = 0, colour = "grey50")+
  geom_boxplot(alpha = .9, aes(fill = clone_class), outlier.shape = NA) +
  scale_fill_manual(values = tx_palette$clone_class, name = "Clone class") +
  theme_cowplot() +
  boxplot_theme +
  stat_compare_means(size = 4)+
  labs(x = "Clone class", y = "#SCNAs / #SNVs") +
  guides(fill = guide_legend(nrow = 2))
ggsave(
  filename = paste0(save_directory, "/seeding_vs_nonseeding_slope_raw_mrcassedintumFALSE.pdf"),
  plot = g,
  width = 3.,
  height = 4
)
ggsave(
  filename = "figures/Suppfig5d_seeding_vs_nonseeding_slope_raw_mrcassedintumFALSE.pdf",
  plot = g,
  width = 3.,
  height = 4
)


#########################################################################
### Are there more SCNAs in polyclonal vs monoclonal seeding tumours? ###
#########################################################################
cat('\nAre there more SCNAs in polyclonal vs monoclonal seeding tumours?\n')
### Figure 4k
# #SCNA events in poly vs monoclonal seeding tumours, by tree level
g <- ggplot(mets_dist,
  aes(tree_level_fac, total_no_events, fill = seedingClonality)) +
  # geom_jitter(aes(group = seedingClonality), width = 0, size = .5, height = 0, colour = 'grey50') +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0), size = .5, colour = 'grey50') +
  geom_boxplot(outlier.shape = NA, alpha = .9) +
  scale_fill_manual(values = unlist(tx_palette$seedingClonality), name = 'Tumour seeding clonality') +
  theme_cowplot() +
  stat_compare_means(aes(group = seedingClonality), label = "p.signif", hide.ns = TRUE, colour = "red", size = 6, tip.length = .01, label.y = 400) +
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)) +
  # guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Tree level", y = "#SCNAs")
ggsave(
  filename = paste0(save_directory, "/seedingClonality_scnaevents_treelevel_boxplot.pdf"),
  plot = g,
  width = 6.5,
  height = 4
)
ggsave(
  filename = "figures/Fig4k_seedingClonality_scnaevents_treelevel_boxplot.pdf",
  plot = g,
  width = 6.5,
  height = 4
)

### CUMULATIVE #cn events in seeding clonalities
g <- ggplot(cumu_df,
  aes(tree_level_fac, cum_total_no_events, fill = seedingClonality)) +
  geom_boxplot() +
  stat_compare_means(aes(group = seedingClonality), label = "p.signif", hide.ns = TRUE, colour = "red", size = 3, tip.length = .01) +
  scale_fill_manual(values = unlist(tx_palette$seedingClonality)) +
  scale_colour_manual(values = unlist(tx_palette$seedingClonality)) +
  theme_cowplot() +
  labs(x = "Tree level", y = "Cumulative no. copy number events")
ggsave(
  filename = paste0(save_directory, "/seedingClonality_cum_total_no_events_treelevel_boxplot.pdf"),
  plot = g,
  width = 6,
  height = 4
)

###############################################################################
### Are there more SCNAs in extrathoracic vs intrathoracic seeding tumours? ###
###############################################################################
cat('\nAre there more SCNAs in extrathoracic vs intrathoracic seeding tumours?\n')

mets_dist_relapse <- mets_dist[`Relapse category` %in% c("Intrathoracic", "Extrathoracic")]
comps <- as.list(data.frame(combn(unique(mets_dist_relapse$`Relapse category`), 2)))

# Supp figure 5h
metscluster_n_events <- melt.data.table(mets_dist, measure.vars = c('total_no_losses', 'total_no_gains'))
metscluster_n_events[, metric := gsub("total_no_", "N. ", variable)]
metscluster_n_events_relapse <- metscluster_n_events[`Relapse category` %in% c("Intrathoracic", "Extrathoracic")]
g <- ggplot(metscluster_n_events_relapse,
    aes(`Relapse category`, value)) +
    facet_wrap(~ met_categories + metric, scales = "free_y", nrow = 1)+
    geom_jitter(size = .5, height = 0, width = .25, colour = "grey50")+
    geom_boxplot(alpha = .9, outlier.shape = NA, aes(fill = `Relapse category`)) +
    theme_cowplot() +
    scale_fill_manual(values = unlist(tx_palette$relapseCat), name = "Tumour relapse category") +
    boxplot_theme +
    stat_compare_means(size = 4, label.y.npc = .8)+
    labs(x = "Tumour relapse category", y = "#SCNAs")
ggsave(
  filename = paste0(save_directory, "/relapse_category_scnaevents_by_met_categories.pdf"),
  plot = g,
  width = 12,
  height = 4.5
)
ggsave(
  filename = "figures/Suppfig5h_relapse_category_scnaevents_by_met_categories.pdf",
  plot = g,
  width = 12,
  height = 4.5
)


### Figure 4l
# Only plot #SCNAs
g <- ggplot(mets_dist_relapse,
    aes(`Relapse category`, total_no_events)) +
    facet_wrap(~met_categories, nrow = 1)+
    geom_jitter(size = .5, width = .25, height = 0, colour = "grey50")+
    geom_boxplot(alpha = .9, outlier.shape = NA, aes(fill = `Relapse category`)) +
    theme_cowplot() +
    scale_fill_manual(values = unlist(tx_palette$relapseCat), name = "Tumour relapse category") +
    boxplot_theme +
    stat_compare_means(size = 4) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "Tumour relapse category", y = "#SCNAs")
ggsave(
  filename = paste0(save_directory, "/relapse_category_SCNABL_by_met_categories.pdf"),
  plot = g,
  width = 6,
  height = 4
)
ggsave(
  filename = "figures/Fig4l_relapse_category_scnaevents_by_met_categories.pdf",
  plot = g,
  width = 12,
  height = 4.5
)


# Test including intra & extra as extra:
any_relapse <- mets_dist[`Relapse category` %in% c("Intrathoracic", "Intra & Extra", "Extrathoracic")]
any_relapse[`Relapse category` == "Intra & Extra", `Relapse category` := "Extrathoracic"]
g <- ggplot(any_relapse,
    aes(`Relapse category`, total_no_events)) +
    facet_wrap(~met_categories, nrow = 1)+
    geom_jitter(size = .5, width = .25, height = 0, colour = "grey50")+
    geom_boxplot(alpha = .9, outlier.shape = NA, aes(fill = `Relapse category`)) +
    theme_cowplot() +
    scale_fill_manual(values = unlist(tx_palette$relapseCat), name = "Tumour relapse category") +
    boxplot_theme +
    stat_compare_means(size = 4) +
    guides(fill = guide_legend(nrow = 2)) +
    labs(x = "Tumour relapse category", y = "#SCNAs")
ggsave(
  filename = paste0(save_directory, "/relapse_category_anyextra_SCNABL_by_met_categories.pdf"),
  plot = g,
  width = 6,
  height = 4
)

###########
### END ###
###########
