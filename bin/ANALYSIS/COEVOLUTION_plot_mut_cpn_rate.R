#!/usr/bin/env Rscript

# To run this script run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/COEVOLUTION_plot_mut_cpn_rate.R \
--clone_distances_path output/mets/mets_default/cohort_outputs/edge_events.csv \
--bin_directory bin \
--save_directory output/mets/mets_default/cohort_outputs/constellation
"
cat("\nPlotting SNV vs SCNA tree plots\n")

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(tidytext))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))

### Command line options
option_list <- list(
  make_option(c("--clone_distances_path"), type="character", action="store", help="File path to combined cohort clone genomic distances file", metavar="character"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--save_directory"), type="character", help="Plot output directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
clone_distances_path <- opt$clone_distances_path
bin_dir <- opt$bin_directory
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

assets_dir <- gsub('bin', '_assets', bin_dir)

### Functions
source(file.path(bin_dir, "/ANALYSIS/analysis_functions.R"))
source(file.path(bin_dir, "/ANALYSIS/plotting_functions.R"))

### Load and process data
# Clinical data
clin_data_path <- file.path(assets_dir, "clinical_data_tumour.tsv")
clin_small <- load_reduced_clinicaldata(clin_data_path)
clin_small[, stage_short := gsub("A|B", "", pathologyTNM)]
clin_small[, stage_short := factor(stage_short, levels = c("I", "II", "III"))]

# Clone genomic distances data
distances_df <- fread(clone_distances_path)

# Mets clone class groups:
mets_cloneInfo <- fread(file.path(assets_dir, "cloneInfo_df.csv")) # change to params.cloneinfo e.g.
setnames(mets_cloneInfo, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
mets_cloneInfo[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]

# Relapse category
patient_df <- fread(file.path(assets_dir, "clinical_data.tsv"))
relapse_df <- unique(patient_df[, .(patient_id = `Patient ID`, `Relapse category`)])

# Dissemination pattern info
seeding_df <- fread(file.path(assets_dir, "seedingTableFull.csv"))

# Colour palette
tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))

########################################################
### Compute cumulative event count on each tree path ###
########################################################

# Process edge event data:
distance_total_df <- process_edge_snv_scna_dist(edge_distance_df = distances_df)
setnames(distance_total_df,
          c('total_interval_gains_bin', 'total_interval_losses_bin', 'total_interval_events_bin'),
          c('total_no_gains', 'total_no_losses', 'total_no_events'))

# Get tree paths from edge data frame:
alpaca_edge_df <- unique(distance_total_df[, .(tumour_id, clone, parent_clone)])
tree_paths_df <- generate_tree_paths_df(alpaca_edge_df, verbose = FALSE)

# Get cumulative events
distance_total_df <- get_cumulative_snv_scna_events(distance_total_df)

# Merge with clinical
distance_total_df[tumour_id == "CRUK0084-Cluster2", tumour_id := "CRUK0084"] 
distance_total_df <- merge.data.table(distance_total_df, clin_small, by = "tumour_id", all.x = TRUE, all.y = FALSE)
distance_total_df[tumour_id == "CRUK0084", tumour_id := "CRUK0084-Cluster2"] 

# Add clone class assignment
mets_cloneInfo[, clone := as.character(clone)]
distance_total_df <- merge.data.table(distance_total_df, mets_cloneInfo, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)
distance_total_df <- assign_clone_classes(edge_clone_info_df = distance_total_df)

# Merge with seeding clonality:
distance_total_df <- merge.data.table(distance_total_df, seeding_df, by = "tumour_id", all.x = TRUE, all.y = FALSE)

# Merge with relapse:
distance_total_df[, patient_id := gsub('-Cluster.*', '', tumour_id)]
distance_total_df <- merge.data.table(distance_total_df, relapse_df, by = "patient_id", all.x = TRUE, all.y = FALSE)


###############################################################################
### Compute whether rate of mets CN higher or primary CN higher in a tumour ###
###############################################################################
### Order the tumours by the difference in average mets path slope to average primary path slope

met_prim_path_slope_diff <- get_met_vs_prim_path_slope_diff(distance_total_df, tree_paths_df)

#####################################################
### Compute the change in gradients at each point ###
#####################################################

parent_slopes <- distance_total_df[, .(tumour_id, clone, parent_slope = slope_raw)]
setnames(parent_slopes, "clone", "parent_clone")
dist_grad_df <- merge.data.table(distance_total_df, parent_slopes, by = c("tumour_id", "parent_clone"), all.x = TRUE, all.y = FALSE)
dist_grad_df[, slope_change := parent_slope - slope_raw]

dist_grad_df[, mean_slope_change_tum := mean(slope_change, na.rm = TRUE), by = "tumour_id"]
# dist_grad_df[, .(tumour_id, clone, parent_clone, slope_raw, parent_slope, slope_change, mean_slope_change_tum)]


########################################################
### Plot mets cohort coloured by seeding clone class ###
########################################################

### Specify plotting parameters 
alpaca_tumours <- unique(distance_total_df$tumour_id)

# max_n_cn_events <- distance_total_df[, max(cum_total_no_events)]
# max_n_gains <- distance_total_df[, max(cum_total_no_gains)]
# max_n_losses <- distance_total_df[, max(cum_total_no_losses)]
# max_n_muts <- distance_total_df[, max(cum_mut_branchlength)]
# max_n_cn_events <- 400
# max_n_muts <- 1000

# Plot formatting 
nrows <- 13
plot_height <- 22
plot_width <- 3
plot_titlesize <- 45

colour_pal <- unlist(tx_palette$clone_class)

nluad <- length(unique(distance_total_df[Histology == 'LUAD', tumour_id]))

#########################
### MUTS VS CN EVENTS ###
#########################
# Plot formatting 
nrows <- 9
plot_height <- 22
plot_width <- 12
plot_titlesize <- 45

# ordered by n clones 
setorder(met_prim_path_slope_diff, n_clones)

for (histo in unique(met_prim_path_slope_diff$Histology)) {
  cat('\n', histo, '\n')
  histo_tums <- unique(met_prim_path_slope_diff[Histology == histo, tumour_id])
  nhisto <- length(histo_tums)
  # Generate plots
  evo_rate_plots <- lapply(histo_tums, function(tum_id) {
    g <- generate_mut_cn_rate_tree(tum_distance_df = distance_total_df[tumour_id == tum_id],
                                  tum_tree_paths_df = tree_paths_df[tumour_id == tum_id],
                                  colour_map = distance_total_df[tumour_id == tum_id, .(clone, colour_label = clone_class)],
                                  colour_pal = colour_pal,
                                  patient_id = gsub("-Cluster.*", "", tum_id))
    return(g)
  })
  # Save plots
  # outpath <- paste0(save_directory, '/', histo, "_constellation_cloneclass.pdf")
  outpath <- paste0('figures/Fig4b_', histo, "_constellation_cloneclass.pdf")
  pdf(outpath, width = plot_width * (nhisto/nluad), height = plot_height)
  print(marrangeGrob(evo_rate_plots, ncol = ceiling(nhisto/nrows), nrow = nrows, top = textGrob(histo, gp = gpar(fontsize = plot_titlesize, font = 8))))
  dev.off()
  # Save plots again
  pdf(outpath, width = plot_width * (nhisto/nluad), height = plot_height)
  print(marrangeGrob(evo_rate_plots, ncol = ceiling(nhisto/nrows), nrow = nrows, top = textGrob(histo, gp = gpar(fontsize = plot_titlesize, font = 8))))
  dev.off()
}

####################################
### Plot trees for supplementary ###
####################################

tums_to_print <- c('CRUK0083', 'CRUK0092', 'CRUK0742', 'CRUK0035')

for (tum_id in tums_to_print) {
  print(tum_id)
  tum_tree <- as.matrix(unique(distance_total_df[tumour_id == tum_id, .(parent_clone, clone)]))
  colour_map <- distance_total_df[tumour_id == tum_id, .(clone, colour_label = clone_class)]
  colour_map[, colour := sapply(colour_label, function(l) colour_pal[as.character(l)])]
  colour_vec <- colour_map[, colour]
  names(colour_vec) <- colour_map[, clone]
  patient_id <- gsub("-Cluster.*", "", tum_id)
  # pdf(paste0(save_directory, '/', tum_id, "treeplot.pdf"), width = 4, height = 4)
  pdf(paste0('figures/Suppfig6e', tum_id, "treeplot.pdf"), width = 4, height = 4)
  plot_tree(tree_edges = tum_tree,
          # vertex_size = vertex_size_vec,
          vertex_size = 25,
          vertex_label_cex = 1.2,
          colour_vector = colour_vec,
          label_tree = FALSE,
          plot_title = tum_id, 
          title_size = 25)
  dev.off()
}


