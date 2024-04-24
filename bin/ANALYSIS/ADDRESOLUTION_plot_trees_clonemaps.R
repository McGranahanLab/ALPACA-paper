#!/usr/bin/env Rscript
"
# Example 1
Rscript bin/ANALYSIS/ADDRESOLUTION_plot_trees_clonemaps.R \
--tumour_id CRUK0048 \
--save_directory output/mets/mets_default/cohort_outputs/example_cases

# Example 2
Rscript bin/ANALYSIS/ADDRESOLUTION_plot_trees_clonemaps.R \
--tumour_id CRUK0022 \
--save_directory output/mets/mets_default/cohort_outputs/example_cases
"

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(igraph))
# suppressPackageStartupMessages(library(cloneMap))

### Command line options 
option_list <- list(
  make_option(c("--tumour_id"), type="character", action="store", metavar="character"),
  make_option(c("--save_directory"), type="character", help="save_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
tum_id <- opt$tumour_id
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

### Functions
source("bin/ANALYSIS/analysis_functions.R")
source("bin/ANALYSIS/plotting_functions.R")

### Mets clone class groups:
mets_cloneInfo <- fread("_assets/cloneInfo_df.csv")
setnames(mets_cloneInfo, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
mets_cloneInfo[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]

############
### MAIN ###
############

# # Example 1
# tum_id <- "CRUK0048"

# # # Example 2
# tum_id <- "CRUK0022"


### PLOT ###
tum_out_dir <- file.path(save_directory, tum_id)
tree_paths <- read_json(path = file.path(tum_out_dir, "tree_paths.json"))
cp_table <- fread(file.path(tum_out_dir, "cp_table.csv"))
tum_mets_cloneInfo <- mets_cloneInfo[tumour_id == tum_id]

# Plot tree
plot_seeding_tree(tree_paths, tum_mets_cloneInfo, tum_out_dir)

# # Plot clonemaps
# colour_vec <- colour_seeding_tree(tree_paths, tum_mets_cloneInfo)
# plot_region_clonemaps(tree_paths, cp_table, tum_out_dir, colour_vec)

###########
### END ###
###########