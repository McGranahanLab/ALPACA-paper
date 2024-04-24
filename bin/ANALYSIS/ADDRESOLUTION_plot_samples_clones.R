#!/usr/bin/env Rscript
"
# Example 1
Rscript bin/ANALYSIS/ADDRESOLUTION_plot_samples_clones.R \
--tumour_id CRUK0048 \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/example_cases

# Example 2
Rscript bin/ANALYSIS/ADDRESOLUTION_plot_samples_clones.R \
--tumour_id CRUK0022 \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/example_cases
"

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(optparse))


### Command line options 
option_list <- list(
  make_option(c("--tumour_id"), type="character", action="store", metavar="character"),
  make_option(c("--output_directory"), type="character", help="output_directory"),
  make_option(c("--save_directory"), type="character", help="save_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
tum_id <- opt$tumour_id
output_directory <- opt$output_directory
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

### Functions
source('bin/ANALYSIS/analysis_functions.R')
source("bin/ANALYSIS/plotting_functions.R")


############
### MAIN ###
############
### EXAMPLE TUMOUR CASES ###
# - Run the following chunks plus the code below separately to produce output for the two example cases
# - Note: The SNP dataframes (asas_df.csv) comprising the genomic location and fractional copy number of each SNP 
# for the example tumour cases is confidential patient information and thus not provided in the code release for this manuscript

if (tum_id == "CRUK0048") {
  
  ### BRCA2 EXAMPLE ###
	cohort <- "mets"
	run_name <- "mets_default"
	# asas_path <- "../publication/output/mets/mets_default/cohort_outputs/example_cases/CRUK0048/asas_df.csv"
  asas_path <- ""
	gene_list <- c("BRCA2", "RB1")
	chromosome <- 13
	min_startpos <- 0
	max_endpos <- 1.4e+08
	sample_levels <- c(`U_CRUK0048_SU_T1.R1` = "Primary R1", `U_CRUK0048_SU_T1.R2` = "Primary R2", `U_CRUK0048_SU_T1.R3` = "Primary R3", `U_CRUK0048_BR_T1.R1` = "Recurrence 1", `U_CRUK0048_BR_T1.R2` = "Recurrence 2")
	clone_levels <- c(13, 4, 8, 1, 2, 9, 3, 12, 11, 6, 10)
	xaxis_breaks <- c(3e+07, 6e+07, 9e+07)
	plot_width <- 3.8
	plot_height <- 7
} else if (tum_id == "CRUK0022") {
  
  ### MDM2 EXAMPLE ###
  cohort <- "mets"
  run_name <- "mets_default"
  # asas_path <- "../publication/output/mets/mets_default/cohort_outputs/example_cases/CRUK0022/asas_df.csv"
  asas_path <- ""
  gene_list <- c("MDM2", "CDK4")
  chromosome <- 12
  min_startpos <- 3e+07
  max_endpos <- 1.4e+08
  sample_levels <- c(`M_CRUK0022_SU_T1.R1` = "Primary R1", `M_CRUK0022_SU_T1.R2` = "Primary R2", `M_CRUK0022_SU_FLN1` = "Lymphnode 1")
  clone_levels <- c(10, 9, 5, 11, 15, 7, 3, 4, 1, 6, 2, 12, 8)
  xaxis_breaks <- c(3e+07, 8e+07, 1.3e+08)
  plot_width <- 3.8
  plot_height <- 5.5
}

### LOAD DATA ###
mets_cloneInfo <- fread("_assets/cloneInfo_df.csv")
setnames(mets_cloneInfo, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
mets_cloneInfo[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]

combined_df <- fread(file.path(output_directory, "cohort_outputs/combined.csv"))
alpaca <- combined_df[tumour_id == tum_id]

tum_out_dir <- file.path(save_directory, tum_id)
tum_seg_gene_path <- file.path(output_directory, "patient_outputs", tum_id, "tum_seg_gene_matched.csv")
tree_paths <- read_json(path = file.path(tum_out_dir, "tree_paths.json"))
tum_mets_cloneInfo <- mets_cloneInfo[tumour_id == tum_id]
colour_vec <- colour_seeding_tree(tree_paths, tum_mets_cloneInfo)

print(tum_out_dir)

### PLOT SAMPLES ###
plot_samples_alleles_gene(tum_id,
                          cohort,
                          run_name,
                          tum_out_dir,
                          tum_seg_gene_path,
                          gene_list,
                          chromosome,
                          min_startpos,
                          max_endpos,
                          asas_path,
                          sample_levels,
                          xaxis_breaks,
                          plot_width,
                          plot_height)


### PLOT CLONES ###
plot_clones_alleles_gene(tum_id,
                          cohort,
                          run_name,
                          alpaca,
                          tum_out_dir,
                          tum_seg_gene_path,
                          gene_list,
                          chromosome,
                          min_startpos,
                          max_endpos,
                          clone_levels,
                          colour_vec,
                          xaxis_breaks,
                          plot_width,
                          plot_height
)

###########
### END ###
###########
