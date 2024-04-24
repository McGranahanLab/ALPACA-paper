#!/usr/bin/env Rscript
# To run this script run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/SEEDING_vs_nonseeding_plot_cross_genome_changes.R \
--edge_cn_change_seg_chrarm_path output/mets/mets_default/cohort_outputs/edge_cn_change_seg_chrarm.csv \
--bin_directory bin \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/seeding_vs_nonseeding
"

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(jsonlite))

### Command line options
option_list <- list(
  make_option(c("--edge_cn_change_seg_chrarm_path"), type="character", action="store", metavar="character"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--output_directory"), type="character", help="alpaca output dir"),
  make_option(c("--save_directory"), type="character", help="save_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
edge_cn_change_seg_chrarm_path <- opt$edge_cn_change_seg_chrarm_path
bin_dir <- opt$bin_directory
output_directory <- opt$output_directory
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

assets_dir <- gsub('bin', '_assets', bin_dir)
patient_dir <- file.path(output_directory, 'patient_outputs')

### Functions
source(file.path(bin_dir, "/ANALYSIS/analysis_functions.R"))

# Colour palette
tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))


### N gains / losses across the genome
edge_change_df <- fread(edge_cn_change_seg_chrarm_path)
seeding_vs_nonseeding_prim_df <- fread(file.path(assets_dir, 'seeding_vs_nonseeding_prim.csv'))

classes <- c("Non-seeding primary", "Seeding")
df <- merge.data.table(edge_change_df, seeding_vs_nonseeding_prim_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)
cn_events_summary <- df[!is.na(seeding_vs_nonseedingprim)]
cn_events_summary[, event_class := factor(seeding_vs_nonseedingprim, levels = classes)]


####################################################################
### Cross-genome plots of # CN CHANGES in seeding vs non-seeding ###
####################################################################
cat('\nCross-genome plots of # CN CHANGES in seeding vs non-seeding primary tumour clones\n')

### Get gene-level CN events
gene_cn_events_df <- rbindlist(lapply(unique(cn_events_summary$tumour_id), function(tum_id){
  print(tum_id)
  # get tumour specific data
  tum_edge_change_arm_df <- cn_events_summary[tumour_id == tum_id]
  
  # overlap with genes for plotting
  tum_seg_gene <- fread(file.path(patient_dir, tum_id, "tum_seg_gene_matched.csv"))
  alpaca_tum_gene <- merge.data.table(tum_edge_change_arm_df, tum_seg_gene, by = c("tumour_id", "segment"), all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE)
  
  # filter for genes contained within segments:
  alpaca_tum_gene <- alpaca_tum_gene[overlap_type == "query_contained_in_subject"]
  alpaca_tum_gene[, chr := as.numeric(chr)]
  alpaca_tum_gene[, start := as.numeric(start)]
  alpaca_tum_gene[, end := as.numeric(end)]
  setorder(alpaca_tum_gene, chr, start, event_class)

  # classify events
	alpaca_tum_gene[, gene_gain := any(any_gain), by = .(tumour_id, clone, gene_name, query_start, query_end)]
	alpaca_tum_gene[, gene_loss := any(any_loss), by = .(tumour_id, clone, gene_name, query_start, query_end)]

	# subset cols
  gene_class_df <- unique(alpaca_tum_gene[, .(tumour_id, 
                                                clone,
                                                gene_name,
                                                chr,
                                                start = query_start,
																								end = query_end,
                                                event_class, 
                                                gene_gain, 
                                                gene_loss, 
                                                oncogene, tumour_suppressor, is_pancan_mut_driver, is_lung_mut_driver, Bailey_TSG_Oncogene, is_cn_driver, 
                                                cn_driver_class, mut_driver_class, is_pancan_driver
  )])
  
  ### Summarise total number of clones in tumour with an event affecting each seeding class
  gene_class_df[, n_tumclones_gain := sum(gene_gain), by = .(event_class, gene_name)]
  gene_class_df[, n_tumclones_loss := sum(gene_loss), by = .(event_class, gene_name)]
  gene_class_df[, tumclone := paste(tumour_id, clone, sep = ":")]
  gene_class_df[, n_tumclones_group := length(unique(tumclone)), by = .(event_class)]
  
  out_df <- unique(gene_class_df[, .(tumour_id, gene_name, chr, start, end, event_class, n_tumclones_gain, n_tumclones_loss, n_tumclones_group)])
  return(out_df)
}))
# print(paste0("nrow(gene_cn_events_df) = ", nrow(gene_cn_events_df))) #nrow = 3347243

### Summarise total number of clones with an event affecting each seeding class
cat('\nSummarise total number of clones with an event affecting each seeding class\n')

gene_cn_events_df[, n_clones_group := sum(n_tumclones_group), by = .(event_class, gene_name)]
gene_cn_events_df[, n_clones_gain := sum(n_tumclones_gain), by = .(event_class, gene_name)]
gene_cn_events_df[, n_clones_loss := sum(n_tumclones_loss), by = .(event_class, gene_name)]
gene_cn_events_df[, frac_clones_gain := n_clones_gain / n_clones_group]
gene_cn_events_df[, frac_clones_loss := n_clones_loss / n_clones_group]

cn_summary <- unique(gene_cn_events_df[, .(gene_name, chr, start, end, event_class, n_clones_gain, n_clones_loss, frac_clones_gain, frac_clones_loss)])
setorder(cn_summary, chr, start, end)

# set colours
seeding_col_vec <- unlist(c(tx_palette$clone_class['Primary'], tx_palette$clone_class['Seeding']))
names(seeding_col_vec) <- classes

cat('\nPlot cross-genome plot\n')
g <- plot_cross_genome_gain_loss_classes(
  assets_dir = assets_dir,
  crossgenome_df = cn_summary,
  max_cn = 40,
  gain_col = "frac_clones_gain",
  loss_col = "frac_clones_loss",
  gain_colour_vec = seeding_col_vec,
  loss_colour_vec = seeding_col_vec,
  plot_title = "% clones with event",
  gain_title = "% Clones with gain",
  loss_title = "% Clones with loss",
  plot_legend = FALSE
)
ggsave(
  plot = g,
  file = paste0(save_directory, "/frac_clones_gainloss_seeding_nonseedingprim.pdf"),
  width = 13, 
  height = 4.5
)
# Also save in figures directory
ggsave(plot = g,
  file = "figures/Fig4d_frac_clones_gainloss_seeding_nonseedingprim.pdf",
  width = 13, 
  height = 4.5
)


### What fraction of the genome is larger in seeding vs non-seeding

cat('\nWhat fraction of the genome is larger in seeding vs non-seeding?\n')
gain_df <- dcast(cn_summary, formula = gene_name + chr + start + end ~ event_class, value.var = 'frac_clones_gain')
loss_df <- dcast(cn_summary, formula = gene_name + chr + start + end ~ event_class, value.var = 'frac_clones_loss')

gain_df[, larger_in_seeding := Seeding > `Non-seeding primary`]
loss_df[, larger_in_seeding := Seeding > `Non-seeding primary`]


setorder(gain_df, chr, start)
setorder(loss_df, chr, start)

gain_df[, genome_size := sum(end - start)]
loss_df[, genome_size := sum(end - start)]

gain_df[, seg_size := end - start]
loss_df[, seg_size := end - start]

cat('\nFraction of genome larger in seeding clones for gains:\n')
unique(gain_df[, .(frac_genome_seeding_larger = sum(seg_size) / genome_size), by = larger_in_seeding])

cat('\nFraction of genome larger in seeding clones for losses:\n')
unique(loss_df[, .(frac_genome_seeding_larger = sum(seg_size) / genome_size), by = larger_in_seeding])

