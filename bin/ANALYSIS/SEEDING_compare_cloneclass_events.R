#!/usr/bin/env Rscript
# To run this script run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/SEEDING_compare_cloneclass_events.R \
--edge_cn_change_seg_chrarm_path output/mets/mets_default/cohort_outputs/edge_cn_change_seg_chrarm.csv \
--n_tums_event_thresh 5 \
--bin_directory bin \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/seeding_vs_nonseeding
"
cat("\nCopy number events in seeding vs non-seeding trajectories\n")

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(jsonlite))

### Command line options
option_list <- list(
  make_option(c("--edge_cn_change_seg_chrarm_path"), type="character", action="store", help="File path to combined cohort edge copy number change file with chrarm added info", metavar="character"),
  make_option(c("--n_tums_event_thresh"), type="character", action="store", metavar="character"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--output_directory"), type="character", help="output_directory"),
  make_option(c("--save_directory"), type="character", help="output_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
edge_cn_change_seg_chrarm_path <- opt$edge_cn_change_seg_chrarm_path
n_tums_event_thresh <- as.numeric(opt$n_tums_event_thresh)
bin_dir <- opt$bin_directory
output_directory <- opt$output_directory
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

assets_dir <- gsub('bin', '_assets', bin_dir)
patient_dir <- file.path(output_directory, "patient_outputs")


### Functions
source(file.path(bin_dir, "ANALYSIS/analysis_functions.R"))

### Load and process data
# Clinical data
clin_data_path <- file.path(assets_dir, "clinical_data_tumour.tsv")
clin_small <- load_reduced_clinicaldata(clin_data_path)

# Gene info
gene_loc_path <- file.path(assets_dir, "hg19refgene_tx421driveranno.csv")
gene_anno <- load_gene_anno(gene_loc_path)
print("Gene annotation file:")
print(head(gene_anno))

# ALPACA edge CN change data
edge_change_arm_df <- fread(edge_cn_change_seg_chrarm_path)
edge_change_arm_df[, clone := as.character(clone)]
edge_change_arm_df <- merge.data.table(edge_change_arm_df, clin_small, by = "tumour_id", all.x = TRUE, all.y = FALSE)

# Mets clone class groups:
cloneInfo_df <- fread(file.path(assets_dir, "cloneInfo_df.csv"))
setnames(cloneInfo_df, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
cloneInfo_df[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]
cloneInfo_df[, clone := as.character(clone)]

# Colour palette
tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))

#######################################################################
### CALL CN EVENTS IN MRCA -> SEEDING VS MRCA -> NON-SEEDING CLONES ###
#######################################################################
cat('\nCALL CN EVENTS IN MRCA -> SEEDING VS MRCA -> NON-SEEDING CLONES\n')

# Merge ALPACA with clone info:
edge_df <- merge.data.table(edge_change_arm_df, cloneInfo_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)

# Add clone class assignment
cat('\nAdd clone class assignment:\n')
edge_df <- assign_clone_classes(edge_clone_info_df = edge_df)

# Get all primary-met tumours where the mrca is not seeding (i.e. a subclone is seeding)
scseed_tumours <- sort(unique(edge_df[mrca_seeding_tum == FALSE, tumour_id]))

### Get gene-level CN events
cat("\n Getting gene level CN events in seeding and non-seeding primary clones\n")
gene_cn_events_df <- rbindlist(lapply(scseed_tumours, function(tum_id){
  print(tum_id)
  # get tumour specific data
  metstum_edge_change_arm_df <- edge_df[tumour_id == tum_id]
  # Compute MRCA -> Seeding and MRCA -> Non-seeding events
  seedingevents <- compute_mrca_seeding_nonseeding_events(metstum_edge_change_arm_df, assets_dir)
  # overlap with genes for plotting
  tum_seg_gene <- fread(file.path(patient_dir, tum_id, "tum_seg_gene_matched.csv"))
  alpaca_tum_gene <- merge.data.table(seedingevents, tum_seg_gene, by = c("tumour_id", "segment"), all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE)
  # filter for genes contained within segments:
  alpaca_tum_gene <- alpaca_tum_gene[overlap_type == "query_contained_in_subject"]
  alpaca_tum_gene[, chr := as.numeric(chr)]
  alpaca_tum_gene[, start := as.numeric(start)]
  alpaca_tum_gene[, end := as.numeric(end)]
  setorder(alpaca_tum_gene, chr, start, event_class)
  # subset cols
  gene_class_df <- unique(alpaca_tum_gene[, .(tumour_id, gene_name, chr, query_start, query_end, event_class, 
    is_amp, is_loh, is_loh_01, is_gain, is_loss, is_gain_1, is_loss_1, is_gain_2, is_loss_2,
    oncogene, tumour_suppressor, is_pancan_mut_driver, is_lung_mut_driver, Bailey_TSG_Oncogene, is_cn_driver, 
    cn_driver_class, mut_driver_class, is_pancan_driver
  )])
  setnames(gene_class_df, c("query_start", "query_end"), c("start", "end"))
  return(gene_class_df)
}))#3347243 rows

write.csv(gene_cn_events_df, file.path(save_directory, "gene_cn_events_df.csv"), row.names = FALSE, quote = FALSE)
# gene_cn_events_df <- fread("gene_cn_events_df.csv")

########################################
### COMPARE SEEDING VS NEVER SEEDING ###
########################################
cat("\n Summarising no. CN events of each type per gene, in seeding and non-seeding primary clones \n")

### Comparing seeding vs non-seeding
classes_to_compare <- c("MRCA --> Non-seeding", "MRCA --> Seeding")
events_to_test <- c("is_amp", "is_loh", "is_loh_01", "is_gain", "is_loss", "is_gain_1", "is_loss_1", "is_gain_2", "is_loss_2")

### Restrict to only patients who have both clone classes (tumours where the trunk is not seeding):
gene_cn_events_df[, pat_both_classes_pres := all(classes_to_compare %in% event_class), by = tumour_id]
gene_class_compare <- gene_cn_events_df[pat_both_classes_pres == TRUE & event_class %in% classes_to_compare] # 3234242 rows 

# summarise number of tumours with an event affecting each seeding class
gene_class_cn_events_summary <- summarise_event_counts(cn_events_df = gene_class_compare, events_to_test = events_to_test)

### Plot number of tumours which have seeding and non-seeding subclones
n_tums_total <- length(unique(edge_df[, tumour_id]))
n_tums_mrcanonseeding <- length(unique(gene_cn_events_df[, tumour_id]))
n_tums_seeding_nonseeding <- length(unique(gene_class_compare[, tumour_id]))

df <- data.table(group = c('Subclones only seeding,\nnon seeding paths present', 'Subclones only seeding,\nno non seeding paths present', 'MRCA seeding'),
                n = c(n_tums_seeding_nonseeding, n_tums_mrcanonseeding - n_tums_seeding_nonseeding, n_tums_total - n_tums_mrcanonseeding)
)

pie <- ggplot(df, aes(x="", y=n, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  theme_void()  +
  theme(legend.text = element_text(size = 16),
        legend.title = element_text(size = 18)) +
  scale_fill_manual(values = brewer.pal(3, "Blues"), name = "Tumour group") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 8)  
pdf(file.path(save_directory, "ntums_mrcanonseeding_pie.pdf"), width = 7, height = 7)
print(pie)
dev.off()

# Additionally save in figures/ directory
pdf("figures/Suppfig6e_ntums_mrcanonseeding_pie.pdf", width = 7, height = 7)
print(pie)
dev.off()

#################################
### STATISTICAL TEST PER GENE ###
#################################
cat("\n Statistical test per gene \n")

# add driver info
oncogenes <- gene_anno[(is_lung_mut_driver == TRUE & mut_driver_OG_TSG == "OG") | cn_driver_OG_TSG == "OG", unique(gene_name)]
tsgs <- gene_anno[(is_lung_mut_driver == TRUE & mut_driver_OG_TSG == "TSG") | cn_driver_OG_TSG == "TSG", unique(gene_name)]
print(oncogenes)
print(tsgs)

gene_class_cn_events_summary <- merge.data.table(gene_class_cn_events_summary, gene_anno, by = c("gene_name", "chr", "start", "end"), all.x = TRUE, all.y = FALSE)
gene_class_cn_events_summary[, gene_class := ifelse(gene_name %in% oncogenes, "oncogene", ifelse(gene_name %in% tsgs, "TSG", "passenger"))]
gene_class_cn_events_summary[, gene_class := factor(gene_class, levels = c("oncogene", "TSG", "passenger"))]
print(head(gene_class_cn_events_summary))

### AMPLIFICATIONS IN ONCOGENES
cat("\nTesting amps in oncogenes")
amps_oncogenes <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary,
                                                            event_col = "no_tums_amp",
                                                            driv_gene_list = oncogenes,
                                                            n_tums_event_thresh = n_tums_event_thresh)
print(amps_oncogenes)
# print(amps_oncogenes[fisher_pval < .05])

### GAINS IN ONCOGENES (> 1 copy)
cat("\nTesting gains > 1 copy in oncogenes")
gains_oncogenes <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_gain_1", driv_gene_list = oncogenes, n_tums_event_thresh = n_tums_event_thresh)
print(gains_oncogenes[fisher_pval < .05])

### GAINS IN ONCOGENES (> 2 copies)
cat("\nTesting gains > 2 copies in oncogenes")
gains2_oncogenes <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_gain_2", driv_gene_list = oncogenes, n_tums_event_thresh = n_tums_event_thresh)
print(gains2_oncogenes[fisher_pval < .05])

### LOSSES IN TSGs
cat("\nTesting losses > 1 copy in tsgs")
loss_tsgs <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_loss_1", driv_gene_list = tsgs, n_tums_event_thresh = n_tums_event_thresh)
print(loss_tsgs[fisher_pval < .2])


### LOH (state 0|1) IN TSGs
cat("\nTesting loh (0|1 state) in tsgs")
loh01_tsgs <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_loh_01", driv_gene_list = tsgs, n_tums_event_thresh = n_tums_event_thresh)
print(loh01_tsgs[fisher_pval < .05])

### LOH IN TSGs
cat("\nTesting loh (0|any state) in tsgs")
loh_tsgs <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_loh", driv_gene_list = tsgs, n_tums_event_thresh = n_tums_event_thresh)
print(loh_tsgs[fisher_pval < .05])


########################################################
### CROSS-GENOME PLOTS + SIGNIFICANT GENE ANNOTATION ###
########################################################
cat("\n Plotting cross-genome plots \n")

amp_sig_genes <- c(amps_oncogenes[fisher_pval < .05, unique(gene_name)])
gain_sig_genes <- c(gains_oncogenes[fisher_pval < .05, unique(gene_name)])
gain2_sig_genes <- c(gains2_oncogenes[fisher_pval < .05, unique(gene_name)])
loh_sig_genes <- c(loh_tsgs[fisher_pval < .05, unique(gene_name)])
loh01_sig_genes <- c(loh01_tsgs[fisher_pval < .05, unique(gene_name)])


cn_events_summary <- copy(gene_class_cn_events_summary)
cn_events_summary[, event_class := factor(event_class, levels = classes_to_compare)]


### Plot
# Plot fraction of tumours with event in each bin between MRCA and seeding:
# set colours
mets_colour_vec <- unlist(c(tx_palette$clone_class['Primary'], tx_palette$clone_class['Seeding']))
names(mets_colour_vec) <- classes_to_compare


# amplification vs loh
gain_1_loh_seeding_non <- plot_cross_genome_gain_loss_classes(
  assets_dir = assets_dir,
  crossgenome_df = cn_events_summary,
  max_cn = 40,
  gain_col = "frac_tums_gain_1",
  loss_col = "frac_tums_loh",
  gain_colour_vec = mets_colour_vec,
  loss_colour_vec = mets_colour_vec,
  plot_title = "% tumours with event",
  gain_anno = gain_sig_genes,
  loss_anno = loh_sig_genes,
  gain_title = "Gain > 1 copy",
  loss_title = "LOH"
)
ggsave(
  plot = gain_1_loh_seeding_non,
  file = file.path(save_directory, "mrca_seeding_nonseeding_cn_events_crossgenome_allgenes_gain_1_loh_genes.pdf"),
  width = 15, 
  height = 4
)
# Additionally save in figures/ directory
ggsave(
  plot = gain_1_loh_seeding_non,
  file = "figures/Suppfig6g_mrca_seeding_nonseeding_cn_events_crossgenome_allgenes_gain_1_loh_genes.pdf",
  width = 15, 
  height = 4
)

# gain vs loss
gain2_loh01_seeding_non <- plot_cross_genome_gain_loss_classes(
  assets_dir = assets_dir,
  crossgenome_df = cn_events_summary,
  max_cn = 30,
  gain_col = "frac_tums_gain_2",
  loss_col = "frac_tums_loh_01",
  gain_colour_vec = mets_colour_vec,
  loss_colour_vec = mets_colour_vec,
  plot_title = "% tumours with event",
  gain_anno = gain2_sig_genes,
  loss_anno = loh01_sig_genes,
  gain_title = "Gain > 2 copies",
  loss_title = "LOH (0|1 state)"
)
ggsave(
  plot = gain2_loh01_seeding_non,
  file = file.path(save_directory, "mrca_seeding_nonseeding_cn_events_crossgenome_allgenes_gain_2_loh01_genes.pdf"),
  width = 15, 
  height = 4
)

########################################
### PLOT ONLY SPECIFIC GENE EXAMPLES ###
########################################
cat("\n Plotting individual chromosome plots \n")

# CHR 11 CCND1 gain > 1 copy
chr11_gain_plot <- plot_chromosome_gene_anno(
  gene_list = c("CCND1"),
  chromosome = 11,
  crossgenome_df = copy(cn_events_summary),
  max_cn = 60,
  event_col = "frac_tums_gain_1",
  colour_vec = mets_colour_vec,
  assets_dir,
  y_title = "Gain > 1 copy",
  plot_title = "%tumours with event"
)
ggsave(
  plot = chr11_gain_plot,
  file = file.path(save_directory, "chr11_mrca_seeding_nonseeding_cn_events_gain1.pdf"),
  width = 3.2, height = 3.4
)
ggsave(
  plot = chr11_gain_plot,
  file = "figures/Fig4g_chr11_mrca_seeding_nonseeding_cn_events_gain1.pdf",
  width = 3.2, height = 3.4
)

# CHR 11 CCND1 gain > 2 copies
chr11_gain2_plot <- plot_chromosome_gene_anno(
  gene_list = c("CCND1"),
  chromosome = 11,
  crossgenome_df = copy(cn_events_summary),
  max_cn = 40,
  event_col = "frac_tums_gain_2",
  colour_vec = mets_colour_vec,
  assets_dir,
  y_title = "Gain > 2 copies",
  plot_title = "%tumours with event"
)
ggsave(
  plot = chr11_gain2_plot,
  file = file.path(save_directory, "chr11_mrca_seeding_nonseeding_cn_events_gain2.pdf"),
  width = 5, height = 2.5
)


# CHR 19
chr19_loh01_plot <- plot_chromosome_gene_anno(
  gene_list = loh01_sig_genes,
  chromosome = 19,
  crossgenome_df = copy(cn_events_summary),
  max_cn = 40,
  event_col = "frac_tums_loh_01",
  colour_vec = mets_colour_vec,
  assets_dir,
  y_title = "LOH (0|1 state)",
  plot_title = "%tumours with event"
)
ggsave(
  plot = chr19_loh01_plot,
  file = file.path(save_directory, "chr19_mrca_seeding_nonseeding_cn_events_loh01.pdf"),
  width = 5, height = 2.5
)

# CHR 18
chr18_loh01_plot <- plot_chromosome_gene_anno(
  gene_list = c("SMAD4"),
  chromosome = 18,
  crossgenome_df = copy(cn_events_summary),
  max_cn = 40,
  event_col = "frac_tums_loh_01",
  colour_vec = mets_colour_vec,
  assets_dir,
  y_title = "LOH (0|1 state)",
  plot_title = "%tumours with event"
)
ggsave(
  plot = chr18_loh01_plot,
  file = file.path(save_directory, "chr18_mrca_seeding_nonseeding_cn_events_loh.pdf"),
  width = 5, height = 2.5
)

##################################
### PLOT DRIVER GENE BAR PLOTS ###
##################################
cat("\n Plotting driver gene barplots \n")

### Run restricted significance testing:
### GAINS IN ONCOGENES
{
  # Restrict to genes with sufficient number of events
  # n_tums_event_thresh <- 10
  onco_test <- gene_class_cn_events_summary[gene_name %in% oncogenes]
  onco_test[, gene_to_test := sum(no_tums_gain_1) >= n_tums_event_thresh, by = gene_name]
  oncogenes_to_test <- onco_test[gene_to_test == TRUE, unique(gene_name)]
  # Restrict to genes with a proportion > bg proportion:
  bg_prop_seeding <- gains_oncogenes[gene_name == "background" & event_class == "MRCA --> Seeding", total_no_tums_gain_1]/gains_oncogenes[gene_name == "background", sum(total_no_tums_gain_1)]
  gain_1_bigprop <- gains_oncogenes[event_class == "MRCA --> Seeding" & prop > bg_prop_seeding, unique(gene_name)]
  oncogenes_to_test <- oncogenes_to_test[oncogenes_to_test %in% gain_1_bigprop]

  gains_oncogenes_final <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_gain_1", driv_gene_list = oncogenes_to_test, n_tums_event_thresh = n_tums_event_thresh)
  gains_oncogenes_final[gene_name != "background", p_adjust := p.adjust(fisher_pval, method = "fdr")]
  setorder(gains_oncogenes_final, fisher_pval, gene_name)
  print(gains_oncogenes_final)
}

### GAINS IN ONCOGENES > 2 copies 
{
  # Restrict to genes with sufficient number of events
  # n_tums_event_thresh <- 5
  onco_test <- gene_class_cn_events_summary[gene_name %in% oncogenes]
  onco_test[, gene_to_test := sum(no_tums_gain_2) >= n_tums_event_thresh, by = gene_name]
  oncogenes_to_test <- onco_test[gene_to_test == TRUE, unique(gene_name)]
  # Restrict to genes with a proportion > bg proportion:
  bg_prop_seeding <- gains2_oncogenes[gene_name == "background" & event_class == "MRCA --> Seeding", total_no_tums_gain_2]/gains2_oncogenes[gene_name == "background", sum(total_no_tums_gain_2)]
  gain_2_bigprop <- gains2_oncogenes[event_class == "MRCA --> Seeding" & prop > bg_prop_seeding, unique(gene_name)]
  oncogenes_to_test <- oncogenes_to_test[oncogenes_to_test %in% gain_2_bigprop]

  gains2_oncogenes_final <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_gain_2", driv_gene_list = oncogenes_to_test, n_tums_event_thresh = n_tums_event_thresh)
  gains2_oncogenes_final[gene_name != "background", p_adjust := p.adjust(fisher_pval, method = "fdr")]
  setorder(gains2_oncogenes_final, fisher_pval, gene_name)
  print(gains2_oncogenes_final)
}

### LOH IN TSGS
{
  # Restrict to genes with sufficient number of events
  # n_tums_event_thresh <- 10
  tsg_test <- gene_class_cn_events_summary[gene_name %in% tsgs]
  tsg_test[, gene_to_test := sum(no_tums_loh) >= n_tums_event_thresh, by = gene_name]
  tsgs_to_test <- tsg_test[gene_to_test == TRUE, unique(gene_name)]
  # Restrict to genes with a proportion > bg proportion:
  bg_prop_seeding <- loh_tsgs[gene_name == "background" & event_class == "MRCA --> Seeding", total_no_tums_loh]/loh_tsgs[gene_name == "background", sum(total_no_tums_loh)]
  loh_bigprop <- loh_tsgs[event_class == "MRCA --> Seeding" & prop > bg_prop_seeding, unique(gene_name)]
  tsgs_to_test <- tsgs_to_test[tsgs_to_test %in% loh_bigprop]

  loh_tsgs_final <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_loh", driv_gene_list = tsgs_to_test, n_tums_event_thresh = n_tums_event_thresh)
  loh_tsgs_final[gene_name != "background", p_adjust := p.adjust(fisher_pval, method = "fdr")]
  setorder(loh_tsgs_final, fisher_pval, gene_name)
  print(loh_tsgs_final)
}
### LOH (0|1 state) IN TSGS
{
  # Restrict to genes with sufficient number of events
  # n_tums_event_thresh <- 5
  tsg_test <- gene_class_cn_events_summary[gene_name %in% tsgs]
  tsg_test[, gene_to_test := sum(no_tums_loh_01) >= n_tums_event_thresh, by = gene_name]
  tsgs_to_test <- tsg_test[gene_to_test == TRUE, unique(gene_name)]
  # Restrict to genes with a proportion > bg proportion:
  bg_prop_seeding <- loh01_tsgs[gene_name == "background" & event_class == "MRCA --> Seeding", total_no_tums_loh_01]/loh01_tsgs[gene_name == "background", sum(total_no_tums_loh_01)]
  loh_bigprop <- loh01_tsgs[event_class == "MRCA --> Seeding" & prop > bg_prop_seeding, unique(gene_name)]
  tsgs_to_test <- tsgs_to_test[tsgs_to_test %in% loh_bigprop]

  loh01_tsgs_final <- run_paired_fishersexact_driv_vs_bg_event(gene_class_cn_events_summary, event_col = "no_tums_loh_01", driv_gene_list = tsgs_to_test, n_tums_event_thresh = n_tums_event_thresh)
  loh01_tsgs_final[gene_name != "background", p_adjust := p.adjust(fisher_pval, method = "fdr")]
  setorder(loh01_tsgs_final, fisher_pval, gene_name)
  print(loh01_tsgs_final)
}



# oncogene gains > 1 copy
plot_event_vs_bg_barplot(
  event_df = gains_oncogenes_final,
  plot_title = "Gains >1 in\noncogenes",
  out_filename = file.path(save_directory, "onco_gains1_gene_bar_proponly.seeding_vs_non.pdf"),
  colour_vec = mets_colour_vec,
  label_col = "no_tums_gain_1",
  plot_width = 3,
  plot_height = 0.3 * length(unique(gains_oncogenes_final$gene_name)) 
)
# oncogene gains > 2 copies
plot_event_vs_bg_barplot(
  event_df = gains2_oncogenes_final,
  plot_title = "Gains >2 in\noncogenes",
  out_filename = file.path(save_directory, "onco_gains2_gene_bar_proponly.seeding_vs_non.pdf"),
  colour_vec = mets_colour_vec,
  label_col = "no_tums_gain_2",
  plot_width = 3,
  plot_height = 0.3 * length(unique(gains2_oncogenes_final$gene_name)) 
)
# tsg lohs 
plot_event_vs_bg_barplot(
  event_df = loh_tsgs_final,
  plot_title = "LOH in\nTSGs",
  out_filename = file.path(save_directory, "tsg_loh_gene_bar_proponly.seeding_vs_non.pdf"),
  colour_vec = mets_colour_vec,
  label_col = "no_tums_loh",
  plot_width = 3,
  plot_height = 0.3 * length(unique(loh_tsgs_final$gene_name)) 
)
# tsg loh 0|1 state
plot_event_vs_bg_barplot(
  event_df = loh01_tsgs_final,
  plot_title = "LOH (0|1 state) in\nTSGs",
  out_filename = file.path(save_directory, "tsg_loh01_gene_bar_proponly.seeding_vs_non.pdf"),
  colour_vec = mets_colour_vec,
  label_col = "no_tums_loh_01",
  plot_width = 3,
  plot_height = 0.3 * length(unique(loh01_tsgs_final$gene_name)) 
)
# Plot legend:
events_legend <- ggplot(loss_tsgs,
  aes(gene_name, prop, fill = event_class)) +
    geom_bar(stat = "identity", position = position_fill()) +
    scale_fill_manual(values = mets_colour_vec, name = "Event class") +
    theme_cowplot() 
legend <- cowplot::get_legend(events_legend)
ggsave(filename = file.path(save_directory, "seeding_vs_non_legend.pdf"), plot = legend, width = 3, height = 1, limitsize = FALSE)

# Additionally save in figures directory:
plot_event_vs_bg_barplot(
  event_df = gains_oncogenes_final,
  plot_title = "Gains >1 in\noncogenes",
  out_filename = "figures/Fig4h_onco_gains1_gene_bar_proponly.seeding_vs_non.pdf",
  colour_vec = mets_colour_vec, 
  label_col = "no_tums_gain_1",
  plot_width = 3,
  plot_height = 0.3 * length(unique(gains_oncogenes_final$gene_name))
)
plot_event_vs_bg_barplot(
  event_df = gains2_oncogenes_final,
  plot_title = "Gains >2 in\noncogenes",
  out_filename = "figures/Suppfig6f_onco_gains2_gene_bar_proponly.seeding_vs_non.pdf",
  colour_vec = mets_colour_vec,
  label_col = "no_tums_gain_2",
  plot_width = 3,
  plot_height = 0.3 * length(unique(gains2_oncogenes_final$gene_name)) 
)
plot_event_vs_bg_barplot(
  event_df = loh_tsgs_final,
  plot_title = "LOH in\nTSGs",
  out_filename = "figures/Suppfig6f_tsg_loh_gene_bar_proponly.seeding_vs_non.pdf",
  colour_vec = mets_colour_vec,
  label_col = "no_tums_loh",
  plot_width = 3,
  plot_height = 0.3 * length(unique(loh_tsgs_final$gene_name)) 
)
plot_event_vs_bg_barplot(
  event_df = loh01_tsgs_final,
  plot_title = "LOH (0|1 state) in\nTSGs",
  out_filename = "figures/Suppfig6f_tsg_loh01_gene_bar_proponly.seeding_vs_non.pdf",
  colour_vec = mets_colour_vec,
  label_col = "no_tums_loh_01",
  plot_width = 3,
  plot_height = 0.3 * length(unique(loh01_tsgs_final$gene_name)) 
)

# extract exact p-values:
print(gains_oncogenes_final[fisher_pval < .05])
print(gains2_oncogenes_final[fisher_pval < .05])
print(loh_tsgs_final[fisher_pval < .05])
print(loh01_tsgs_final[fisher_pval < .05])

###########################
### SPLIT LUAD AND LUSC ###
###########################
cat("\n Testing splitting histologies \n")

histologies <- c("LUAD", "LUSC")
event_cols <- gsub("is", "no_tums", events_to_test)

histo_results <- lapply(histologies, function(h) {
  h_tumours <- clin_small[Histology == h, tumour_id]
  
  # Summarise event counts
  h_gene_class_cn_events_summary <- summarise_event_counts(gene_class_compare[tumour_id %in% h_tumours], events_to_test = events_to_test)
  h_gene_class_cn_events_summary <- merge.data.table(h_gene_class_cn_events_summary, gene_anno, by = c("gene_name", "chr", "start", "end"), all.x = TRUE, all.y = FALSE)
  h_gene_class_cn_events_summary[, gene_class := ifelse(gene_name %in% oncogenes, "oncogene", ifelse(gene_name %in% tsgs, "TSG", "passenger"))]
  h_gene_class_cn_events_summary[, gene_class := factor(gene_class, levels = c("oncogene", "TSG", "passenger"))]

  # Statistical test per gene
  events <- lapply(event_cols, function(event_col) {
    if(grepl("gain|amp", event_col)) {
      driv_gene_list <- oncogenes
    } else driv_gene_list <- tsgs  
    fish_results <- run_paired_fishersexact_driv_vs_bg_event(h_gene_class_cn_events_summary, event_col = event_col, driv_gene_list = driv_gene_list)
    return(fish_results[fisher_pval < .05])
  })
  names(events) <- event_cols
  return(list(h_gene_class_cn_events_summary = h_gene_class_cn_events_summary, sig_events = events))
})
names(histo_results) <- histologies


############################################
### Compare total # events in each group ###
############################################
cat("\n Plotting total #events across all genes for each event type in seeding vs non-seeding paths \n")
event_cols <- gsub("is", "no_tums", events_to_test)
genemelt <- melt.data.table(gene_class_cn_events_summary, measure.vars = event_cols, variable.name = 'event_type_var', value.name = 'no_tums_gene_event_inclass')
genemelt[, total_no_tums_gene_event := sum(no_tums_gene_event_inclass), by = .(gene_name, event_type_var)]
genemelt[, gene_event_proportion := no_tums_gene_event_inclass / total_no_tums_gene_event]
# head(genemelt)

### Filter for genes with >= 5 total events
genemelt_filt <- genemelt[total_no_tums_gene_event >= n_tums_event_thresh]
genemelt_filt[, event_type := gsub("no_tums_", "", event_type_var)]
genemelt_filt[, event_type := gsub("_1", " > 1 copy", event_type)]
genemelt_filt[, event_type := gsub("_2", " > 2 copies", event_type)]
genemelt_filt[event_type == "loh", event_type := "LOH (0|any state)"]
genemelt_filt[event_type == "loh_01", event_type := "LOH (0|1 state)"]
genemelt_filt[, event_type := factor(event_type, levels = c("amp", "gain > 2 copies", "gain > 1 copy", "gain", "LOH (0|any state)", "LOH (0|1 state)", "loss", "loss > 1 copy", "loss > 2 copies"))]
# head(genemelt_filt)

# event types to exclude:
gene_seeding_odds <- genemelt_filt[event_class == "MRCA --> Seeding"]
event_type_exclude <- c('loss > 1 copy')
gene_seeding_odds <- gene_seeding_odds[!event_type %in% event_type_exclude]

# Wilcoxon unpaired unparametric comparison between means of odds for events seeding branch in driver genes vs passengers
comps <- as.list(as.data.frame(combn(as.character(unique(gene_seeding_odds$gene_class)), 2)))
g <- ggplot(gene_seeding_odds, aes(gene_class, gene_event_proportion, fill = gene_class)) +
  geom_boxplot(outlier.size = .5, size = .3) +
  facet_wrap(~event_type, ncol = 1, strip.position = "left") +
  theme_cowplot() +
  labs(x = "SCNA event", y = "Odds of event occurring\non seeding vs non-seeding branch") +
  coord_flip() +
  geom_hline(yintercept = 0) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        panel.spacing = unit(.5, "lines"),
        # axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_manual(values = brewer.pal(3, "Paired"), name = 'Gene class') +
  stat_compare_means(comparisons = comps, label = "p.signif", hide.ns = TRUE, colour = "red", size = 3, tip.length = .01)
ggsave(
  filename = file.path(save_directory, "seedingBranchOdds_geneclassFilteredStats.pdf"),
  plot = g,
  width = 5.5,
  height = 4.5
)

#########################################################
### Line above plot up with
### Odds ratio SUMMING all genes for each gene class: ###
#########################################################
genemelt_filt[, n_eventtype_cloneclass_geneclass := sum(no_tums_gene_event_inclass), by = .(event_class, gene_class, event_type)]
genemelt_filt[, n_eventtype_cloneclass_all := sum(no_tums_gene_event_inclass), by = .(event_class, event_type)]

df <- unique(genemelt_filt[, .(event_class, gene_class, event_type, n_eventtype_cloneclass_geneclass, n_eventtype_cloneclass_all)])
df[, gene_class_to_compare := "TSG"]
df[grepl("amp|gain", event_type), gene_class_to_compare := "oncogene"]
df <- df[gene_class == gene_class_to_compare]
event_type_odds <- levels(df$event_type)[!levels(df$event_type) %in% c("loss > 1 copy", "loss > 2 copies")]

event_fish_results <- rbindlist(lapply(event_type_odds, function(e) {
  print(e)
  cont <- df[event_type == e, .(event_class, n_eventtype_cloneclass_all, n_eventtype_cloneclass_geneclass)]
  geneclasscompare <- df[event_type == e, unique(gene_class_to_compare)]
  cont_table <- as.matrix(cont)
  rownames(cont_table) <- cont_table[, 1]
  cont_table <- cont_table[, 2:ncol(cont_table)]
  cont_table <- apply(cont_table, c(1, 2), as.numeric)
  if (nrow(cont_table) == 2) {
    print(cont_table)
    ft <- fisher.test(cont_table)
    out_df <- data.table(event_type = e,
                          geneclasscompare = geneclasscompare,
                          odds_ratio = ft$estimate,
                          ci_lower = ft$conf.int[1],
                          ci_upper = ft$conf.int[2],
                          conf_level = attributes(ft$conf.int)[['conf.level']],
                          pval = ft$p.value)
    return(out_df)
  }
}))
print(event_fish_results)

event_fish_results[, event_type := factor(event_type, levels = event_fish_results[, event_type])]
event_fish_results[, result_signif := (ci_lower > 1) | (ci_upper < 1)]

## arrange both gene-level and summed results next to each other
gene_seed_odds_g <- ggplot(gene_seeding_odds, aes(gene_class, gene_event_proportion * 100, fill = gene_class)) +
  geom_boxplot(outlier.size = .5, size = .3) +
  facet_wrap(~event_type, ncol = 1, strip.position = "left") +
  theme_cowplot() +
  labs(x = "SCNA event", y = "% events affecting\ngene occuring\non seeding path") +
  coord_flip() +
  geom_hline(yintercept = 0) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        # strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.text.y.left = element_blank(),
        panel.spacing = unit(.5, "lines"),
        # axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values = brewer.pal(3, "Paired"), name = 'Gene class') +
  stat_compare_means(comparisons = comps, label = "p.format", hide.ns = TRUE, label_colour = "red", size = 4, tip.length = .01)

g_total_odds <- ggplot(event_fish_results,
  aes(forcats::fct_rev(event_type), odds_ratio, colour = geneclasscompare)) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width=.2), size = .4)+
  theme_cowplot() +
  scale_colour_manual(values = brewer.pal(3, "Paired"), name = 'Gene class') +
  geom_hline(yintercept = 1) +
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        legend.position = "none") +
  labs(x = "SCNA event", y = "Odds of event occurring\non seeding path in gene\nclass vs background") +
  coord_flip() +
  geom_text(aes(forcats::fct_rev(event_type), ci_upper + .2, label = "*"),  colour = "red", size = 7, data = event_fish_results[result_signif == TRUE])

g <- ggarrange(
  g_total_odds,
  gene_seed_odds_g,
  nrow = 1,
  widths = c(1.6, 2)
)
# Combine plots and save
ggsave(filename = file.path(save_directory, "seedingBranchOdds_total_genelevelGeneClass.pdf"), plot = g, width = 8, height = 5, limitsize = FALSE)
ggsave(filename = "figures/Fig4i_seedingBranchOdds_total_genelevelGeneClass.pdf", plot = g, width = 8, height = 5, limitsize = FALSE)

# Plot legend:
odds_legend <- ggplot(event_fish_results,
  aes(forcats::fct_rev(event_type), odds_ratio, colour = geneclasscompare)) +
  geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width=.2), size = .4)+
  theme_cowplot() +
  scale_colour_manual(values = brewer.pal(3, "Paired"), name = 'Gene class') +
  theme(legend.margin=margin(c(0,0,0,0)))
legend <- cowplot::get_legend(odds_legend)
ggsave(filename = file.path(save_directory, 'odds_legend.pdf'), plot = legend, width = 3, height = 3, limitsize = FALSE, bg = 'transparent')



# plot barplot summarising the background of all events:
df[, bg_sum := sum(n_eventtype_cloneclass_all), by = event_type]
df[, bg_frac := n_eventtype_cloneclass_all / bg_sum]
g <- ggplot(df,
  aes(event_type, bg_frac, fill = event_class)) +
  geom_col() +
  theme_cowplot() +
  geom_text(aes(label = n_eventtype_cloneclass_all), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = mets_colour_vec, name = "Tree path") +
  labs(x = "SCNA event", y = "Subclonal SCNA events affecting genes in each path") +
  coord_flip()
ggsave(
  filename = file.path(save_directory, "seedingBranchOdds_background.pdf"),
  plot = g,
  width = 8,
  height = 4.5
)

###########
### END ###
###########
