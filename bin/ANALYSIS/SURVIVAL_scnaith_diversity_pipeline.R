#!/usr/bin/env Rscript
# To run this script run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/SURVIVAL_scnaith_diversity_pipeline.R \
--clone_distances_path output/primary/primary_default/cohort_outputs/edge_events.csv \
--pairwise_dist_path output/primary/primary_default/cohort_outputs/cn_dist_pairwise.csv \
--bin_directory bin \
--mets_output_directory output/mets/mets_default \
--save_directory output/primary/primary_default/cohort_outputs/scna_ith
"

#############################################################
### SCRIPT TO ANALYSE SCNA-ITH ASSOCIATIONS WITH SURVIVAL ###
#############################################################
cat('\nSCRIPT TO ANALYSE SCNA-ITH ASSOCIATIONS WITH SURVIVAL\n')

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(scales))

### Command line options
option_list <- list(
  make_option(c("--clone_distances_path"), type="character", action="store", help="File path to combined cohort clone genomic distances file", metavar="character"),
  make_option(c("--pairwise_dist_path"), type="character", action="store", help="File path to combined cohort clone genomic distances file", metavar="character"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--mets_output_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--save_directory"), type="character", help="output plots location")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)


### Parse args
clone_distances_path <- opt$clone_distances_path
pairwise_dist_path <- opt$pairwise_dist_path
bin_dir <- opt$bin_directory
mets_output_directory <- opt$mets_output_directory
save_directory <- opt$save_directory

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

assets_dir <- gsub('bin', '_assets', bin_dir)

### Functions
source(file.path(bin_dir, "/ANALYSIS/analysis_functions.R"))
source(file.path(bin_dir, "/ANALYSIS/SURVIVAL_functions.R"))

### Data
# Clinical data
all_patient_df <- fread(file.path(assets_dir, "20221109_TRACERx421_all_patient_df.tsv"))

# Clone genomic distances data
distances_df <- fread(clone_distances_path)
pairwise_dist <- fread(pairwise_dist_path)

# Non-metastatic patients
non_mets_pats <- fread(file.path(assets_dir, 'nonMetastaticPatients.txt'), header = FALSE)
non_mets_pats <- as.character(non_mets_pats$V1)

# ln clones
ln_prims <- fread(file.path(assets_dir, 'lymphnode_primaries.csv'))

# evo metrics
evo_df <- fread(file.path(assets_dir, "20221110_TRACERx421_evolutionary_metrics.tsv"))

tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))

boxplot_theme <- theme(strip.background = element_rect(fill = "transparent", colour = "black"),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))

# Seeding info
seeding_df <- fread(file.path(assets_dir, "seedingTableFull.csv"))

# Relapse category
patient_df <- fread(file.path(assets_dir, "clinical_data.tsv"))
relapse_df <- unique(patient_df[, .(patient_id = `Patient ID`, `Relapse category`)])

############
### MAIN ###
############
cat('\nProcess data for survival analysis:\n')

### Process edge event data for survival analysis:
distances_df <- process_edge_snv_scna_dist(distances_df)
dist_df <- process_edgeevents_for_survival_analysis(distances_df, evo_df)

### Process pairwise CN distance data for survival analysis:
pairwise_dist_sc <- process_pairwisecndist_for_survival_analysis(pairwise_dist, ln_prims)


### Combine all metrics and summarise on a tumour level
tum_df_allsizes <- merge.data.table(
  x = unique(dist_df[parent_clone != "diploid", .(tumour_id,
                                                  `SCNA-ITH (sample)` = frac_abberant_genom_subcl,
                                                  `Mean #SCNAs` = n_scnas_per_clone,
                                                  sc_GD,
                                                  genom_frac_event,
                                                  slope_norm_meantum,
                                                  n_amplitude_scnas_per_clone,
                                                  n_amplitude_gains_per_clone,
                                                  n_amplitude_losses_per_clone,
                                                  frac_subclonal_scnas,
                                                  frac_subclonal_gains,
                                                  frac_subclonal_losses,
                                                  purity_bootstrapped,
                                                  ploidy_bootstrapped,
                                                  n_clones,
                                                  n_segs,
                                                  slope_variance)]),
  y = unique(pairwise_dist_sc[, .(tumour_id,
                                  seg_size,
                                  maxtum_l1_dist_pairwise_allelespecific,
                                  maxtum_l1_w_dist_pairwise_allelespecific,
                                  maxtum_l2_dist_pairwise_allelespecific,
                                  maxtum_l2_w_dist_pairwise_allelespecific,
                                  maxtum_l2_dist_pairwise_sumallele,
                                  CCD = maxtum_l2_dist_pairwise_AB)]),
  by = "tumour_id", all.x = TRUE, all.y = FALSE
)
tum_df_raw <- copy(tum_df_allsizes)

### Summarise for all segments: 
tum_df <- tum_df_allsizes[seg_size == "all"]
# write.csv(tum_df, file.path(assets_dir, 'scna_div_metrics.csv'), row.names = FALSE, quote = FALSE)

###########################################
### Multi-Cox full model -> forest plot ### 
###########################################
cat('\nRun survival analysis for full primary cohort\n')

metrics_list <- c("SCNA-ITH (sample)", "CCD")
in_data_surv <- preprocess_survival_data(all_patient_df, tum_df, metrics_list)

# Plot K-M by median
out_pdf <- file.path(save_directory, "kmplot.pdf")
group_vec <- paste0(metrics_list, "_cat")
plot_km_bymedian(in_data_surv, group_vec, out_pdf)
plot_km_bymedian(in_data_surv, group_vec, out_pdf = "figures/Fig5e_Suppfig6a_kmplot.pdf")

# Plot K-M by tertiles
out_pdf <- file.path(save_directory, "kmplot.pdf")
group_vec <- paste0(metrics_list, "_tertiles")
plot_km_bytertiles(in_data_surv, group_vec, out_pdf)
plot_km_bytertiles(in_data_surv, group_vec, out_pdf = "figures/Fig5f_Suppfig6b_kmplot.pdf")

### Run MV analysis
metrics_list <- c("SCNA-ITH (sample)", "CCD")
in_data_surv <- preprocess_survival_data(all_patient_df, tum_df, metrics_list)
setnames(in_data_surv, 'CCD_by_sd', "Clone CN diversity (per s.d.)")
setnames(in_data_surv, 'pTNM_stage', 'pTNM stage')
setnames(in_data_surv, 'histology_LUAD', 'Histology')
setnames(in_data_surv, 'adjuvant_tx', 'Adjuvant treatment')
setnames(in_data_surv, 'SCNA-ITH (sample)_by_sd', 'SCNA-ITH (sample, per s.d.)')
variables_cox <- c("Age (+10 years)", "pTNM stage", "PackYears (+10)", "Histology", "sex", "Adjuvant treatment", "Clone CN diversity (per s.d.)", "SCNA-ITH (sample, per s.d.)")
out_pdf <- file.path(save_directory, "forestplot_CCD_by_sd_SCNAITHsample.noLNclones.TRACERxcovariates.pdf")
run_multivariate_analysis(in_data_surv, variables_cox, out_pdf)
run_multivariate_analysis(in_data_surv, variables_cox, out_pdf = "figures/Fig5g_forestplot_CCD_by_sd_SCNAITHsample.noLNclones.TRACERxcovariates.pdf")


#######################################################
### Run survival analysis with LUAD/LUSC separately ###
#######################################################
cat('\nRun survival analysis with LUAD/LUSC separately\n')

metrics_list <- c("CCD")

### LUAD
luad_in_data_surv <- preprocess_survival_data(all_patient_df[histology_lesion1 == "Invasive adenocarcinoma"], tum_df, metrics_list)
out_pdf <- "figures/Suppfig6c_kmplot_LUAD.pdf"
group_vec <- paste0(metrics_list, "_cat")
plot_km_bymedian(luad_in_data_surv, group_vec, out_pdf, title_start = "DFS - LUAD")
# run_multivariate_analysis(in_data_surv = luad_in_data_surv, variables_cox = c("Age (+10 years)", "pTNM_stage", "PackYears (+10)", "sex", "adjuvant_tx", "CCD"), out_pdf = file.path(save_directory, "kmplot_CCD.noLNclones_ptnm.LUAD.pdf"))


### LUSC
lusc_in_data_surv <- preprocess_survival_data(all_patient_df[histology_lesion1 == "Squamous cell carcinoma"], tum_df, metrics_list)
out_pdf <- "figures/Suppfig6d_kmplot_LUSC.pdf"
group_vec <- paste0(metrics_list, "_cat")
plot_km_bymedian(lusc_in_data_surv, group_vec, out_pdf, title_start = "DFS - LUSC")
# run_multivariate_analysis(in_data_surv = lusc_in_data_surv, variables_cox = c("Age (+10 years)", "pTNM_stage", "PackYears (+10)", "sex", "adjuvant_tx", "CCD"), out_pdf = file.path(save_directory, "kmplot_CCD.noLNclones_ptnm.LUSC.pdf"))


########################
### Further analyses ###
########################

### How are CCD and tumour relapse related? 
tum_df[, cruk_id := gsub("-Cluster.*", "", tumour_id)]
tum_df[, is_non_metastatic_patient := cruk_id %in% non_mets_pats]
tum_df[, is_metastatic := tumour_id %in% list.files(file.path(mets_output_directory, "patient_outputs"))]
tum_df[is_metastatic == TRUE, recurrence_status := "Metastatic"]
tum_df[is_non_metastatic_patient == TRUE, recurrence_status := "Non-metastatic"]
recur_colour_vec <- unlist(c(tx_palette$clone_class['Primary'], tx_palette$clone_class['Seeding']))
names(recur_colour_vec) <- c("Non-metastatic", "Metastatic")
# Figure 6c
g <- ggplot(tum_df[!is.na(recurrence_status)], aes(recurrence_status, CCD)) +
  geom_boxplot(aes(fill = recurrence_status)) +
  theme_cowplot() +
  scale_fill_manual(values = recur_colour_vec, name = "Patient group") +
  labs(x = "Recurrence status") +
  stat_compare_means(size = 6, label.y = 155) +
  theme(legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))
pdf(file.path(save_directory, 'recurrence_status_CCD.pdf'), width = 3.5, height = 4)
print(g)
dev.off()
pdf('figures/Fig5c_recurrence_status_CCD.pdf', width = 3.5, height = 4)
print(g)
dev.off()
# How many patients of each type?
table(tum_df[!is.na(recurrence_status), recurrence_status])


### Do truncal seeding tumours have higher CCD?
# Figure 6d
# Mets clone class groups:
mets_cloneInfo <- fread(file.path(assets_dir, "cloneInfo_df.csv"))
setnames(mets_cloneInfo, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
mets_cloneInfo[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]
mets_cloneInfo[, truncal_seeding := ifelse(any(clonalClones == TRUE & seedingClones == TRUE), "Truncal", "Subclonal"), by = tumour_id]
tum_df_seed <- merge.data.table(tum_df, unique(mets_cloneInfo[, .(tumour_id, truncal_seeding)]), by = "tumour_id", all = FALSE)
trunkseed_colour_vec <- unlist(c(tx_palette$clone_class['Seeding'], tx_palette$clone_class['Shared']))
names(trunkseed_colour_vec) <- c("Subclonal", "Truncal")

g <- ggplot(tum_df_seed, aes(truncal_seeding, CCD)) +
  geom_boxplot(aes(fill = truncal_seeding)) +
  theme_cowplot() +
  scale_fill_manual(values = trunkseed_colour_vec) +
  labs(x = "Tumour seeding status") +
  stat_compare_means(size = 6, label.y = 155) +
  theme(legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))
pdf(file.path(save_directory, 'truncal_seeding_CCD.pdf'), width = 3.5, height = 4)
print(g)
dev.off()
pdf('figures/Fig5d_truncal_seeding_CCD.pdf', width = 3.5, height = 4)
print(g)
dev.off()
# How many patients of each type?
table(tum_df_seed[!is.na(truncal_seeding), truncal_seeding])


### How do CCD and SCNA-ITH density across the cohort compare? 
cmap <- colorRampPalette(c("white", "lightblue", "blue2"))
g <- ggplot(in_data_surv, aes(`SCNA-ITH (sample)`, CCD)) +
  geom_density_2d_filled(contour_var = 'ndensity', alpha = 1)+ 
  geom_point(size = .8) +
  # facet_wrap(~TNM_stage, nrow = 1, scales = "free") +
  geom_smooth() +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "transparent", colour = "black")) +
  scale_fill_manual(values = cmap(13)) +
  stat_cor()
pdf(file.path(save_directory, "scnaith_vs_CCD.pdf"), width = 4, height = 3.8)
print(g)
dev.off()


# How do they correlate across stages?
# Supp figure 6e
cmap <- colorRampPalette(c("white", "lightblue", "blue2"))
g <- ggplot(in_data_surv, aes(`SCNA-ITH (sample)`, CCD)) +
  geom_density_2d_filled(contour_var = 'ndensity', alpha = 1)+ 
  geom_point(size = .8) +
  facet_wrap(~TNM_stage, nrow = 1, scales = "free") +
  geom_smooth(method = "lm") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "transparent", colour = "black")) +
  scale_fill_manual(values = cmap(13)) +
  stat_cor()
pdf(file.path(save_directory, "scnaith_vs_CCD_stage.pdf"), width = 8, height = 3.8)
print(g)
dev.off()
pdf("figures/Suppfig6e_scnaith_vs_CCD_stage.pdf", width = 8, height = 3.8)
print(g)
dev.off()



# Compare the density plots of the two metrics
# First normalise metrics by doing z-score
# Supp fig 6f
in_data_surv[, `SCNA-ITH (sample)_norm` := (`SCNA-ITH (sample)` - mean(`SCNA-ITH (sample)`)) / `SCNA-ITH (sample)_sd`]
in_data_surv[, CCD_norm := (CCD - mean(CCD)) / CCD_sd]
metrics_dens_compare <- melt.data.table(in_data_surv, measure.vars = c('SCNA-ITH (sample)_norm', 'CCD_norm'))
metrics_dens_compare[, metric := gsub("_norm", "", variable)]
metrics_dens_compare[metric == "SCNA-ITH (sample)", metric := "SCNA-ITH\n(sample)"]
metrics_dens_compare[, metric := factor(metric, levels = c("SCNA-ITH\n(sample)", "CCD"))]
g <- ggplot(metrics_dens_compare, aes(value, metric, fill = TNM_stage)) +
  geom_density_ridges(alpha = .7, linewidth=.6) +
  theme_cowplot() +
  theme(axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "bottom") +
  scale_fill_manual(values = unlist(tx_palette$stage_short), name = "NSCLC stage") +
  scale_y_discrete(expand = c(0, .1)) +
  labs(x = "Normalised z-score", y= 'Metric')
pdf(file.path(save_directory, "density_scnaith_CCD_stage_ridges.pdf"), width = 5, height = 5)
print(g)
dev.off()
pdf("figures/Suppfig6f_density_scnaith_CCD_stage_ridges.pdf", width = 5, height = 5)
print(g)
dev.off()


# Compare N scnas per clone in mono vs polyclonal:
# Figure 4j
# Merge with seeding clonality:
tum_df_seeding <- merge.data.table(tum_df, seeding_df, by = "tumour_id", all.x = TRUE, all.y = FALSE)
# Merge with relapse:
tum_df_seeding[, patient_id := gsub('-Cluster.*', '', tumour_id)]
tum_df_seeding <- merge.data.table(tum_df_seeding, relapse_df, by = "patient_id", all.x = TRUE, all.y = FALSE)

tum_df_seeding_melt1 <- melt.data.table(tum_df_seeding, measure.vars = c("SCNA-ITH (sample)", "Mean #SCNAs"))#, "CCD"))
tum_df_seeding_melt1[variable == "SCNA-ITH (sample)", variable := "SCNA-ITH\n(sample)"]
tum_df_seeding_melt1[, variable := factor(variable, levels = c("SCNA-ITH\n(sample)", "Mean #SCNAs"))]
g <- ggplot(tum_df_seeding_melt1[!is.na(seedingClonality)], 
  aes(seedingClonality, value, fill = seedingClonality)) +
  facet_wrap(~variable, nrow = 1, scales = "free_y") + 
  theme_cowplot() +
  scale_fill_manual(values = unlist(tx_palette$seedingClonality), name = "Tumour seeding clonality") +
  labs(x = "Tumour seeding clonality") +
  geom_jitter(size = 0.5, width = .25, height = 0, colour = "grey50")+
  geom_boxplot(alpha = .9, outlier.shape = NA) +
  boxplot_theme +
  guides(fill = guide_legend(nrow = 2)) +
  stat_compare_means(size = 4)
pdf("figures/Fig4j_seedingClonality_vs_MeanNoSCNAs.scnaITH.pdf", width = 5, height = 4)
print(g)
dev.off()

###########
### END ###
###########