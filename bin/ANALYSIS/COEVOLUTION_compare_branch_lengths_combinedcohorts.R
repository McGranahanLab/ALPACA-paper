#!/usr/bin/env Rscript

# To run this script run the following command from the command line (make sure you are in the root directory of this repo): 
"
Rscript bin/ANALYSIS/COEVOLUTION_compare_branch_lengths_combinedcohorts.R \
--mets_clone_distances_path output/mets/mets_default/cohort_outputs/edge_events.csv \
--prim_clone_distances_path output/primary/primary_default/cohort_outputs/edge_events.csv \
--bin_directory bin \
--save_directory output/mets/mets_default/cohort_outputs/n_events
"
cat("\nComparing SCNA and SNV branch lengths\n")

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(lmerTest))


### Command line options
option_list <- list(
  make_option(c("--mets_clone_distances_path"), type="character", action="store", help="File path to combined mets cohort clone genomic distances file", metavar="character"),
  make_option(c("--prim_clone_distances_path"), type="character", action="store", help="File path to combined primary cohort clone genomic distances file", metavar="character"),
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--save_directory"), type="character", help="Plot output directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
mets_clone_distances_path <- opt$mets_clone_distances_path
prim_clone_distances_path <- opt$prim_clone_distances_path
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
mets_edge_distances <- fread(mets_clone_distances_path)
prim_edge_distances <- fread(prim_clone_distances_path)

### Colour palette
tx_palette <- read_json(path = file.path(assets_dir, "publication_palette.json"))

############
### MAIN ###
############

prim_edge_counts_df <- process_edge_snv_scna_dist(edge_distance_df = prim_edge_distances)
mets_edge_counts_df <- process_edge_snv_scna_dist(edge_distance_df = mets_edge_distances)

### Merge mets and prim:
prim_edge_counts_df[, cohort := "Tx421 primary"]
mets_edge_counts_df[, cohort := "Tx421 primary-metastasis"]

prim_met_edge_df <- rbind(prim_edge_counts_df, mets_edge_counts_df, fill = TRUE)

### Merge with clinical info:
prim_met_edge_df[tumour_id == "CRUK0084-Cluster2", tumour_id := "CRUK0084"] 
prim_met_edge_df <- merge.data.table(prim_met_edge_df,
                                    clin_small,
                                    by = "tumour_id",
                                    all.x = TRUE,
                                    all.y = FALSE)

### Is mut branch length correlated with CN branch length?
g <- ggplot(prim_met_edge_df[is_mrca == 'subclone'],
            aes(mut_branchlength_norm_subclonal, total_cn_events_norm_subclonal,
                colour = Histology)) +
  geom_point(size = .8) +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_cowplot() +
  stat_cor(label = "p") +
  facet_wrap(~ cohort) +
  scale_colour_manual(values = tx_palette$histology_subtypes) +
  geom_abline(linetype = 'dashed', intercept = 0, slope = 1) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        strip.background = element_rect(fill = "white", colour = "black"))+
  labs(x = 'Subclonal #SNVs',
        y = 'Subclonal #SCNAs',
        title = NULL
    )
pdf(file.path(save_directory, "prim_mets_snv_scna_branch_length.pdf"), width = 7, height = 3.5)
print(g)
dev.off()

# Also save in figures directory
pdf("figures/Fig4c_prim_mets_snv_scna_branch_length.pdf", width = 7, height = 3.5)
print(g)
dev.off()

# calculate p-values for correlation between mut and cn branch lengths separately to avoid R p<2.2e−16 values
library(dplyr)
LUAD = prim_met_edge_df[prim_met_edge_df$Histology == "LUAD"]
LUSC = prim_met_edge_df[prim_met_edge_df$Histology == "LUSC"]
OTHER = prim_met_edge_df[prim_met_edge_df$Histology == "Other"]


pvals_cor_LUAD <- LUAD %>%
  filter(is_mrca == 'subclone') %>%
  group_by(cohort) %>%
  summarise(
    cor_test = list(cor.test(mut_branchlength_norm_subclonal, total_cn_events_norm_subclonal)),
    .groups = 'drop'
  ) %>%
  mutate(
    estimate = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) format(x$p.value, digits = 3, scientific = TRUE))
  ) %>%
  select(cohort, estimate, p_value)
print(pvals_cor_LUAD)

pvals_cor_LUSC <- LUSC %>%
  filter(is_mrca == 'subclone') %>%
  group_by(cohort) %>%
  summarise(
    cor_test = list(cor.test(mut_branchlength_norm_subclonal, total_cn_events_norm_subclonal)),
    .groups = 'drop'
  ) %>%
  mutate(
    estimate = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) format(x$p.value, digits = 3, scientific = TRUE))
  ) %>%
  select(cohort, estimate, p_value)
print(pvals_cor_LUSC)

pvals_cor_OTHER <- OTHER %>%
  filter(is_mrca == 'subclone') %>%
  group_by(cohort) %>%
  summarise(
    cor_test = list(cor.test(mut_branchlength_norm_subclonal, total_cn_events_norm_subclonal)),
    .groups = 'drop'
  ) %>%
  mutate(
    estimate = sapply(cor_test, function(x) x$estimate),
    p_value = sapply(cor_test, function(x) format(x$p.value, digits = 3, scientific = TRUE))
  ) %>%
  select(cohort, estimate, p_value)
print(pvals_cor_OTHER)

#############################
### Linear model analysis ###
#############################

### Predicting SCNA branch length per clone in\nTRACERx421 primary-metastasis cohort (Linear mixed-effects model)
# Reading in mets clone class groups:
cloneInfo_df <- fread(file.path(assets_dir, "cloneInfo_df.csv"))
setnames(cloneInfo_df, c("SampleID", "PyCloneCluster_SC"), c("tumour_id", "clone"))
cloneInfo_df[, tumour_id := gsub("_Cluster", "-Cluster", tumour_id)]
cloneInfo_df[, clone := as.character(clone)]
# Merge ALPACA with clone info:
mets_edge_df <- merge.data.table(mets_edge_counts_df, cloneInfo_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)

# Add clone class assignment
mets_edge_df <- assign_clone_classes(edge_clone_info_df = mets_edge_df)
mets_edge_df <- merge.data.table(mets_edge_df, clin_small, by = "tumour_id", all.x = TRUE, all.y = FALSE)

# Linear mixed effects model
mets_lm_full <- lmer(total_interval_events_bin ~ mut_branchlength + clone_class + Histology + (1|tumour_id), data = mets_edge_df[is_mrca == 'subclone'])
mets_lm_reduced <- lmer(total_interval_events_bin ~ Histology + clone_class + (1|tumour_id), data = mets_edge_df[is_mrca == 'subclone'])
anova(mets_lm_full, mets_lm_reduced)
summary(mets_lm_full)

# Specify lme model with normalised clone #events
mets_lme_simple <- lmerTest::lmer(total_cn_events_norm_subclonal ~ mut_branchlength_norm_subclonal
                                                        + clone_class
                                                        + Histology
                                                        + (1 | tumour_id)
                    , data = mets_edge_df[is_mrca == 'subclone'])
summary(mets_lme_simple)

model_df <- as.data.table(coef(summary(mets_lme_simple)), keep.rownames = 'Term')
model_df <- setNames(model_df, c("Term", "Estimate", "std_error", "df", "t_value", "p_value"))
model_df[p_value < 0.05, p_signif := "*"]
model_df[p_value < 0.01, p_signif := "**"]
model_df[p_value < 0.001, p_signif := "***"]

confint_df <- as.data.table(confint(mets_lme_simple), keep.rownames = 'Term')
model_df <- merge.data.table(model_df, confint_df, by = 'Term', all = FALSE)

model_df[Term == 'mut_branchlength_norm_subclonal', Term := '#SNVs']
model_df[grepl("clone_class", Term), Term := paste0(Term, " clone")]
model_df[grepl("clone_class", Term), Term := gsub("clone_class", "", Term)]

### Plot linear model results
colour_pal <- c('#SNVs' = brewer.pal(10, 'Paired')[1], 
                'HistologyLUSC' = as.character(unlist(tx_palette$histology)['LUSC']), 
                'HistologyOther' = as.character(unlist(tx_palette$histology)['Other']),
                'Seeding clone' = as.character(unlist(tx_palette$clone_class)['Seeding']),
                'Metastatic clone' = as.character(unlist(tx_palette$clone_class)['Metastatic']),
                'Primary clone' = as.character(unlist(tx_palette$clone_class)['Primary']),
                '(Intercept)' = "grey50")

g <- ggplot(model_df,
    aes(Term, Estimate, colour = Term)) +
    geom_pointrange(aes(ymin = Estimate - 2 * std_error, ymax = Estimate + 2 * std_error)) +
    # geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`)) + # same as above
    geom_text(aes(label = p_signif, y = Estimate + 3 * std_error), data = model_df[!is.na(p_signif)], colour = 'red') +
    theme_cowplot() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    coord_flip() +
    scale_colour_manual(values = colour_pal) +
    labs(title = "Predicting #SCNAs on a branch in\nTRACERx421 primary-metastasis cohort")
pdf(file.path(save_directory, "mets_lm_snv_scna_eventsnorm.pdf"), width = 9, height = 6)
print(g)
dev.off()
# Also save in figures directory
pdf("figures/Suppfig6a_mets_lm_snv_scna_eventsnorm.pdf", width = 9, height = 6)
print(g)
dev.off()

### Logistic regression predicting whether clone is metastatic clone 
# Predicting met clones: Variance explained by gains, losses and mutations, correcting for tumour id

### Specify and fit a generalized linear mixed model (GLMM) with tumour ID as random effect
mets_lm_full <- glmer(metClones ~ total_gains_norm_subclonal + total_losses_norm_subclonal + mut_branchlength_norm_subclonal + Histology + (1|tumour_id), data = mets_edge_df[is_mrca == 'subclone'], family = binomial)
summary(mets_lm_full)

model_df <- as.data.table(coef(summary(mets_lm_full)), keep.rownames = 'Term')
model_df <- setNames(model_df, c("Term", "Estimate", "std_error", "z_value", "p_value"))
model_df[p_value < 0.05, p_signif := "*"]
model_df[p_value < 0.01, p_signif := "**"]
model_df[p_value < 0.001, p_signif := "***"]

model_df[Term == 'mut_branchlength_norm_subclonal', Term := '#SNVs']
model_df[grepl("total_", Term), Term := gsub("total_", 'N. ', Term)]
model_df[grepl("_norm_subclonal", Term), Term := gsub("_norm_subclonal", '', Term)]


### Plot linear model results
colour_pal <- c('#SNVs' = brewer.pal(10, 'Paired')[1], 
                'N. gains' = brewer.pal(11, 'RdBu')[2], 
                'N. losses' = brewer.pal(11, 'RdBu')[10],
                'HistologyLUSC' = as.character(unlist(tx_palette$histology)['LUSC']), 
                'HistologyOther' = as.character(unlist(tx_palette$histology)['Other']),
                '(Intercept)' = "grey50")

g <- ggplot(model_df,
    aes(Term, Estimate, colour = Term)) +
    geom_pointrange(aes(ymin = Estimate - 2 * std_error, ymax = Estimate + 2 * std_error)) +
    # geom_pointrange(aes(ymin = `2.5 %`, ymax = `97.5 %`)) + # same as above
    geom_text(aes(label = p_signif, y = Estimate + 3 * std_error), data = model_df[!is.na(p_signif)], colour = 'red') +
    theme_cowplot() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    coord_flip() +
    scale_colour_manual(values = colour_pal) +
    labs(title = "Predicting metastatic\nclones in TRACERx421\nprimary-metastasis cohort")
pdf(file.path(save_directory, "mets_lm_snv_scna_predictingMetClones.pdf"), width = 7, height = 5)
print(g)
dev.off()
# Also save in figures directory
pdf("figures/Suppfig6c_mets_lm_snv_scna_predictingMetClones.pdf", width = 7, height = 5)
print(g)
dev.off()

###########
### END ###
###########