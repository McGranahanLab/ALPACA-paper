#!/bin/bash

# This script generates all figures in the ALPACA manuscript.
# Please ensure you run `cd ALPACA-paper` before running this script.
# Activate evnironment with: conda activate alpaca_conda
# Run this script using the following command:
# ./run_all_figures.sh

if [[ ! "${PWD##*/}" == "ALPACA-paper" ]]; then
  echo "Error: Script must be run from the 'ALPACA-paper' directory."
  exit 1
fi
mkdir -p figures
### FIGURE 1 ###

papermill bin/ANALYSIS/Fig1.ipynb bin/ANALYSIS/Fig1_executed.ipynb

# Supplementary Figure 1

papermill bin/ANALYSIS/Fig1.ipynb bin/ANALYSIS/Supplemenrary_fig1_executed.ipynb

### FIGURE 2 ###

# Figure 2a
Rscript bin/ANALYSIS/ADDRESOLUTION_plot_samples_clones.R \
--tumour_id CRUK0048 \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/example_cases

Rscript bin/ANALYSIS/ADDRESOLUTION_plot_trees_clonemaps.R \
--tumour_id CRUK0048 \
--save_directory output/mets/mets_default/cohort_outputs/example_cases

# Figure 2b
Rscript bin/ANALYSIS/ADDRESOLUTION_plot_samples_clones.R \
--tumour_id CRUK0022 \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/example_cases

Rscript bin/ANALYSIS/ADDRESOLUTION_plot_trees_clonemaps.R \
--tumour_id CRUK0022 \
--save_directory output/mets/mets_default/cohort_outputs/example_cases


# Supplementary Figure 2

papermill bin/ANALYSIS/Supplementary_fig2.ipynb bin/ANALYSIS/Supplementary_fig2_executed.ipynb

### FIGURE 3 ###

# Figure 3b.i, 3c.i
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUAD_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin

# Figure 3b.ii, 3c.ii
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUSC_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin

# Supplementary Figure 3

papermill bin/ANALYSIS/Supplementary_fig3.ipynb bin/ANALYSIS/Supplementary_fig3_executed.ipynb

# Supplementary Figure 4a

# LUAD
Rscript bin/ANALYSIS/DRIVERORDERING_supplementaryanalysis.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUAD_all_topfreqevents.csv \
--all_events_path output/primary/primary_default/cohort_outputs/eventordering/LUAD_all_allevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin \
--save_directory output/primary/primary_default/cohort_outputs/eventordering

# LUSC
Rscript bin/ANALYSIS/DRIVERORDERING_supplementaryanalysis.R \
--top_events_path output/primary/primary_default/cohort_outputs/eventordering/LUSC_all_topfreqevents.csv \
--all_events_path output/primary/primary_default/cohort_outputs/eventordering/LUSC_all_allevents.csv \
--cohort_tree_directory _assets/tree_data/primary \
--bin_directory bin \
--save_directory output/primary/primary_default/cohort_outputs/eventordering

# Supplementary Figure 4b.i, 4c.i
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/mets/mets_default/cohort_outputs/eventordering/LUAD_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/mets \
--bin_directory bin

# Supplementary Figure 4b.ii, 4c.ii
Rscript bin/ANALYSIS/DRIVERORDERING_run_btm.R \
--top_events_path output/mets/mets_default/cohort_outputs/eventordering/LUSC_all_topfreqevents.csv \
--cohort_tree_directory _assets/tree_data/mets \
--bin_directory bin

### FIGURE 4 ###

# Figure 4b, Supplementary Figure 5e
Rscript bin/ANALYSIS/COEVOLUTION_plot_mut_cpn_rate.R \
--clone_distances_path output/mets/mets_default/cohort_outputs/edge_events.csv \
--bin_directory bin \
--save_directory output/mets/mets_default/cohort_outputs/constellation

# Figure 4c, Supplementary Figures 5a, 5c
Rscript bin/ANALYSIS/COEVOLUTION_compare_branch_lengths_combinedcohorts.R \
--mets_clone_distances_path output/mets/mets_default/cohort_outputs/edge_events.csv \
--prim_clone_distances_path output/primary/primary_default/cohort_outputs/edge_events.csv \
--bin_directory bin \
--save_directory output/mets/mets_default/cohort_outputs/n_events

# Figure 4d
Rscript bin/ANALYSIS/SEEDING_vs_nonseeding_plot_cross_genome_changes.R \
--edge_cn_change_seg_chrarm_path output/mets/mets_default/cohort_outputs/edge_cn_change_seg_chrarm.csv \
--bin_directory bin \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/seeding_vs_nonseeding

# Figures 4e, 4k, 4l, Supplementary Figures 5b, 5d, 5h
Rscript bin/ANALYSIS/SEEDING_compare_seeding_non_seeding.R \
--clone_distances_path output/mets/mets_default/cohort_outputs/edge_events.csv \
--clone_metrics_path output/mets/mets_default/cohort_outputs/clone_ploidy_floh.csv \
--bin_directory bin \
--save_directory output/mets/mets_default/cohort_outputs/n_events

# Figures 4g, 4h, 4i, Supplementary Figures 5e, 5f, 5g
Rscript bin/ANALYSIS/SEEDING_compare_cloneclass_events.R \
--edge_cn_change_seg_chrarm_path output/mets/mets_default/cohort_outputs/edge_cn_change_seg_chrarm.csv \
--n_tums_event_thresh 5 \
--bin_directory bin \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/seeding_vs_nonseeding

### FIGURE 5 ###

# Figures 4j, 5b, 5c, 5d, 5e, 5f, 5g, Supplementary Figures 6a, 6b, 6c, 6d, 6e, 6f

papermill bin/ANALYSIS/Fig5b.ipynb bin/ANALYSIS/Fig5b_executed.ipynb

Rscript bin/ANALYSIS/SURVIVAL_scnaith_diversity_pipeline.R \
--clone_distances_path output/primary/primary_default/cohort_outputs/edge_events.csv \
--pairwise_dist_path output/primary/primary_default/cohort_outputs/cn_dist_pairwise.csv \
--bin_directory bin \
--mets_output_directory output/mets/mets_default \
--save_directory output/primary/primary_default/cohort_outputs/scna_ith

### END ###
