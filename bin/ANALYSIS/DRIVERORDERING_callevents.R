#!/usr/bin/env Rscript

# To run this script for both cohorts 'primary' and 'mets',
# run the following command from the command line (make sure you are in the root directory of this repo): 
# Primary cohort
"
Rscript bin/ANALYSIS/DRIVERORDERING_callevents.R \
--edge_cn_change_seg_chrarm_path output/primary/primary_default/cohort_outputs/edge_cn_change_seg_chrarm.csv \
--n_chrarm_drivclass_events 10 \
--n_sc_chrarm_drivclass_events 5 \
--n_focal_drivclass_events 10 \
--n_sc_focal_drivclass_events 5 \
--n_snv_driver_events 10 \
--n_sc_snv_driver_events 10 \
--min_n_tums_w_event 5 \
--bin_directory bin \
--mut_table_path _assets/mutTable/primary_mutTable.csv \
--output_directory output/primary/primary_default \
--save_directory output/primary/primary_default/cohort_outputs/eventordering
"

# Mets cohort
"
Rscript bin/ANALYSIS/DRIVERORDERING_callevents.R \
--edge_cn_change_seg_chrarm_path output/mets/mets_default/cohort_outputs/edge_cn_change_seg_chrarm.csv \
--n_chrarm_drivclass_events 10 \
--n_sc_chrarm_drivclass_events 5 \
--n_focal_drivclass_events 10 \
--n_sc_focal_drivclass_events 5 \
--n_snv_driver_events 10 \
--n_sc_snv_driver_events 10 \
--min_n_tums_w_event 5 \
--bin_directory bin \
--mut_table_path _assets/mutTable/mets_mutTable.csv \
--output_directory output/mets/mets_default \
--save_directory output/mets/mets_default/cohort_outputs/eventordering
"
cat("\nSelecting top most frequent driver events in the cohort\n")

options(stringsAsFactors = FALSE)

### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

### Command line options 
option_list <- list(
  make_option(c("--edge_cn_change_seg_chrarm_path"), type="character", action="store", help="File path to combined cohort edge copy number change file with chrarm added info", metavar="character"),
  make_option(c("--n_chrarm_drivclass_events"), type="numeric", action="store", help="No. top chr arm events to use in Bradley Terry model (amps and LOHs)", metavar="numeric"),
  make_option(c("--n_sc_chrarm_drivclass_events"), type="numeric", action="store", help="No. top subclonal SCNA arm events to use in Bradley Terry model", metavar="numeric"),
  make_option(c("--n_focal_drivclass_events"), type="numeric", action="store", help="No. top focal SCNA events to use in Bradley Terry model (amps and LOHs)", metavar="numeric"),
  make_option(c("--n_sc_focal_drivclass_events"), type="numeric", action="store", help="No. top SUBCLONAL focal SCNA events to use in Bradley Terry model (amps and LOHs)", metavar="numeric"),
  make_option(c("--n_snv_driver_events"), type="numeric", action="store", help="No. top SNV driver events to use in Bradley Terry model", metavar="numeric"),
  make_option(c("--n_sc_snv_driver_events"), type="numeric", action="store", help="No. top SUBCLONAL SNV driver events to use in Bradley Terry model", metavar="numeric"), 
  make_option(c("--min_n_tums_w_event"), type="numeric", action="store", help="Minimum no. tumours with event required to assign event as one of the top most frequent", metavar="numeric"), 
  make_option(c("--bin_directory"), type="character", help="Directory where R scripts are saved"),
  make_option(c("--mut_table_path"), type="character", help="mut_table_path"),
  make_option(c("--output_directory"), type="character", help="output_directory"),
  make_option(c("--save_directory"), type="character", help="save_directory")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)

### Parse args
edge_cn_change_seg_chrarm_path <- opt$edge_cn_change_seg_chrarm_path
n_chrarm_drivclass_events <- opt$n_chrarm_drivclass_events
n_sc_chrarm_drivclass_events <- opt$n_sc_chrarm_drivclass_events
n_focal_drivclass_events <- opt$n_focal_drivclass_events
n_sc_focal_drivclass_events <- opt$n_sc_focal_drivclass_events
n_snv_driver_events <- opt$n_snv_driver_events 
n_sc_snv_driver_events <- opt$n_sc_snv_driver_events 
min_n_tums_w_event <- opt$min_n_tums_w_event
bin_dir <- opt$bin_directory
mut_table_path <- opt$mut_table_path
output_directory <- opt$output_directory
save_directory <- opt$save_directory

assets_dir <- gsub('bin', '_assets', bin_dir)
patient_dir <- file.path(output_directory, "patient_outputs")
print(patient_dir)

if (!file.exists(save_directory)) dir.create(save_directory, recursive = TRUE)

### Functions
source(file.path(bin_dir, "ANALYSIS/analysis_functions.R"))
source(file.path(bin_dir, "ANALYSIS/DRIVERORDERING_functions.R"))


### Load and process data
# Clinical data
clin_data_path <- file.path(assets_dir, "clinical_data_tumour.tsv")
clin_small <- load_reduced_clinicaldata(clin_data_path)

# Gene info
gene_loc_path <- file.path(assets_dir, "hg19refgene_tx421driveranno.csv")
gene_anno <- load_gene_anno(gene_loc_path)

# SNV driver muttable
muttable <- fread(mut_table_path)
gene_drivclass <- unique(gene_anno[is_pancan_driver == TRUE, .(gene_name, oncogene, tumour_suppressor, cn_driver_class, mut_driver_class, is_lung_mut_driver, is_pancan_driver)])
driver_muttab_small <- process_driver_muttable_for_btm(muttable, gene_drivclass)
driver_muttab_small <- merge.data.table(driver_muttab_small, clin_small, by = "tumour_id", all.x = TRUE, all.y = FALSE)
driver_muttab_small[, clone := as.character(cluster_with_driver)]
# only including lung mut drivers:
driver_muttab_small <- driver_muttab_small[is_lung_mut_driver == TRUE]

# ALPACA edge CN change data
edge_change_arm_df <- fread(edge_cn_change_seg_chrarm_path)
edge_change_arm_df[, clone := as.character(clone)]
edge_change_arm_df <- merge.data.table(edge_change_arm_df, clin_small, by = "tumour_id", all.x = TRUE, all.y = FALSE)

############
### MAIN ###
############

### Classify CN events per tumour 
# Get tumour level event table for each edge, including:
#   - Driver SNVs
#   - Amplifications affecting oncogenes
#   - Amplifications and LOH affecting chr arms
#   - if arm_amp == TRUE -> gene_amp = FALSE

### Merge genes back to CN change per edge, by gene-tumour segment mapping
cat('\nMerge genes back to CN change per edge, by gene-tumour segment mapping\n')
all_tumours <- unique(edge_change_arm_df$tumour_id)
gene_seg_mapping_df <- rbindlist(lapply(all_tumours, function(tum_id) {
  print(tum_id)
  tum_seg_gene <- fread(file.path(patient_dir, tum_id, "tum_seg_gene_matched.csv"))
  drivergene_seg_mapping_tum <- tum_seg_gene[is_lung_mut_driver == TRUE]
  return(drivergene_seg_mapping_tum)
}))
gene_edge_cn <- merge.data.table(gene_seg_mapping_df, edge_change_arm_df, by = c("tumour_id", "segment"), all.x = TRUE, all.y = FALSE, allow.cartesian = TRUE)

### Classify events in driver genes specifically 
cat('\nClassify events in driver genes specifically\n')
# Only consider an oncogene as amplified if the gene is inside a CN segment
# Consider an TSG as lost if the gene overlaps with a lost CN segment
gene_edge_cn[, onco_amps := overlap_type == "query_contained_in_subject" & cn_driver_class == "cn_oncogene" & any_amp == TRUE & (is_lung_mut_driver == TRUE | is_cn_driver == TRUE)]
gene_edge_cn[, tsg_homodel := cn_driver_class == "cn_tsg" & homodel == TRUE & (is_lung_mut_driver == TRUE | is_cn_driver == TRUE)]
gene_edge_cn[, tsg_loh := cn_driver_class == "cn_tsg" & any_loh == TRUE & (is_lung_mut_driver == TRUE | is_cn_driver == TRUE)]

### Classify amplifications on tree edges:
cat('\nClassify amplifications on tree edges\n')
# - if arm_amp == TRUE -> gene_amp = FALSE
gene_edge_cn[, onco_amps_final := onco_amps]
gene_edge_cn[any_arm_amp_acquired_edge == TRUE, onco_amps_final := FALSE]


### Filter drivers based on prevalence across the cohort 
cat('\nFilter drivers based on prevalence across the cohort\n')
# Get top most frequent events, including most frequent subclonal
# Make sure to remove events from BT model if there are not enough tumours with those events 
# (e.g. only use one homozygous deletion in LUSC)

histologies <- c("LUAD", "LUSC")

for (histo in histologies) {
  cat('\n\n', histo, '\n')

  ### Combined truncal and subclonal events
  print("Combined truncal and subclonal events")
  clonality <- "all"
  
  all_edges <- unique(edge_change_arm_df[Histology == histo, .(tumour_id, clone, parent_clone)])
  setorder(all_edges, tumour_id, clone)

  # Extract distinct dataframes for each event
  edge_chrarm_amp_df <- unique(edge_change_arm_df[Histology == histo & any_arm_amp_acquired_edge == TRUE, .(tumour_id, clone, event_id = paste0("+", chr_arm))])
  edge_chrarm_loh_df <- unique(edge_change_arm_df[Histology == histo & any_arm_loh_acquired_edge == TRUE, .(tumour_id, clone, event_id = paste(chr_arm, "LOH"))])
  edge_onco_amp_df <- gene_edge_cn[Histology == histo & onco_amps_final == TRUE, .(tumour_id, clone, event_id = paste0("+", gene_name))]
  driver_muttab <- driver_muttab_small[Histology == histo, .(tumour_id, clone, event_id = gene_name)]
  edge_chrarm_amp_sc_df <- unique(edge_change_arm_df[Histology == histo & any_arm_amp_acquired_edge == TRUE & parent_clone != "diploid", .(tumour_id, clone, event_id = paste0("+", chr_arm))])
  edge_chrarm_loh_sc_df <- unique(edge_change_arm_df[Histology == histo & any_arm_loh_acquired_edge == TRUE & parent_clone != "diploid", .(tumour_id, clone, event_id = paste(chr_arm, "LOH"))])

  # Save in a list for input into filtering function
  selected_events_list <- list(
    "chr arm amp" = list(name = "chr arm amp", n_events = n_chrarm_drivclass_events, df = edge_chrarm_amp_df, full_df = edge_chrarm_amp_df),
    "chr arm LOH" = list(name = "chr arm LOH", n_events = n_chrarm_drivclass_events, df = edge_chrarm_loh_df, full_df = edge_chrarm_loh_df),
    "gene amp" = list(name = "gene amp", n_events = n_focal_drivclass_events, df = edge_onco_amp_df, full_df = edge_onco_amp_df),
    "gene snv" = list(name = "gene snv", n_events = n_snv_driver_events, df = driver_muttab, full_df = driver_muttab),
    "chr arm amp sc" = list(name = "chr arm amp sc", n_events = n_sc_chrarm_drivclass_events, df = edge_chrarm_amp_sc_df, full_df = edge_chrarm_amp_df),
    "chr arm LOH sc" = list(name = "chr arm LOH sc", n_events = n_sc_chrarm_drivclass_events, df = edge_chrarm_loh_sc_df, full_df = edge_chrarm_loh_df)
  )

  # Quantify event frequency
  histo_all_events <- get_event_frequency_fromlist(selected_events_list = selected_events_list, min_n_tums_w_event = min_n_tums_w_event)
  
  # Extract top most frequent events, and extract all tumours harbouring any clone with those top events present.
  histo_top_events <- unique(histo_all_events[top_event == TRUE, .(tumour_id, clone, event_id)])
  tumour_events <- merge.data.table(histo_top_events, all_edges, by = c("tumour_id", "clone"), all = TRUE)
  tumour_events <- tumour_events[!is.na(parent_clone)]
  
  # Save
  write.csv(histo_all_events, paste0(save_directory, '/', paste(histo, clonality, "allevents.csv", sep = "_")), row.names = FALSE, quote = FALSE)
  write.csv(tumour_events, paste0(save_directory, '/', paste( histo, clonality, "topfreqevents.csv", sep = "_")), row.names = FALSE, quote = FALSE)

  ### Subclonal events only
  print("Subclonal events only")
  clonality <- "subclonal"

  # Extract distinct dataframes for each event (now create the subclonal only dataframes we haven't yet filtered)
  edge_onco_amp_sc_df <- gene_edge_cn[Histology == histo & onco_amps_final == TRUE & parent_clone != "diploid", .(tumour_id, clone, event_id = paste0("+", gene_name))]
  driver_muttab_sc <- driver_muttab_small[Histology == histo & cluster_with_driver_clonality == "S", .(tumour_id, clone, event_id = gene_name)]

  # Save in a list for input into filtering function
  selected_events_list_sc <- list(
    "gene amp sc" = list(name = "gene amp", n_events = n_sc_focal_drivclass_events, df = edge_onco_amp_sc_df, full_df = edge_onco_amp_sc_df),
    "gene snv sc" = list(name = "gene snv", n_events = n_sc_snv_driver_events, df = driver_muttab_sc, full_df = driver_muttab_sc),
    "chr arm amp sc" = list(name = "chr arm amp sc", n_events = n_sc_chrarm_drivclass_events, df = edge_chrarm_amp_sc_df, full_df = edge_chrarm_amp_sc_df),
    "chr arm LOH sc" = list(name = "chr arm LOH sc", n_events = n_sc_chrarm_drivclass_events, df = edge_chrarm_loh_sc_df, full_df = edge_chrarm_loh_sc_df)
  )

  # Quantify event frequency
  histo_all_events_sc <- get_event_frequency_fromlist(selected_events_list = selected_events_list_sc, min_n_tums_w_event = min_n_tums_w_event)
  
  # Extract top most frequent events, and extract all tumours harbouring any clone with those top events present.
  histo_top_events_sc <- unique(histo_all_events_sc[top_event == TRUE, .(tumour_id, clone, event_id)])
  tumour_events_sc <- merge.data.table(histo_top_events_sc, all_edges, by = c("tumour_id", "clone"), all = TRUE)
  tumour_events_sc <- tumour_events_sc[!is.na(parent_clone)]
  
  # Save
  write.csv(histo_all_events_sc, paste0(save_directory, '/', paste(histo, clonality, "allevents.csv", sep = "_")), row.names = FALSE, quote = FALSE)
  write.csv(tumour_events_sc, paste0(save_directory, '/', paste(histo, clonality, "topfreqevents.csv", sep = "_")), row.names = FALSE, quote = FALSE)
}

