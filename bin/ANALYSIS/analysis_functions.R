### R functions for SCNA evolution analysis from ALPACA output data
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
options(stringsAsFactors = FALSE)

### Function to read in tree paths for one tumour
read_tree_json <- function(tree_json_path) {
  suppressPackageStartupMessages(require(jsonlite))
  tree_paths <- read_json(path = tree_json_path)
  return(tree_paths)
}

### Function to load genomic loci of chromosome arms (currently using hg19)
load_chr_arm_loci <- function(cytoband_loc_path) {
  chr_info <- fread(cytoband_loc_path)
  chr_info[, chromosome := gsub("chr", "", chrom)]
  chr_info[, chr_arm := substr(band, 1, 1)]
  chr_info[, chr_arm := paste0(chromosome, chr_arm)]
  chr_info[, arm_start := min(start), by = chr_arm]
  chr_info[, arm_end := max(end), by = chr_arm]
  chr_info[, arm_length := arm_end - arm_start]
  chr_arms <- unique(chr_info[chromosome %in% 1:22, 
      .(chromosome, arm_start, arm_end, chr_arm, arm_length)]
  )
  return(chr_arms)
}


### Function to load gene annotation file
load_gene_anno <- function(gene_loc_path) {
  require(data.table)
  print("Reading gene annotation file")
  gene_anno_df <- fread(gene_loc_path)
  
  ### Preprocess gene annotation file:
  # if genes have >1 chromosome assigned,
  # exclude these from analysis (i.e. don"t include in background model)
  gene_anno_df[, n := length(unique(chromosome)), by = gene_name]
  gene_anno <- gene_anno_df[n == 1]
  gene_anno$n <- NULL
  
  # summarise multiple gene positions to widest interval per gene:
  gene_anno[, chr := as.numeric(chromosome)]
  gene_anno[, start := min(as.numeric(start)), by = gene_name]
  gene_anno[, end := max(as.numeric(end)), by = gene_name]
  gene_anno[, chromosome := NULL]
  
  # add cytoband & transcripts annotation (if some genes span multiple cytobands)
  gene_anno[, cytoband_anno := paste(sort(cytoband_name), collapse = ";"), by = gene_name]
  gene_anno[, cytoband_name := NULL]
  gene_anno[, transcripts_anno := paste(sort(transcripts), collapse = ";"), by = gene_name]
  gene_anno[, transcripts := NULL]
  # remove strand information
  gene_anno[, strand := NULL]
  gene_anno <- unique(gene_anno)
  
  n_genes <- length(unique(gene_anno[, gene_name]))
  print(paste0("Number of genes with unique chr: ", n_genes))

  # note, sometimes > 1 gene has the same genomic position - take the first gene name listed in these cases
  gene_anno[, n := length(unique(gene_name)), by = c("chr", "start", "end")]
  gene_anno[, gene_label := gene_name[1], by = c("chr", "start", "end")]
  gene_anno <- gene_anno[gene_label == gene_name]
  gene_anno[, gene_label := NULL]
  gene_anno[, n := NULL]
  gene_anno <- unique(gene_anno)

  n_genes <- length(unique(gene_anno[, gene_name]))
  print(paste0("Number of genes at a unique genomic locus: ", n_genes))    

  ### Additionally, obtain gene driver classifications
  #    - Define oncogenes and TSGs for SNVs and SCNAs
  driver_class_dt <- gene_anno[is_pancan_mut_driver == TRUE | is_cn_driver == TRUE]
  driver_class_dt <- driver_class_dt[cn_driver_OG_TSG %in% c("OG", "TSG") | mut_driver_OG_TSG %in% c("OG", "TSG")]

  # Add flag for gene annotation as oncogene or tsg in either SNV or SCNA space
  driver_class_dt[, oncogene := FALSE]
  driver_class_dt[cn_driver_OG_TSG == "OG" | mut_driver_OG_TSG == "OG", oncogene := TRUE]
  driver_class_dt[, tumour_suppressor := FALSE]
  driver_class_dt[cn_driver_OG_TSG == "TSG" | mut_driver_OG_TSG == "TSG", tumour_suppressor := TRUE]

  ### Copy number drivers
  # note: NFE2L2, TERT, MET: copy number oncogene and mutational TSG --> classing them as copy number oncogenes
  # define final SCNA driver class
  driver_class_dt[oncogene == TRUE, cn_driver_class := "cn_oncogene"]
  driver_class_dt[tumour_suppressor == TRUE, cn_driver_class := "cn_tsg"]
  driver_class_dt[oncogene == TRUE & tumour_suppressor == TRUE, cn_driver_class := "cn_oncogene"]

  ### SNV drivers
  driver_class_dt[is_pancan_mut_driver == TRUE & mut_driver_OG_TSG == "TSG", mut_driver_class := "mut_tsg"]
  driver_class_dt[is_pancan_mut_driver == TRUE & mut_driver_OG_TSG == "OG", mut_driver_class := "mut_oncogene"]

  driver_class_dt <- unique(driver_class_dt[,
      .(gene_name, chr, start, end, cytoband_anno, is_pancan_mut_driver, is_lung_mut_driver, is_cn_driver,
      cn_driver_OG_TSG, mut_driver_OG_TSG, oncogene, tumour_suppressor, cn_driver_class, mut_driver_class
  )])

  # merge back in with full gene annotation df
  driver_class_dt[, is_pancan_driver := TRUE]
  col_overlap <- colnames(driver_class_dt)[colnames(driver_class_dt) %in% colnames(gene_anno)]
  gene_anno_out <- merge.data.table(gene_anno, driver_class_dt, by = col_overlap, all = TRUE)
  gene_anno_out[is.na(is_pancan_driver), is_pancan_driver := FALSE]
  setorder(gene_anno_out, chr, start)

  return(gene_anno_out)
}

### Function to load clinical data
load_reduced_clinicaldata <- function(clin_data_path) {
  clin_df <- fread(clin_data_path)
  clin_small <- unique(clin_df[, .(tumour_id = tumour_id_muttable_cruk,
                                  Histology = Histology_per_tumour_id_muttable,
                                  # LUAD_pred_subtype_with.IMA_per_tumour,
                                  smoking_status_merged,
                                  pathologyTNM
                                  )])

  clin_small[, Histology := ifelse(Histology == "Invasive adenocarcinoma", "LUAD",
      ifelse(Histology == "Squamous cell carcinoma", "LUSC", "Other"))
  ]
  return(clin_small)
}


### Function to compute clone ploidy and FLOH
get_clone_ploidy_floh <- function(alpaca_tum_df) {
  alpaca_tum_df[, clone := gsub("clone", "", clone)]
  
  alpaca_tum_clone <- unique(alpaca_tum_df[, .(segment, clone, pred_CN_A, pred_CN_B)])
  alpaca_tum_clone[, start := sapply(segment, function(x) as.numeric(strsplit(x, "_")[[1]][2]))]
  alpaca_tum_clone[, end := sapply(segment, function(x) as.numeric(strsplit(x, "_")[[1]][3]))]
  alpaca_tum_clone[, seg_len := end - start]

  alpaca_tum_clone[, total_pred_CN := pred_CN_A + pred_CN_B]

  genome_len <- sum(unique(alpaca_tum_clone[, .(segment, seg_len)])$seg_len)
  alpaca_tum_clone[, total_pred_CN_w := total_pred_CN * seg_len / genome_len]
  alpaca_tum_clone[, seg_fraction := seg_len / genome_len]

  ### clone ploidy
  alpaca_tum_clone[, clone_ploidy := sum(total_pred_CN_w), by = clone]

  ### clone floh
  alpaca_tum_clone[, is_loh := pred_CN_A == 0 | pred_CN_B == 0]
  alpaca_tum_clone[, clone_floh := sum((is_loh == TRUE) * seg_fraction), by = clone]

  alpaca_tum_clone <- unique(alpaca_tum_clone[, .(clone,
                                                  clone_ploidy,
                                                  clone_floh)])
  return(alpaca_tum_clone)
}

### Function to convert tree paths to a tree edge matrix
convert_tree_paths_to_edgemat <- function(tree_paths) {
  tree_edges <- lapply(tree_paths, function(p) {
    p <- gsub("clone", "", p)
    path_edges <- cbind(p[-length(p)], p[-1])
    return(path_edges)
  })
  tree_edges <- unique(do.call(rbind, tree_edges))
  return(tree_edges)
}

### Function that extracts all end-to-end tree paths from an edge matrix
convert_edgemat_to_tree_paths <- function(tree_graph) {
  clones_in_tree = unique(c(as.matrix(tree_graph)))
  if (length(clones_in_tree)==1) paths = list(1)
  else {
    # tree_graph = apply(tree_graph, c(1,2), as.character)
    tree_graph = as.data.frame(tree_graph)
    colnames(tree_graph) = c("Parent", "Child")
    paths = list()
    trunk = unique(tree_graph$Parent[!tree_graph$Parent %in% tree_graph$Child])
    terminal.clones = tree_graph[!(tree_graph$Child %in% tree_graph$Parent), "Child"]
    for (terminal in terminal.clones) {
      p = c(terminal)
      current_clone = terminal
      while (current_clone != trunk) {
        parent = tree_graph[tree_graph$Child==current_clone, "Parent"]
        p = c(p,parent)
        current_clone = parent
      }
      p = rev(p)
      paths = append(paths, list(p))
    }
    return(paths)
  }
}

### Function to load clones with recent subclonal expansion
get_sc_exp_metrics <- function(tree_structure, tree_object, clone_props) {
  # recent subclonal expansion
  clonality_table <- tree_object$clonality_out$clonality_table_corrected
  clonality_table <- clonality_table[rownames(clonality_table) %in% c(tree_structure), , drop = FALSE]
  is_terminal <- sapply(rownames(clonality_table), function(x) x %in% tree_structure[, 2] & !(x %in% tree_structure[, 1]))
  recent_sc_exp_clones <- (clonality_table == "clonal") * (is_terminal)
  recent_sc_exp_clones <- apply(recent_sc_exp_clones, 1, function(x) any(x == 1))
  # clonal illusion
  clonalillu_clones <- apply(clonality_table, 1, function(x) any(x == "clonal"))
  # max clone size in any region
  clone_props <- as.data.frame(clone_props)
  rownames(clone_props) <- gsub("clone", "", clone_props$clone)
  clone_props$clone <- NULL
  max_clone_size <- apply(clone_props, 1, max)
  # max clone size in any primary region
  primary_regs <- colnames(clone_props)[grepl('\\.R',  colnames(clone_props))]
  max_prim_clone_size <- apply(clone_props[, primary_regs, drop = FALSE], 1, max)
  # max CCF in any region
  ccf_cluster_table <- tree_object$nested_pyclone$ccf_cluster_table
  ccf_cluster_table <- ccf_cluster_table[rownames(ccf_cluster_table) %in% c(tree_structure), , drop = FALSE]
  max_clone_ccf <- apply(ccf_cluster_table, 1, max)
  # max CCF in any primary region 
  primary_regs_ccf <- colnames(ccf_cluster_table)[grepl('\\.R',  colnames(ccf_cluster_table))]
  max_prim_clone_ccf <- apply(ccf_cluster_table[, primary_regs_ccf, drop = FALSE], 1, max)  
  # out_df <- data.table(clone = names(recent_sc_exp_clones), is_recent_sc_exp = recent_sc_exp_clones)
  
  return(data.table(clone = names(recent_sc_exp_clones),
                    is_terminal = is_terminal,
                    is_recent_sc_exp = recent_sc_exp_clones,
                    is_clonal_illusion = clonalillu_clones,
                    max_clone_size = max_clone_size,
                    max_prim_clone_size = max_prim_clone_size,
                    max_clone_ccf = max_clone_ccf,
                    max_prim_clone_ccf = max_prim_clone_ccf,
                    max_terminal_ccf_tum = max(is_terminal * max_clone_size)))
}

### Function to convert distance matrix to pairwise dataframe
dist_to_df <- function(dist_matrix) {
  dist_matrix <- as.matrix(dist_matrix)
  col_pairs <- t(combn(colnames(dist_matrix), 2))
  pairs_df <- data.frame(col_pairs, dist = dist_matrix[col_pairs])
  return(pairs_df)
}

### Function to extract cn distance between normal diploid clone and MRCA
get_diploid_mrca_cn_distance <- function(alpaca_mrca, mrca_clone, dist_method = "euclidean") {
  # get cn distance from normal diploid --> mrca clone:
  require(dplyr)

  # alpaca_mrca <- alpaca_clone[alpaca_clone$cluster == mrca_clone]
  alpaca_mrca[, pred_CN_A_w := pred_CN_A * seg_fraction]
  alpaca_mrca[, pred_CN_B_w := pred_CN_B * seg_fraction]

  # The weighted diploid CN is just equal to the segment fraction of the
  # genome for each allele

  # alpaca_mrca[, diploid_CN_A_w := 1 * seg_fraction]
  # alpaca_mrca[, diploid_CN_B_w := 1 * seg_fraction]

  mrca_cn_dist_A <- as.numeric(dist(t(alpaca_mrca[, .(pred_CN_A_w, seg_fraction)])))
  mrca_cn_dist_B <- as.numeric(dist(t(alpaca_mrca[, .(pred_CN_B_w, seg_fraction)])))

  distance_df <- data.frame(clone1 = "diploid",
                            clone2 = unique(alpaca_mrca$clone),
                            cn_distance_A = mrca_cn_dist_A,
                            cn_distance_B = mrca_cn_dist_B,
                            mean_dist_allele = mean(c(mrca_cn_dist_A, mrca_cn_dist_B)))
  return(distance_df)
}

### Function to extract pairwise cn distance
pairwise_cn_distance_alpaca <- function(alpaca_clone, large_seg_thresh = NULL, small_seg_thresh = NULL) {

  alpaca_clone_copy <- copy(alpaca_clone)
  alpaca_clone_copy[, chr := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][1]), segment)]
  alpaca_clone_copy[, start := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][2]), segment)]
  alpaca_clone_copy[, end := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][3]), segment)]
  alpaca_clone_copy[, seg_len := end - start]
  
  # Filter ALPACA output by only large/only small segments
  if (!is.null(large_seg_thresh) & !is.null(small_seg_thresh)) {
    stop("Please only filter output by large OR small segments, not both")
  }
  if (!is.null(large_seg_thresh)) {
    alpaca_clone_largesegs <- alpaca_clone_copy[seg_len > large_seg_thresh]
    if (nrow(alpaca_clone_largesegs) > 0) {
      alpaca_clone_copy <- copy(alpaca_clone_largesegs)
    } else return(NULL)
  }
  if (!is.null(small_seg_thresh)) {
    alpaca_clone_smallsegs <- alpaca_clone_copy[seg_len < small_seg_thresh]
    if (nrow(alpaca_clone_smallsegs) > 0) {
      alpaca_clone_copy <- copy(alpaca_clone_smallsegs)
    } else return(NULL)
  }

  alpaca_clone_copy[, genome_len := sum(seg_len), by = "clone"]
  alpaca_clone_copy[, seg_fraction := seg_len / genome_len]
  # Weighted copy number
  alpaca_clone_copy[, pred_CN_A_w := pred_CN_A * seg_fraction]
  alpaca_clone_copy[, pred_CN_B_w := pred_CN_B * seg_fraction]
  setorder(alpaca_clone_copy, chr, start)

  # widen df to get new column for each clone
  alpaca_A <- dcast(alpaca_clone_copy, chr + start + segment ~ clone, value.var = "pred_CN_A")
  alpaca_B <- dcast(alpaca_clone_copy, chr + start + segment ~ clone, value.var = "pred_CN_B")
  alpaca_A_w <- dcast(alpaca_clone_copy, chr + start + segment ~ clone, value.var = "pred_CN_A_w")
  alpaca_B_w <- dcast(alpaca_clone_copy, chr + start + segment ~ clone, value.var = "pred_CN_B_w")
  setorder(alpaca_A, chr, start)
  setorder(alpaca_B, chr, start)
  setorder(alpaca_A_w, chr, start)
  setorder(alpaca_B_w, chr, start)
  
  alpaca_AB <- rbind(alpaca_A, alpaca_B)
  
  # Compute cn distance between each pair of clones:
  cols_to_drop <- c("chr", "start", "segment")
  cols_to_keep <- colnames(alpaca_A)[!colnames(alpaca_A) %in% cols_to_drop]
  alpaca_A_mat <- apply(as.matrix(alpaca_A[, ..cols_to_keep, drop = FALSE]), c(1, 2), as.numeric)
  alpaca_B_mat <- apply(as.matrix(alpaca_B[, ..cols_to_keep, drop = FALSE]), c(1, 2), as.numeric)
  alpaca_A_w_mat <- apply(as.matrix(alpaca_A_w[, ..cols_to_keep, drop = FALSE]), c(1, 2), as.numeric)
  alpaca_B_w_mat <- apply(as.matrix(alpaca_B_w[, ..cols_to_keep, drop = FALSE]), c(1, 2), as.numeric)
  alpaca_AB_mat <- apply(as.matrix(alpaca_AB[, ..cols_to_keep, drop = FALSE]), c(1, 2), as.numeric)


  l1_clone_cn_dist_A <- dist_to_df(dist(t(alpaca_A_mat), method = "manhattan"))
  l1_clone_cn_dist_B <- dist_to_df(dist(t(alpaca_B_mat), method = "manhattan"))
  l2_clone_cn_dist_A <- dist_to_df(dist(t(alpaca_A_mat), method = "euclidean"))
  l2_clone_cn_dist_B <- dist_to_df(dist(t(alpaca_B_mat), method = "euclidean"))
  l2_clone_cn_dist_AB <- dist_to_df(dist(t(alpaca_AB_mat), method = "euclidean"))
  l1_clone_cn_dist_A_w <- dist_to_df(dist(t(alpaca_A_w_mat), method = "manhattan"))
  l1_clone_cn_dist_B_w <- dist_to_df(dist(t(alpaca_B_w_mat), method = "manhattan"))
  l2_clone_cn_dist_A_w <- dist_to_df(dist(t(alpaca_A_w_mat), method = "euclidean"))
  l2_clone_cn_dist_B_w <- dist_to_df(dist(t(alpaca_B_w_mat), method = "euclidean"))
  

  # Set colnames
  colnames(l1_clone_cn_dist_A) <- c("clone1", "clone2", "l1_cn_distance_A")
  colnames(l1_clone_cn_dist_B) <- c("clone1", "clone2", "l1_cn_distance_B")
  colnames(l2_clone_cn_dist_A) <- c("clone1", "clone2", "l2_cn_distance_A")
  colnames(l2_clone_cn_dist_B) <- c("clone1", "clone2", "l2_cn_distance_B")
  colnames(l2_clone_cn_dist_AB) <- c("clone1", "clone2", "l2_cn_distance_AB")
  colnames(l1_clone_cn_dist_A_w) <- c("clone1", "clone2", "l1_cn_distance_A_w")
  colnames(l1_clone_cn_dist_B_w) <- c("clone1", "clone2", "l1_cn_distance_B_w")
  colnames(l2_clone_cn_dist_A_w) <- c("clone1", "clone2", "l2_cn_distance_A_w")
  colnames(l2_clone_cn_dist_B_w) <- c("clone1", "clone2", "l2_cn_distance_B_w")

  distance_df <- as.data.table(Reduce(function(...) merge(..., by = c("clone1", "clone2"), all = TRUE),
    list(l1_clone_cn_dist_A, l1_clone_cn_dist_B, l2_clone_cn_dist_A, l2_clone_cn_dist_B,
          l1_clone_cn_dist_A_w, l1_clone_cn_dist_B_w, l2_clone_cn_dist_A_w, l2_clone_cn_dist_B_w,
          l2_clone_cn_dist_AB)))

  distance_df[, l1_mean_dist_allele := (l1_cn_distance_A + l1_cn_distance_B) / 2]
  distance_df[, l2_mean_dist_allele := (l2_cn_distance_A + l2_cn_distance_B) / 2]
  distance_df[, l1_mean_dist_allele_w := (l1_cn_distance_A_w + l1_cn_distance_B_w) / 2]
  distance_df[, l2_mean_dist_allele_w := (l2_cn_distance_A_w + l2_cn_distance_B_w) / 2]

  return(as.data.table(distance_df))
}

### Function to compute interval events between a pair of clones based on
# their cross genome copy number profile
get_interval_events <- function(parent_cn, child_cn, change_threshold_cap = NA) {
  cn_change <- child_cn - parent_cn

  if (!is.na(change_threshold_cap)) {
    cn_change[cn_change > change_threshold_cap] <- change_threshold_cap
    cn_change[cn_change < -change_threshold_cap] <- -change_threshold_cap
  }

  # count interval events
  interval_events <- c(cn_change[1])
  
  if (length(cn_change) > 1) {
    for (seg_change_id in 2:length(cn_change)) {
      seg_change <- cn_change[seg_change_id]
      if (interval_events[length(interval_events)] != seg_change) {
        interval_events <- c(interval_events, seg_change)
      }
    }
  }

  gains <- interval_events[interval_events > 0]
  losses <- interval_events[interval_events < 0] * -1

  # add binary events:
  gains_bin <- as.numeric(gains > 0)
  losses_bin <- as.numeric(losses > 0)

  interval_events_df <- data.frame(interval_gains = sum(gains),
                                   interval_losses = sum(losses),
                                   interval_gains_bin = sum(gains_bin),
                                   interval_losses_bin = sum(losses_bin))

  interval_events_df$interval_events = interval_events_df$interval_gains + interval_events_df$interval_losses
  interval_events_df$interval_events_bin = interval_events_df$interval_gains_bin + interval_events_df$interval_losses_bin

  return(interval_events_df)
}

### Function to compute interval events
compute_interval_events_alpaca <- function(edge_df, alpaca_tum_out) {  
  alpaca_clone_cn <- unique(alpaca_tum_out[, .(segment, clone = gsub('clone', '', clone), pred_CN_A, pred_CN_B)])
  alpaca_clone_cn[, chr := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][1]), segment)]
  alpaca_clone_cn[, start := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][2]), segment)]

  setorder(alpaca_clone_cn, chr, start)
  all_clones <- unique(edge_df[, clone])
  chromosomes <- sort(unique(as.numeric(alpaca_clone_cn$chr)))
  chr_events <- rbindlist(lapply(chromosomes, function(chromosome) {
    int_events <- rbindlist(lapply(all_clones, function(child) {
      parent <- edge_df[clone == child, parent_clone]
      child_cn_A <- alpaca_clone_cn[chr == chromosome & clone == child, pred_CN_A]
      child_cn_B <- alpaca_clone_cn[chr == chromosome & clone == child, pred_CN_B]
      parent_cn_A <- alpaca_clone_cn[chr == chromosome & clone == parent, pred_CN_A]
      parent_cn_B <- alpaca_clone_cn[chr == chromosome & clone == parent, pred_CN_B]
      interval_events_A <- get_interval_events(parent_cn = parent_cn_A, child_cn = child_cn_A)
      interval_events_B <- get_interval_events(parent_cn = parent_cn_B, child_cn = child_cn_B)
      interval_events_A$allele <- "A"
      interval_events_B$allele <- "B"
      interval_events_df <- rbind(interval_events_A, interval_events_B)
      # sum interval events across alleles:
      cols_to_sum <- colnames(interval_events_df)[colnames(interval_events_df) != "allele"]
      total_events <- paste0("total_", cols_to_sum)
      interval_events_df <- as.data.table(interval_events_df)
      interval_events_df[, (total_events) := lapply(.SD, sum), .SDcols = cols_to_sum]
      interval_events_df[, clone := child]
      interval_events_df[, chr := chromosome]
      return(interval_events_df)
    }))
    return(int_events)
  }))
  # Now sum across chromosomes
  cols_to_keep <- c("allele", "clone", "chr")
  cols_to_sum <- colnames(chr_events)[!colnames(chr_events) %in% cols_to_keep]
  fullgenome_events <- copy(chr_events)
  for (c in cols_to_sum) {
    fullgenome_events[, (c) := sum(get(c)), by = .(clone, allele)]
  }
  final_cols <- c(cols_to_sum, "clone", "allele")
  fullgenome_events <- unique(fullgenome_events[, ..final_cols])
  return(fullgenome_events)
}


### Function to compute distance events on all edges of a tree from:
# 1. mutational distance
# 2. cn interval events
# (note: doing this in long form for both alleles)
get_clone_genomic_distances <- function(alpaca_tum, tree_rlist_object, tree_paths) {
  ### 1. Mutational distance
  # SNV tree data:
  tree_structure <- convert_tree_paths_to_edgemat(tree_paths)
  trunk <- unique(tree_structure[, 1][!tree_structure[, 1] %in% tree_structure[, 2]])
  tree_structure <- rbind(tree_structure, c("diploid", trunk))
  tree_df <- data.table(tree_structure)
  colnames(tree_df) <- c("parent_clone", "clone")
  # only take clusters in the tree:
  edgelength <- tree_rlist_object$graph_pyclone$edgelength
  edgelength <- edgelength[names(edgelength) %in% unique(c(as.matrix(tree_structure)))]
  branchlength_df <- as.data.table(edgelength)
  colnames(branchlength_df) <- c("clone", "mut_branchlength")
  # merge edgelength with tree structure:
  tree_df <- merge.data.table(tree_df, branchlength_df, by = "clone", all.x = T)
  tree_df[, tmb := sum(mut_branchlength)]

  ### 2. Interval events
  # ALPACA copy number data:
  # Compute interval events for every tree edge
  int_events_df <- compute_interval_events_alpaca(edge_df = tree_df, alpaca_tum_out = alpaca_tum)  
  
  ### 4. Merge mutational distance, cn distance and interval events
  tree_df_final <- tree_df[, .(clone, parent_clone, mut_branchlength, tmb)]
  tree_df_final <- merge.data.table(tree_df_final, int_events_df, by = "clone", all = TRUE)
  tree_df_final[, n_segs := length(unique(alpaca_tum[, segment]))]
  return(tree_df_final)
}


### Function to get all ancestors of a cluster from a tree graph
get_cluster_ancestors <- function(tree_graph, clone, include_trunk = FALSE, include_cluster = FALSE) {
  allpaths <- convert_edgemat_to_tree_paths(tree_graph)
  paths_clone_ispresent_ids <- grep(paste0("\\b", clone, "\\b"), allpaths)
  ancestors <- list()
  for (p in paths_clone_ispresent_ids) {
    path <- allpaths[[p]]
    path_pos <- which(path == clone)
    if (path_pos == 1) ancests <- NULL
    else ancests <- path[1:(path_pos - 1)]
    ancestors <- append(ancestors, list(ancests))
  }
  ancestors <- unique(unlist(ancestors))
  if (include_trunk == F) {
    trunk <- allpaths[[1]][1]
    ancestors <- ancestors[ancestors != trunk]
  }
  if (include_cluster == T) {
    ancestors <- c(ancestors, clone)
  }
  return(ancestors)
}


### Function to get all descendants of a cluster from a tree graph
get_cluster_descendants <- function(tree_graph, clone, include_clone = FALSE) {
  allpaths = convert_edgemat_to_tree_paths(tree_graph)
  paths_clone_ispresent_ids = grep(paste0("\\b",clone,"\\b"),allpaths)

  descendants = list()
  for (p in paths_clone_ispresent_ids) {
    path = allpaths[[p]]
    path_pos = which(path==clone)
    if (path_pos==length(path)) descends = NULL
    else descends = path[(path_pos+1):length(path)]
    descendants = append(descendants, list(descends))
  }
  final_descendants <- unique(unlist(descendants))
  if (include_clone == TRUE) {
    final_descendants <- c(clone, final_descendants)
  }
  return(final_descendants)
}


cp_to_ccf <- function(tree_structure, cp_matrix) {
  clones_in_tree <- unique(c(tree_structure))
  cp_matrix <- cp_matrix[rownames(cp_matrix) %in% clones_in_tree, , drop = FALSE]
  ccf_table <- cp_matrix
  ccf_table[] <- 0

  regions <- colnames(cp_matrix)
  for (r in regions) {
		for (c in clones_in_tree) {
			descendents <- get_cluster_descendants(tree_structure, c)
			descendents <- c(c, descendents)
			ccf_table[rownames(ccf_table) == c, r] <- sum(cp_matrix[rownames(cp_matrix) %in% descendents, r])
		}
	}
  ccf_table <- ccf_table * 100
  return(ccf_table)
}


### Function to sum all the edge events (mutations, cn distance, cn interval events)
# along all tree paths from the MRCA to each clone.
# This function takes in the output of generate_distance_data_alpaca()
get_cumulative_events <- function(dist_df) {
  dist_out <- copy(dist_df)
  tree_structure <- as.matrix(unique(dist_df[, .(parent_clone, clone)]))
  clusters <- unique(c(as.matrix(tree_structure)))
  int_event_cols <- c("interval_gains", "interval_losses", "interval_gains_bin", "interval_losses_bin", "interval_events", "interval_events_bin")
  edge_event_cols <- c("mut_branchlength", int_event_cols, paste0("total_", int_event_cols))
  cum_cols <- paste0("cum_", edge_event_cols)
  out_cols <- c("clone", "allele", cum_cols)
  for (c in clusters) {
    # get ancestors of cluster
    ancestors <- get_cluster_ancestors(tree_graph = tree_structure, clone = c, include_trunk = TRUE, include_cluster = TRUE)
    dist_out[clone == c, cluster_ancestors := paste(ancestors, collapse = ";")]
    ancestors_df <- dist_out[clone %in% ancestors]
    ancestors_df[, (cum_cols) := lapply(.SD, sum), .SDcols = edge_event_cols, by = allele]
    dist_out[clone == c, (cum_cols) := ancestors_df[clone == c, ..cum_cols]]
  }
  return(dist_out)
}


### Function to generate dataframe of cn change on each edge on each genomic segment
get_edge_cn_df <- function(alpaca_tum, tree_paths) {
  alpaca_tum[, cluster := gsub("clone", "", clone)]
  all_tum_clones <- unique(alpaca_tum$cluster)
  all_tum_clones <- all_tum_clones[all_tum_clones != "diploid"]
  tree_structure <- convert_tree_paths_to_edgemat(tree_paths)
  trunk <- unique(tree_structure[, 1][!tree_structure[, 1] %in% tree_structure[, 2]])
  tree_structure <- rbind(tree_structure, c("diploid", trunk))
  tree_df <- data.table(tree_structure)
  colnames(tree_df) <- c("parent_clone", "clone")

  # get edge cn table
  edge_cn <- rbindlist(lapply(all_tum_clones, function(child) {
    parent <- tree_df[clone == child, parent_clone]
    alpaca_tum_child <- unique(alpaca_tum[cluster == child, .(segment, cluster, pred_CN_A, pred_CN_B)])
    alpaca_tum_parent <- unique(alpaca_tum[cluster == parent, .(segment, cluster, pred_CN_A, pred_CN_B)])
    setnames(alpaca_tum_child, c("cluster", "pred_CN_A", "pred_CN_B"), c("clone", paste0("child_", c("pred_CN_A", "pred_CN_B"))))
    setnames(alpaca_tum_parent, c("cluster", "pred_CN_A", "pred_CN_B"), c("parent_clone", paste0("parent_", c("pred_CN_A", "pred_CN_B"))))
    alpaca_edge <- merge.data.table(alpaca_tum_child, alpaca_tum_parent, "segment", all = TRUE)
    # add cluster ancestors and descendants
    anc <- get_cluster_ancestors(tree_structure, clone = child, include_trunk = TRUE, include_cluster = FALSE)
    dec <- get_cluster_descendants(tree_structure, clone = child)
    alpaca_edge[, cluster_ancestors := paste(anc, collapse = ";")]
    alpaca_edge[, cluster_descendants := paste(dec, collapse = ";")]
    return(alpaca_edge)
  }))

  edge_cn[, cn_change_A := child_pred_CN_A - parent_pred_CN_A]
  edge_cn[, cn_change_B := child_pred_CN_B - parent_pred_CN_B]
  edge_cn[, cn_change_meanallele := (cn_change_A + cn_change_B)/2]

  edge_cn[, chr := as.numeric(sapply(segment, function(s) strsplit(s, "_")[[1]][1]))]
  edge_cn[, start := as.numeric(sapply(segment, function(s) strsplit(s, "_")[[1]][2]))]
  edge_cn[, end := as.numeric(sapply(segment, function(s) strsplit(s, "_")[[1]][3]))]
  edge_cn[, seg_len := end - start]
  edge_cn[, genome_len := sum(seg_len), by = "clone"]
  setorder(edge_cn, chr, start)

  return(edge_cn)
}

### Function to compute (allele-specific) copy number event classifications for each edge, each genomic segment
classify_edge_cn_events_seglevel <- function(edge_cn_cohort, clone_metrics_df) {
  edge_change_df <- merge.data.table(edge_cn_full, clone_metrics_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)
  # Get the CN ratio from child to parent:
  edge_change_df[, child_parent_ratio_A := child_pred_CN_A / parent_pred_CN_A]
  edge_change_df[, child_parent_ratio_B := child_pred_CN_B / parent_pred_CN_B]
  # fix 0/0 = 1
  edge_change_df[child_pred_CN_A == 0 & parent_pred_CN_A == 0, child_parent_ratio_A := 1]
  edge_change_df[child_pred_CN_B == 0 & parent_pred_CN_B == 0, child_parent_ratio_B := 1]

  # Classify events: amplifications, gains, LOH
  edge_change_df[, amp_A := child_parent_ratio_A > 2 & child_pred_CN_A >= 4 & child_pred_CN_A > clone_ploidy]
  edge_change_df[, amp_B := child_parent_ratio_B > 2 & child_pred_CN_B >= 4 & child_pred_CN_B > clone_ploidy]
  edge_change_df[, gain_A := cn_change_A > 0]
  edge_change_df[, gain_B := cn_change_B > 0]
  edge_change_df[, loss_A := cn_change_A < 0]
  edge_change_df[, loss_B := cn_change_B < 0]
  edge_change_df[, loh_A := parent_pred_CN_A > 0 & child_pred_CN_A == 0]
  edge_change_df[, loh_B := parent_pred_CN_B > 0 & child_pred_CN_B == 0]
  edge_change_df[, loh_otherallele1_A := parent_pred_CN_A > 0 & child_pred_CN_A == 0 & child_pred_CN_B == 1]
  edge_change_df[, loh_otherallele1_B := parent_pred_CN_B > 0 & child_pred_CN_B == 0 & child_pred_CN_A == 1]
  # Either allele event classification
  edge_change_df[, any_amp := amp_A | amp_B]
  edge_change_df[, any_gain := gain_A | gain_B]
  edge_change_df[, any_loss := loss_A | loss_B]
  edge_change_df[, any_loh := loh_A | loh_B]
  edge_change_df[, any_loh_otherallele1 := loh_otherallele1_A | loh_otherallele1_B]
  # Homozygous dels
  edge_change_df[, homodel := (child_pred_CN_A == 0 & child_pred_CN_B == 0) & (parent_pred_CN_A > 0 | parent_pred_CN_B > 0)]
  return(edge_change_df)
}

### Tom Watkins' master GenomicRanges function
getOverlapGRanges <- function(query, subject, output="gr", ...) {
  require(GenomicRanges); require(IRanges)
  ## Input check
  if(!((class(query) == "GRanges" & class(subject) == "GRanges") | (class(query) == "IRanges" & class(subject) == "IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }

  ## Find overlapping ranges
  if(class(query) == "GRanges") {
    seqlengths(query) <- rep(NA, length(seqlengths(query)))
    seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
  }
  overlap.index.mat <- as.matrix(findOverlaps(query, subject, ...))
  query <- query[overlap.index.mat[, 1]]
  subject <- subject[overlap.index.mat[, 2]]
  overlap.mat <- cbind(query_start = start(query), query_end = end(query), subject_start = start(subject), subject_end = end(subject))

  ## Pre-queries for overlaps
  startup <- overlap.mat[,"subject_start"] < overlap.mat[,"query_start"]
  enddown <- overlap.mat[,"subject_end"] > overlap.mat[,"query_end"]
  startin <- overlap.mat[,"subject_start"] >= overlap.mat[,"query_start"] & overlap.mat[,"subject_start"] <= overlap.mat[, "query_end"]
  endin <- overlap.mat[,"subject_end"] >= overlap.mat[,"query_start"] & overlap.mat[,"subject_end"] <= overlap.mat[, "query_end"]
  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin
  contained <- startup & enddown

  ## Overlap types in one vector
  overlap_type <- rep("", length(overlap.mat[,"query_start"]))
  overlap_type[olup] <- "upstream_overlap"
  overlap_type[oldown] <- "downstream_overlap"
  overlap_type[inside] <- "subject_inside_query"
  overlap_type[contained] <- "query_contained_in_subject"

  ## Overlap positions
  overlap_start <- rep(0, length(overlap.mat[, "query_start"]))
  overlap_end <- rep(0, length(overlap.mat[, "query_start"]))
  overlap_start[olup] <- overlap.mat[, "query_start"][olup]
  overlap_end[olup] <- overlap.mat[, "subject_end"][olup]
  overlap_start[oldown] <- overlap.mat[, "subject_start"][oldown]
  overlap_end[oldown] <- overlap.mat[, "query_end"][oldown]
  overlap_start[inside] <- overlap.mat[, "subject_start"][inside]
  overlap_end[inside] <- overlap.mat[, "subject_end"][inside]
  overlap_start[contained] <- overlap.mat[, "query_start"][contained]
  overlap_end[contained] <- overlap.mat[, "query_end"][contained]
  ## Absolute and relative length of overlaps
  overlap_length <- (overlap_end - overlap_start) + 1
  overlap_prop_query <- overlap_length / width(query) * 100
  overlap_prop_subject <- overlap_length / width(subject) * 100

  ## Output type
  oldf <- data.frame(query_hit_idx=overlap.index.mat[,1], subject_hit_idx=overlap.index.mat[,2], overlap.mat, overlap_start, overlap_end, overlap_length, overlap_prop_query, overlap_prop_subject, overlap_type)
  if(class(query) == "GRanges") {
    oldf <- cbind(chrom=as.character(seqnames(query)), oldf)
  }
  if(output=="df") {
    return(cbind(oldf, as.data.frame(mcols(subject))))
  }
  if(output=="gr") {
    if(class(query)=="GRanges") {
      elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf, mcols(subject))
    }
    if(class(query)=="IRanges") {#Todo add the metadata in th eiranges
      query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"query_start"], end=oldf[,"query_end"]), strand = Rle(strand(rep("+", length(query)))), oldf)
    }
    return(query)
  }
}

### Function to generate dataframe of cn changes and cn events for each chromosome arm on each tree edge, based on input thresholds
get_chrarm_cn_events_tumour <- function(edge_change_tum, chr_arm_df, arm_event_threshold, armlohparent_threshold) {
  edge_change_tum[, pred_CN_A := child_pred_CN_A]
  edge_change_tum[, pred_CN_B := child_pred_CN_B]
  
  ### add average CN per chromosome arm
  # create genomic bin granges objects
  bin_ir <- IRanges(start = chr_arm_df$arm_start, end = chr_arm_df$arm_end)
  bin_gr <- GRanges(seqnames = chr_arm_df$chromosome, ranges = bin_ir)
  mcols(bin_gr) <- chr_arm_df

  # create segments granges objects
  seg_ir <- IRanges(start = edge_change_tum$start, end = edge_change_tum$end)
  seg_gr <- GRanges(seqnames = edge_change_tum$chr, ranges = seg_ir)
  mcols(seg_gr) <- edge_change_tum

  # get overlaps between genomic bins and segments
  overlap <- getOverlapGRanges(query = seg_gr, subject = bin_gr)
  overlap_dt <- as.data.table(overlap)
  setorder(overlap_dt, chr, start)
  cols <- c(colnames(edge_change_tum), "chr_arm", "arm_length")
  alpaca_clone_arm_cn <- unique(overlap_dt[, ..cols])
  alpaca_clone_arm_cn[, arm_fraction := seg_len / arm_length]

  # if segment spans two chrarms, split into two:
  large_segs <- unique(alpaca_clone_arm_cn[arm_fraction > 1, .(segment, chr_arm)])
  if (nrow(large_segs) > 0) {
    large_seg_cn <- rbindlist(lapply(seq(nrow(large_segs)), function(i) {
      s <- large_segs[i, segment]
      chrarm <- large_segs[i, chr_arm]
      seg_start <- alpaca_clone_arm_cn[segment == s & chr_arm == chrarm]
      original_seg_end <- seg_start[, unique(end)]
      arm_l <- seg_start[, unique(arm_length)]
      seg_start[, end := arm_length]
      seg_start[, segment := paste(chr, start, end, sep = "_")]
      seg_start[, arm_fraction := (end - start)/arm_length]
      # add additional segment that is overlapping the q arm
      if (chrarm == "p") {
        seg_end <- alpaca_clone_arm_cn[segment == s & chr_arm == "q"]
        seg_end[, start := arm_l]
        seg_end[, end := original_seg_end]
        seg_end[, segment := paste(chr, start, end, sep = "_")]
        seg_end[, arm_fraction := (end - start)/arm_length]
        seg_to_add <- rbind(seg_start, seg_end)
      } else {
        seg_to_add <- seg_start
      }
      return(seg_to_add) 
    }))
    small_seg_cn <- unique(alpaca_clone_arm_cn[arm_fraction <= 1])
    alpaca_clone_arm_cn <- rbind(small_seg_cn, large_seg_cn)
    setorder(alpaca_clone_arm_cn, chr, start)
  }
  alpaca_clone_arm_cn[, arm_size_covered := sum(seg_len), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, total_arm_fraction_covered := arm_size_covered / arm_length]
  alpaca_clone_arm_cn[, arm_fraction_covered := seg_len / arm_size_covered]
  # Add average CN change per chromosome arm 
  alpaca_clone_arm_cn[, pred_CN_A_w_mean_chr_arm := sum(pred_CN_A * arm_fraction_covered), by = .(clone, chr, chr_arm)]
  alpaca_clone_arm_cn[, pred_CN_B_w_mean_chr_arm := sum(pred_CN_B * arm_fraction_covered), by = .(clone, chr, chr_arm)]
  alpaca_clone_arm_cn[, total_pred_CN_w_mean_chr_arm := pred_CN_A_w_mean_chr_arm + pred_CN_B_w_mean_chr_arm]
  ### Call arm events like in Tx100: if 90% of the arm is lost/gained/amp
  # AMP
  alpaca_clone_arm_cn[, size_arm_amp_A := sum(seg_len * amp_A), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, size_arm_amp_B := sum(seg_len * amp_B), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, frac_arm_amp_A := size_arm_amp_A / arm_size_covered]
  alpaca_clone_arm_cn[, frac_arm_amp_B := size_arm_amp_B / arm_size_covered]
  alpaca_clone_arm_cn[, arm_amp_acquired_edge_A := frac_arm_amp_A > arm_event_threshold]
  alpaca_clone_arm_cn[, arm_amp_acquired_edge_B := frac_arm_amp_A > arm_event_threshold]
  # GAIN
  alpaca_clone_arm_cn[, size_arm_gain_A := sum(seg_len * gain_A), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, size_arm_gain_B := sum(seg_len * gain_B), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, frac_arm_gain_A := size_arm_gain_A / arm_size_covered]
  alpaca_clone_arm_cn[, frac_arm_gain_B := size_arm_gain_B / arm_size_covered]
  alpaca_clone_arm_cn[, arm_gain_acquired_edge_A := frac_arm_gain_A > arm_event_threshold]
  alpaca_clone_arm_cn[, arm_gain_acquired_edge_B := frac_arm_gain_B > arm_event_threshold]
  # LOSS
  alpaca_clone_arm_cn[, size_arm_loss_A := sum(seg_len * loss_A), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, size_arm_loss_B := sum(seg_len * loss_B), by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, frac_arm_loss_A := size_arm_loss_A / arm_size_covered]
  alpaca_clone_arm_cn[, frac_arm_loss_B := size_arm_loss_B / arm_size_covered]
  alpaca_clone_arm_cn[, arm_loss_acquired_edge_A := frac_arm_loss_A > arm_event_threshold]
  alpaca_clone_arm_cn[, arm_loss_acquired_edge_B := frac_arm_loss_B > arm_event_threshold]
  # LOH called differently -> if >= 90% of the arm is at 0, but <=20% of the parent arm was at 0
  alpaca_clone_arm_cn[, child_frac_arm_0_A := sum(seg_len * (child_pred_CN_A == 0)) / arm_size_covered, by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, child_frac_arm_0_B := sum(seg_len * (child_pred_CN_B == 0)) / arm_size_covered, by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, child_frac_arm_0 := sum(seg_len * ((child_pred_CN_A == 0) | (child_pred_CN_B == 0))) / arm_size_covered, by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, parent_frac_arm_0_A := sum(seg_len * (parent_pred_CN_A == 0)) / arm_size_covered, by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, parent_frac_arm_0_B := sum(seg_len * (parent_pred_CN_B == 0)) / arm_size_covered, by = .(chr, chr_arm, clone)]
  alpaca_clone_arm_cn[, arm_loh_acquired_edge_A := (parent_frac_arm_0_A <= armlohparent_threshold) & (child_frac_arm_0_A >= arm_event_threshold)]
  alpaca_clone_arm_cn[, arm_loh_acquired_edge_B := (parent_frac_arm_0_B <= armlohparent_threshold) & (child_frac_arm_0_B >= arm_event_threshold)]
  ### Summary - events called in any allele
  alpaca_clone_arm_cn[, any_arm_loh_acquired_edge := arm_loh_acquired_edge_A | arm_loh_acquired_edge_B]
  alpaca_clone_arm_cn[, any_arm_amp_acquired_edge := arm_amp_acquired_edge_A | arm_amp_acquired_edge_B]
  alpaca_clone_arm_cn[, any_arm_gain_acquired_edge := arm_gain_acquired_edge_A | arm_gain_acquired_edge_B]
  alpaca_clone_arm_cn[, any_arm_loss_acquired_edge := arm_loss_acquired_edge_A | arm_loss_acquired_edge_B]
  return(alpaca_clone_arm_cn)
}


### Function to overlap alpaca cohort tumour segments with all genes locations:
overlap_genes_segments_cohort <- function(genes_df, segs_df) {
  require(data.table)
  require(GenomicRanges)
  
  # Get genes overlapped with edge events
  # create gene granges objects
  genes_gr <- GRanges(seqnames = genes_df$chr, ranges = IRanges(start = genes_df$start, end = genes_df$end))
  mcols(genes_gr) <- genes_df
  # create segments granges objects
  segs_gr <- GRanges(seqnames = segs_df$chr, ranges = IRanges(start = segs_df$start, end = segs_df$end))
  mcols_alpaca <- copy(segs_df)
  mcols_alpaca$chr <- NULL
  mcols_alpaca$start <- NULL
  mcols_alpaca$end <- NULL
  mcols(segs_gr) <- mcols_alpaca

  # Get overlaps between genomic genes and segments
  print("Setting query as genes df, subject as segments df in overlap")
  overlap <- getOverlapGRanges(query = genes_gr, subject = segs_gr)
  overlap_dt <- as.data.table(overlap)

  # Choose segment that overlaps each gene the most
  print("Assigning genes to segments with biggest overlap")
  overlap_dt[, max_overlap_seg := max(overlap_prop_query), by = .(tumour_id, gene_name)]
  overlap_dt[, best_overlap_seg := overlap_prop_query == max_overlap_seg]
  overlaps_genesubset <- overlap_dt[best_overlap_seg == TRUE]
  overlaps_genesubset[, n_segs_overlapping_gene := .N, by = .(tumour_id, gene_name)]
  # If >1 segments equally overlap the gene, pick the biggest segment
  overlaps_genesubset[, subject_width := subject_end - subject_start]
  overlaps_genesubset[, max_seg_size_per_gene := max(subject_width), by = .(tumour_id, gene_name)]
  overlaps_genesubset[, biggest_overlapping_seg := subject_width == max_seg_size_per_gene]
  overlaps_genesubset <- overlaps_genesubset[biggest_overlapping_seg == TRUE]
  return(overlaps_genesubset)
}

### Function to process scna distance vs #snvs per tumour edge
process_edge_snv_scna_dist <- function(edge_distance_df) {
  
  ### edge distances
  edge_distance_df[, clone := as.character(clone)]
  cols <- c("tumour_id", "clone", "parent_clone", "mut_branchlength", 
            "tmb", "n_segs", "cluster_ancestors", 
            "total_interval_gains", "total_interval_losses", "total_interval_events",
            "total_interval_gains_bin", "total_interval_losses_bin", "total_interval_events_bin"
            # ,"l1_mean_cn_branchlength", "l2_mean_cn_branchlength", "l1_mean_cn_distance_MRCA", "l2_mean_cn_distance_MRCA"
            )
  distance_total_df <- unique(edge_distance_df[, ..cols])
  distance_total_df[, is_mrca := ifelse(parent_clone == "diploid", "mrca", "subclone")]

  # add total no. copy number events per tumour
  distance_total_df[, sum_cn_amplitude_events_tumour := sum(total_interval_events), by = tumour_id]
  distance_total_df[, sum_cn_events_tumour := sum(total_interval_events_bin), by = tumour_id]
  distance_total_df[, sum_subclonal_cn_events_tumour := sum(total_interval_events_bin), by = .(tumour_id, is_mrca)]
  distance_total_df[, subclonal_tmb := sum(mut_branchlength), by = .(tumour_id, is_mrca)]

  distance_total_df[, sum_amplitude_gains_tumour := sum(total_interval_gains), by = tumour_id]
  distance_total_df[, sum_amplitude_losses_tumour := sum(total_interval_losses), by = tumour_id]
  distance_total_df[, sum_gains_tumour := sum(total_interval_gains_bin), by = tumour_id]
  distance_total_df[, sum_subclonal_gains_tumour := sum(total_interval_gains_bin), by = .(tumour_id, is_mrca)]
  distance_total_df[, sum_losses_tumour := sum(total_interval_losses_bin), by = tumour_id]
  distance_total_df[, sum_subclonal_losses_tumour := sum(total_interval_losses_bin), by = .(tumour_id, is_mrca)]

  # normalise branch lengths by subclonal tmb/cn events
  distance_total_df[, branchlength_by_tmb := mut_branchlength / tmb]
  distance_total_df[, total_cn_events_norm := total_interval_events_bin / sum_cn_events_tumour]
  distance_total_df[, total_gains_norm := total_interval_gains_bin / sum_cn_events_tumour]
  distance_total_df[, total_losses_norm := total_interval_losses_bin / sum_cn_events_tumour]
  
  distance_total_df[, total_cn_events_norm_subclonal := total_interval_events_bin / sum_subclonal_cn_events_tumour]
  distance_total_df[, total_gains_norm_subclonal := total_interval_gains_bin / sum_subclonal_cn_events_tumour]
  distance_total_df[, total_losses_norm_subclonal := total_interval_losses_bin / sum_subclonal_cn_events_tumour]

  # distance_total_df[, sum_subclonal_l1_cn_branchlength_tumour := sum(l1_mean_cn_branchlength), by = .(tumour_id, is_mrca)]
  # distance_total_df[, sum_subclonal_l2_cn_branchlength_tumour := sum(l2_mean_cn_branchlength), by = .(tumour_id, is_mrca)]
  # distance_total_df[, total_l1_cn_branchlength_norm_subclonal := l1_mean_cn_branchlength / sum_subclonal_l1_cn_branchlength_tumour]
  # distance_total_df[, total_l2_cn_branchlength_norm_subclonal := l2_mean_cn_branchlength / sum_subclonal_l2_cn_branchlength_tumour]

  distance_total_df[, mut_branchlength_norm := mut_branchlength / tmb]
  distance_total_df[, mut_branchlength_norm_subclonal := mut_branchlength / subclonal_tmb]

  # n clones
  distance_total_df[, n_clones := length(unique(clone)), by = tumour_id]

  # Compute SCNA/SNV events slope
  distance_total_df[, slope_raw := total_interval_events_bin / mut_branchlength]
  distance_total_df[, slope_norm_all := total_cn_events_norm / mut_branchlength_norm]
  distance_total_df[, slope_norm := total_cn_events_norm_subclonal / mut_branchlength_norm_subclonal]
  distance_total_df[, slope_norm_meantum := mean(slope_norm), by = tumour_id]
  distance_total_df[, patient_id := gsub("_Cluster.*", "", tumour_id)]
  distance_total_df[, clone := as.character(clone)]

  return(distance_total_df)
}


### Function to generate a dataframe for each tree path per tumour
generate_tree_paths_df <- function(alpaca_edge_df, verbose = TRUE) {
  all_tumours <- sort(unique(alpaca_edge_df$tumour_id))
  
  full_tree_paths_df <- rbindlist(lapply(all_tumours, function(tum_id) {
    if (verbose == TRUE) print(tum_id)
    tree_structure <- as.matrix(alpaca_edge_df[tumour_id == tum_id & parent_clone != "diploid", .(parent_clone, clone)])
    tree_paths <- convert_edgemat_to_tree_paths(tree_structure)
    tree_paths_df <- rbindlist(lapply(tree_paths, function(p){
        path_dt <- data.table(path_id = paste(p, collapse = '_'),
                                clone = as.character(p),
                                child_clones_path = p[c(2:length(p), length(p))])
        return(path_dt)
    }))
    tree_paths_df[, tumour_id := tum_id]
  }))
  full_tree_paths_df[, clone := as.character(clone)]
  return(full_tree_paths_df)
}

### Function to add on columns calculating cumulative branch lengths
get_cumulative_snv_scna_events <- function(distance_total_df) {
  working_df <- copy(distance_total_df)
  edge_event_cols <- c("mut_branchlength", 
                        # "l1_mean_cn_branchlength",
                        # "l2_mean_cn_branchlength",
                        "total_no_events", "total_no_gains", "total_no_losses")
  cum_cols <- paste0("cum_", edge_event_cols)
  working_df[, c("cum_mut_branchlength",
                  # "cum_l1_cn_branchlength",
                  # "cum_l2_cn_branchlength",
                  "cum_total_no_events",
                  "cum_total_no_gains",
                  "cum_total_no_losses") :=
  rbindlist(lapply(seq(nrow(working_df)), function(i) {
    # print(i)
    tum_id <- working_df[i, tumour_id]
    clone_id <- working_df[i, clone]

    tree_structure <- as.matrix(unique(working_df[tumour_id == tum_id & parent_clone != "diploid", .(parent_clone, clone)]))
    ancestors_vec <- get_cluster_ancestors(tree_structure, clone = clone_id, include_trunk = TRUE, include_cluster = TRUE)

    pcdist <- working_df[tumour_id == tum_id]

    ancestors_df <- pcdist[clone %in% ancestors_vec]
    ancestors_df[, (cum_cols) := lapply(.SD, sum), .SDcols = edge_event_cols]

    clone_out <- ancestors_df[clone == clone_id, ..cum_cols]
    # print(clone_out)
    return(clone_out)
  }))]
  return(working_df)
}

### Function to get the mean slope per tree paths ending in metastatic clones vs ending in primary clones

get_met_vs_prim_path_slope_diff <- function(distance_total_df, tree_paths_df) {
  dist_paths <- merge.data.table(distance_total_df, tree_paths_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)
  # Assign each path to a class based on the class of its terminal nodes
  path_classes <- dist_paths[clone == child_clones_path, .(tumour_id, path_id, clone_class)]
  path_classes[, path_class := as.character(clone_class)]
  # Consider only two path classes: Metastatic and Primary
  path_classes[path_class == "Seeding", path_class := "Metastatic"]
  dist_path_class <- merge.data.table(dist_paths, path_classes[, .(tumour_id, path_id, path_class)], by = c("tumour_id", "path_id"), all.x = TRUE, all.y = FALSE)

  # get mean slope per path class
  slope_class_path <- unique(dist_path_class[parent_clone != "diploid", .(tumour_id, clone, path_class, slope_raw)])
  slope_class_path[, mean_slope_clone_class_path := mean(slope_raw), by = .(tumour_id, path_class)]
  slope_class_path <- unique(slope_class_path[, .(tumour_id, path_class, mean_slope_clone_class_path)])
  slope_class_wide <- dcast(slope_class_path, formula = tumour_id ~ path_class, value.var = "mean_slope_clone_class_path")
  colnames(slope_class_wide)[2:3] <- paste0(colnames(slope_class_wide)[2:3], "_path_mean_slope")
  # Get difference in path class between Mets and Primary path per tumour
  slope_class_wide[, path_class_slope_diff := Metastatic_path_mean_slope - Primary_path_mean_slope]

  # Merge with class path average slope
  dist_path_class <- merge.data.table(dist_path_class, slope_class_path, by = c("tumour_id", "path_class"), all.x = TRUE)
  # Merge with class path average slope - wide version for comparison
  dist_path_class <- merge.data.table(dist_path_class, slope_class_wide, by = "tumour_id", all.x = TRUE)
  dist_path_class[, met_slope_higher := path_class_slope_diff > 0]
  return(dist_path_class)
}



####################################################################
### Functions for generating a tumour tree paths edge count plot ###
####################################################################

generate_mut_cn_rate_tree <- function(tum_distance_df, tum_tree_paths_df, colour_map, colour_pal, patient_id) {
  require(data.table)
  require(ggplot2)
  # merge with tree paths
  pcdist_paths <- merge.data.table(tum_distance_df, tum_tree_paths_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)
  pcdist_paths <- merge.data.table(pcdist_paths, colour_map, by = "clone", all.x = TRUE, all.y = FALSE)
  
  # Colour by the child clone colour
  child_colour_map <- colour_map[, .(child_clones_path = clone, child_colour_label = colour_label)]
  pcdist_paths <- merge.data.table(pcdist_paths, child_colour_map, by = "child_clones_path", all.x = TRUE, all.y = FALSE)
  g <- ggplot(pcdist_paths, aes(cum_mut_branchlength, cum_total_no_events))+
    geom_line(aes(group = path_id, colour = as.character(child_colour_label)), size = 1) +
    geom_point(size = 2.5, aes(colour = as.character(colour_label))) +
    # geom_text(size = 8, aes(label = clone)) +
    theme_void() +
    labs(title = patient_id) +
    scale_colour_manual(values = colour_pal, name = 'Clone class') +
    theme(legend.position = "none", plot.title = element_text(size = 19))
  return(g)
}

generate_mut_cn_rate_tree_loss_gain_plot <- function(tum_distance_df, tum_tree_paths_df, colour_map, colour_pal, xmax=NULL, ymax=NULL){
  require(data.table)
  require(ggplot2)
  # merge with tree paths
  pcdist_paths <- merge.data.table(tum_distance_df,
                                   tum_tree_paths_df,
                                   by = c("tumour_id", "clone"),
                                   all.x = TRUE,
                                   all.y = FALSE
  )
  pcdist_paths <- merge.data.table(pcdist_paths,
                                    colour_map,
                                    by = "clone",
                                    all.x = TRUE,
                                    all.y = FALSE)
  
  # Colour by the child clone colour
  child_colour_map <- colour_map[, .(child_clones_path = clone,
                                    child_colour_label = colour_label)]
  pcdist_paths <- merge.data.table(pcdist_paths,
                                    child_colour_map,
                                    by = "child_clones_path",
                                    all.x = TRUE,
                                    all.y = FALSE)
  g <- ggplot(pcdist_paths,
          aes(cum_total_no_gains, cum_total_no_losses))+
    geom_line(aes(group = path_id, colour = as.character(child_colour_label)), size = 1) +
    geom_point(size = 2.5, aes(colour = as.character(colour_label))) +
    # geom_text(size = 8, aes(label = clone)) +
    theme_void() +
    scale_colour_manual(values = colour_pal) +
    theme(legend.position = "none")
  if (!is.null(xmax) & !is.null(ymax)) {
    g <- g +
     xlim(c(min(pcdist_paths$cum_total_no_gains, na.rm = TRUE), xmax)) +
     ylim(c(min(pcdist_paths$cum_total_no_losses, na.rm = TRUE), xmax))
  }
  return(g)
}

generate_mut_cnchange_rate_tree <- function(tum_distance_grad_df, tum_tree_paths_df, colour_map, colour_pal, patient_id) {
  require(data.table)
  require(ggplot2)
  # merge with tree paths
  # if ()
  # pcdist_paths <- merge.data.table(tum_distance_grad_df, tum_tree_paths_df, by = c("tumour_id", "clone"), all.x = TRUE, all.y = FALSE)
  pcdist_paths <- merge.data.table(tum_distance_grad_df, colour_map, by = "clone", all.x = TRUE, all.y = FALSE)
  
  # Colour by the child clone colour
  child_colour_map <- colour_map[, .(child_clones_path = clone, child_colour_label = colour_label)]
  pcdist_paths <- merge.data.table(pcdist_paths, child_colour_map, by = "child_clones_path", all.x = TRUE, all.y = FALSE)
  g <- ggplot(pcdist_paths, aes(cum_mut_branchlength, slope_raw))+
    geom_line(aes(group = path_id, colour = as.character(child_colour_label)), size = 1) +
    geom_point(size = 2.5, aes(colour = as.character(colour_label))) +
    geom_hline(yintercept = 0) +
    # geom_text(size = 8, aes(label = clone)) +
    theme_void() +
    labs(title = patient_id) +
    scale_colour_manual(values = colour_pal, name = 'Clone class') +
    theme(legend.position = "none", plot.title = element_text(size = 19))
  return(g)
}

######################################
### Functions for seeding analysis ###
######################################

### Function to assign clone class to clones in metastatic transition 
assign_clone_classes <- function(edge_clone_info_df) {
  tmp_df <- copy(edge_clone_info_df)
  # Add clone class assignment
  print('Assigning clone class')
  tmp_df[, clone_class := ifelse(seedingClones == TRUE, 'Seeding', 
                              ifelse(parent_clone == "diploid", 'Non-seeding MRCA',
                                ifelse(sharedClones == TRUE, 'Shared',
                                  ifelse(primaryClones == TRUE, 'Primary',
                                    ifelse(metClones == TRUE, 'Metastatic', NA)))))]
  tmp_df[, clone_class := factor(clone_class, levels = c('Non-seeding MRCA', 'Shared', 'Primary', 'Metastatic', 'Seeding'))]
  tmp_df[, mrca_seeding_tum := any(parent_clone == "diploid" & seedingClones == TRUE), by = tumour_id]

  # Classify whether any ancestors of the clone are seeding
  edge_cloneInfo_df <- unique(tmp_df[, .(tumour_id, clone, parent_clone, clone_class, clonalClones, seedingClones,
                                          sharedClones, primaryClones, metClones, cluster_ancestors
  )])
  print('Computing seeding ancestors')
  edge_cloneInfo_df[, any_ancestors_seeding := sapply(seq(nrow(edge_cloneInfo_df)), function(i) {
    tum_id <- edge_cloneInfo_df[i, tumour_id]
    ancestors <- strsplit(edge_cloneInfo_df[i, cluster_ancestors], ";")[[1]]
    seeding_clones <- edge_cloneInfo_df[tumour_id == tum_id & seedingClones == TRUE, clone]
    return(any(ancestors %in% seeding_clones))
  })]

  # Get pre-post-non seeding classification
  print('Assinging pre-post seeding class')
  edge_cloneInfo_df[, seeding_class := ifelse(primaryClones == TRUE & any_ancestors_seeding == FALSE, "Non-seeding",
                                        ifelse(sharedClones == TRUE, "Pre-seeding",
                                          ifelse(metClones == TRUE, "Post-seeding, metastasis",
                                            ifelse(primaryClones == TRUE & any_ancestors_seeding == TRUE, "Post-seeding, primary", NA)
  )))]
  edge_cloneInfo_df[, nonseeding_path_present := any(grepl("Non-seeding", seeding_class)), by = tumour_id]
  
  # Merge back into full edge df:
  mergecols <- colnames(edge_cloneInfo_df)[colnames(edge_cloneInfo_df) %in% colnames(tmp_df)]
  tmp_df <- merge.data.table(tmp_df, edge_cloneInfo_df, by = mergecols, all.x = TRUE, all.y = FALSE)
  return(tmp_df)
}

### Function to compare events to mrca for each subclone:
compare_events_to_mrca <- function(tum_edge_df){
  out_df <- copy(tum_edge_df)
  # Extract MRCA:
  cn_mrca <- out_df[parent_clone == "diploid", .(segment, pred_CN_A, pred_CN_B)]
  setnames(cn_mrca, c("pred_CN_A", "pred_CN_B"), c("MRCA_pred_CN_A", "MRCA_pred_CN_B"))
  # Merge mrca cn with full cn bin table
  out_df <- merge.data.table(out_df, cn_mrca, by = "segment", all = TRUE, allow.cartesian = TRUE)
  # For seeding clones, for each segment, get difference between MRCA clone and seeding clone cn_A, cn_B:
  out_df[, pred_CN_diff_to_mrca_A := pred_CN_A - MRCA_pred_CN_A]
  out_df[, pred_CN_diff_to_mrca_B := pred_CN_B - MRCA_pred_CN_B]
  return(out_df)
}


### Function to get events from mrca --> [a clone class]
# MRCA --> seeding/nonseeding:
# Get average predicted CN difference between MRCA and all the clones of class X [i.e. seeding clones/nonseeding clones]
# e.g. alpaca_cloneclass <- alpaca_gene_seeding[clone_class_full == "seeding"]
get_mrca_cloneclass_events <- function(alpaca_cloneclass) {

  cn_mrca_cloneclass <- copy(alpaca_cloneclass)

  # Classify event from MRCA --> cloneclass clones:
  # (include a threshold to not over-count events)

  cn_mrca_cloneclass[, CN_ratio_to_mrca_A := pred_CN_A / MRCA_pred_CN_A]
  cn_mrca_cloneclass[, CN_ratio_to_mrca_B := pred_CN_B / MRCA_pred_CN_B]
  
  # FIX ZEROS
  cn_mrca_cloneclass[MRCA_pred_CN_A == 0, CN_ratio_to_mrca_A := 1]
  cn_mrca_cloneclass[MRCA_pred_CN_B == 0, CN_ratio_to_mrca_B := 1]

  cn_mrca_cloneclass[, amp_A := any(CN_ratio_to_mrca_A >= 2 & pred_CN_A >= 4 & pred_CN_A > clone_ploidy), by = segment]
  cn_mrca_cloneclass[, amp_B := any(CN_ratio_to_mrca_B >= 2 & pred_CN_B >= 4 & pred_CN_B > clone_ploidy), by = segment]
  cn_mrca_cloneclass[, is_amp := amp_A | amp_B]
  
  cn_mrca_cloneclass[, gain_A := any(pred_CN_diff_to_mrca_A > 0), by = segment]
  cn_mrca_cloneclass[, gain_B := any(pred_CN_diff_to_mrca_B > 0), by = segment]
  cn_mrca_cloneclass[, is_gain := gain_A | gain_B]
  
  cn_mrca_cloneclass[, loss_A := any(pred_CN_diff_to_mrca_A < 0), by = segment]
  cn_mrca_cloneclass[, loss_B := any(pred_CN_diff_to_mrca_B < 0), by = segment]
  cn_mrca_cloneclass[, is_loss := loss_A | loss_B]
  
  cn_mrca_cloneclass[, gain_1_A := any(pred_CN_diff_to_mrca_A > 1), by = segment]
  cn_mrca_cloneclass[, gain_1_B := any(pred_CN_diff_to_mrca_B > 1), by = segment]
  cn_mrca_cloneclass[, is_gain_1 := gain_1_A | gain_1_B]
  
  cn_mrca_cloneclass[, loss_1_A := any(pred_CN_diff_to_mrca_A < -1), by = segment]
  cn_mrca_cloneclass[, loss_1_B := any(pred_CN_diff_to_mrca_B < -1), by = segment]
  cn_mrca_cloneclass[, is_loss_1 := loss_1_A | loss_1_B]
  
  cn_mrca_cloneclass[, gain_2_A := any(pred_CN_diff_to_mrca_A > 2), by = segment]
  cn_mrca_cloneclass[, gain_2_B := any(pred_CN_diff_to_mrca_B > 2), by = segment]
  cn_mrca_cloneclass[, is_gain_2 := gain_2_A | gain_2_B]
  
  cn_mrca_cloneclass[, loss_2_A := any(pred_CN_diff_to_mrca_A < -2), by = segment]
  cn_mrca_cloneclass[, loss_2_B := any(pred_CN_diff_to_mrca_B < -2), by = segment]
  cn_mrca_cloneclass[, is_loss_2 := loss_2_A | loss_2_B]

  cn_mrca_cloneclass[, loh_A := any((pred_CN_diff_to_mrca_A < 0 & pred_CN_A == 0)), by = segment]
  cn_mrca_cloneclass[, loh_B := any((pred_CN_diff_to_mrca_B < 0 & pred_CN_B == 0)), by = segment]
  cn_mrca_cloneclass[, is_loh := loh_A | loh_B]
  
  cn_mrca_cloneclass[, loh_01_A := any((pred_CN_diff_to_mrca_A < 0 & pred_CN_A == 0) & (pred_CN_B < 1.1)), by = segment]
  cn_mrca_cloneclass[, loh_01_B := any((pred_CN_diff_to_mrca_B < 0 & pred_CN_B == 0) & (pred_CN_A < 1.1)), by = segment]
  cn_mrca_cloneclass[, is_loh_01 := loh_01_A | loh_01_B]

  # Extract summarised mrca --> cloneclass dataframe for each segment:
  cn_mrca_cloneclass_events <- unique(cn_mrca_cloneclass[,
    .(segment,
    is_amp,
    is_loh,
    is_loh_01,
    is_gain,
    is_loss,
    is_gain_1,
    is_loss_1,
    is_gain_2,
    is_loss_2
  )])

  return(cn_mrca_cloneclass_events)
}

  
compute_mrca_seeding_nonseeding_events <- function(metstum_edge_change_arm_df, assets_dir) {
  # Compare events to mrca
  alpaca_seedingclass <- compare_events_to_mrca(metstum_edge_change_arm_df)
  
  # 2) Compute mrca --> seeding events
  alpaca_seeding_clones <- alpaca_seedingclass[clone_class == 'Seeding']
  cn_mrca_seeding_events <- get_mrca_cloneclass_events(alpaca_cloneclass = alpaca_seeding_clones)
  cn_mrca_seeding_events[, event_class := "MRCA --> Seeding"]

  # 3) get mrca --> non-seeding events
  # restrict to terminal nodes
  alpaca_nonseeding_clones <- alpaca_seedingclass[seeding_class == "Non-seeding"]
  cn_mrca_nonseeding_events <- get_mrca_cloneclass_events(alpaca_cloneclass = alpaca_nonseeding_clones)
  cn_mrca_nonseeding_events[, event_class := "MRCA --> Non-seeding"]

  # merge events into one df
  seeding_class_event_df <- rbindlist(list(cn_mrca_seeding_events, cn_mrca_nonseeding_events))
  
  # add on tumour specific columns:
  tum_cols <- unique(alpaca_seedingclass[, .(tumour_id, segment, chr, start, end)], by = "segment")
  seeding_class_event_df <- merge.data.table(seeding_class_event_df, tum_cols, by = "segment", all.x = TRUE, all.y = FALSE)
  
  # save output
  return(seeding_class_event_df)
}


### Functions to get summary counts of events
compute_event_total_counts <- function(cn_events_df, events_to_test) {
  cn_events_copy <- copy(cn_events_df)
  total_no_tumours <- length(unique(cn_events_copy$tumour_id))
  # Summarise events across all tumours
  for (event in events_to_test) {
    print(event)
    n_col_name <- gsub("is_", "no_tums_", event)
    ntotal_col_name <- gsub("is_", "total_no_tums_", event)
    frac_col_name <- gsub("is_", "frac_tums_", event)
    cn_events_copy[, eval(n_col_name) := sum(get(event)), by = .(gene_name, event_class)]
    cn_events_copy[, eval(frac_col_name) := get(n_col_name) / total_no_tumours]
    cn_events_copy[, eval(ntotal_col_name) := sum(get(event)), by = .(gene_name)]
  }
  return(cn_events_copy)
}

summarise_event_counts <- function(cn_events_df, events_to_test) {
  cn_events_total_df <- compute_event_total_counts(cn_events_df, events_to_test)
  final_cols <- c("event_class", "chr", "start", "end", "gene_name", colnames(cn_events_total_df)[grepl("total_no_|no_|frac_", colnames(cn_events_total_df))])
  cn_events_summary <- unique(cn_events_total_df[, ..final_cols])
  return(cn_events_summary)
}

### Function to compute Fisher"s exact test of gene event compared to all genes
# between two cluster classes
compute_fishersexact <- function(gene_df){
  cont_table <- as.matrix(gene_df)
  rownames(cont_table) <- cont_table[, 1]
  cont_table <- cont_table[, 2:ncol(cont_table)]
  cont_table <- apply(cont_table, c(1, 2), as.numeric)
  ft <- fisher.test(cont_table)
  gene_df[, fisher_pval := ft$p.value]
  return(gene_df)
}

### Function to run fishers exact test for driver gene list against background for certain event type
run_paired_fishersexact_driv_vs_bg_event <- function(gene_class_cn_events_summary, event_col, driv_gene_list, verbose = TRUE, n_tums_event_thresh = NULL) {
  tmp <- copy(gene_class_cn_events_summary)
  print("Restricting testing to only genes with at least n_tums_event_thresh events")

  if (is.null(n_tums_event_thresh)) {  
    n_tums_event_thresh <- 1
  }
  n_tums_event_thresh <- as.numeric(n_tums_event_thresh)
  print("n_tums_event_thresh:")
  print(n_tums_event_thresh)
  
  total_gene_event_col <- paste0("total_gene_", event_col)
  tmp[, eval(total_gene_event_col) := sum(get(event_col)), by = .(gene_name)]
  tmp <- tmp[get(total_gene_event_col) >= n_tums_event_thresh]
  
  print("Enumerating total number of background events in each class")
  total_event_col <- paste0("total_", event_col)
  tmp[, eval(total_event_col) := sum(get(event_col)), by = event_class]
  
  print("Testing each gene against background")
  driv_gene_list <- driv_gene_list[driv_gene_list %in% unique(tmp$gene_name)]
  
  driv_events <- rbindlist(lapply(driv_gene_list, function(gene) {
    if (verbose == TRUE) print(gene)
    gene_event <- tmp[gene_name == gene, .(event_class, total_event_col = get(total_event_col), event_col = get(event_col))]
    # print("Computing fisher's exact test")
    gene_event_df <- compute_fishersexact(gene_df = gene_event)
    setnames(gene_event_df, c("total_event_col", "event_col"), c(total_event_col, event_col))
    gene_event_df[, gene_name := gene]
    # print(gene_event_df)
    return(gene_event_df)
  }))
  all_events <- tmp[, .(N = sum(get(event_col))), by = event_class]
  all_events[, gene_name := "background"]
  setnames(all_events, "N", total_event_col)
  all_events[, eval(event_col) := get(total_event_col)]
  driv_event_bg <- rbind(driv_events, all_events, fill = TRUE)
  no_tum_event_gene <- paste0(event_col, "_gene")
  driv_event_bg[, eval(no_tum_event_gene) := sum(get(event_col)), by = gene_name]
  driv_event_bg[, prop := get(event_col) / get(no_tum_event_gene)]
  return(driv_event_bg)
}

### Function to melt genomic position table for cross-genome plotting:
melt_df_for_crossgenome_plot <- function(df, assets_dir){
  require(data.table)
  chromosome_lengths <- fread(file.path(assets_dir, "chromosome_lengths.csv"))
  full_genome_summary_melt <- melt.data.table(df,
                                             measure.vars = c("start", "end"),
                                             value.name = "chromosome_position")
  chromosome_lengths[, chr := as.integer(gsub("chr", "", chr))]

  full_genome_summary_melt <- merge.data.table(full_genome_summary_melt, chromosome_lengths, by = "chr", all.x = T)
  full_genome_summary_melt[, genome_position := chromosome_position + shift]
  return(full_genome_summary_melt)
}

### Function to plot cross-genome plot of gains and loss for > 1 event class:
plot_cross_genome_gain_loss_classes <- function(assets_dir,
                                                crossgenome_df,
                                                max_cn,
                                                gain_col = "frac_tums_amp",
                                                loss_col = "frac_tums_LOH",
                                                colour_vec = NULL,
                                                gain_colour_vec = NULL,
                                                loss_colour_vec = NULL,
                                                plot_title = "% tumours with event",
                                                gain_anno = NULL,
                                                loss_anno = NULL,
                                                gain_title = NULL,
                                                loss_title = NULL,
                                                legend_name = 'Event class',
                                                plot_legend = TRUE){
  require(ggplot2)
  require(cowplot)
  

  
  # melt table for plotting
  crossgenome_df <- melt_df_for_crossgenome_plot(df = crossgenome_df, assets_dir = assets_dir)
  setorder(crossgenome_df, chr, genome_position)

  classes_to_compare <- levels(crossgenome_df$event_class)
  print("Plotting cross-genome plot, comparing classes:")
  print(classes_to_compare)

  gain_col <- sym(gain_col)
  loss_col <- sym(loss_col)
  
  if (!is.null(colour_vec)) {
    gain_colour_vec = colour_vec
    loss_colour_vec = colour_vec
  }

  frac_gains <- ggplot(crossgenome_df,
                       aes(genome_position, !!gain_col * 100)) +
    geom_line(aes(group = interaction(event_class, chr), colour = event_class), alpha = .9, size = .4) +
    theme_classic() +
    geom_vline(aes(xintercept = cumsum), linetype = 2, color = "grey") +
    geom_hline(aes(yintercept = 0), alpha = .4) +
    geom_text(aes(x = chr_midpoint, label = chr), y = max_cn - 5, color = "black", data = unique(crossgenome_df[, .(chr, chr_midpoint)]), size = 5) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 17),
      axis.text.y = element_text(size = 14),
      plot.title = element_text(size = 20),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )+
    labs(x = NULL, y = "Gain", title = plot_title) +
    scale_y_continuous(limits = c(0, max_cn))+
    scale_colour_manual(values = gain_colour_vec, name = legend_name)


  frac_losses <- ggplot(crossgenome_df,
                        aes(genome_position, !!loss_col * -100)) +
    geom_line(aes(group = interaction(event_class, chr), colour = event_class), alpha = .9, size = .4) +
    theme_classic() +
    geom_vline(aes(xintercept = cumsum), linetype = 2, color = "grey") +
    geom_hline(aes(yintercept = 0), alpha = .4) +
    geom_text(aes(x = chr_midpoint, label = chr), y = -(max_cn - 5), color = "black", data = unique(crossgenome_df[, .(chr, chr_midpoint)]), size = 5) +
    theme(
      title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.title.y = element_text(size = 17),
      axis.text.y = element_text(size = 14),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )+
    labs(x = "Genome position", y = "Loss") +
    scale_y_continuous(limits = c(-max_cn, 0))+
    scale_colour_manual(values = loss_colour_vec, name = legend_name)
  
  
  ### Add axis title labels, if specified
  if (!is.null(gain_title)) {
    frac_gains <- frac_gains +
      ylab(gain_title)
  } 
  if (!is.null(loss_title)) {
    frac_losses <- frac_losses +
      ylab(loss_title)
  }
  
  ### Add gene annotations, if specified
  if (!is.null(gain_anno)) {
      for (gene in gain_anno){
        gene_start <- unique(crossgenome_df[gene_name == gene & variable == "start", genome_position])
        gene_end <- unique(crossgenome_df[gene_name == gene & variable == "end", genome_position])
        
        frac_gains <- frac_gains +
          annotate(geom = 'rect', xmin = gene_start, xmax = gene_end, ymin = 0, ymax = max_cn * 1/3 + 11, fill = 'grey') + 
          annotate('text', label = gene, x = gene_start, y = max_cn * 1/3 + 12, colour = 'grey40')
    }
  }
  if (!is.null(loss_anno)) {
      for (gene in loss_anno){
        gene_start <- unique(crossgenome_df[gene_name == gene & variable == "start", genome_position])
        gene_end <- unique(crossgenome_df[gene_name == gene & variable == "end", genome_position])
        
        frac_losses <- frac_losses +
          annotate(geom = 'rect', xmin = gene_start, xmax = gene_end, ymax = 0, ymin = -max_cn * 1/3 - 11, fill = 'grey') + 
          annotate('text', label = gene, x = gene_start, y = -max_cn * 1/3 - 12, colour = 'grey40')
    }
  }
  
  # remove legend
  if (plot_legend == FALSE) {
    frac_gains <- frac_gains + 
      theme(legend.position = "none")
    frac_losses <- frac_losses + 
      theme(legend.position = "none")
  }

  rel_heights_vec = c(8, 8)
  crossgenome_plots <- plot_grid(frac_gains, frac_losses, align="v", rel_heights=rel_heights_vec, ncol=1)

  return(crossgenome_plots)
}

### Function to plot event count vs background for seeding vs another class
plot_event_vs_bg_barplot <- function(event_df, plot_title, out_filename, colour_vec, label_col, plot_width=3, plot_height=6.5){
  # event_df <- copy(event_df)
  # order genes
  event_df[, prop_seeding := sum(prop * (event_class == "MRCA --> Seeding")), by = gene_name]
  event_df[, gene_name_sig := ""]
  event_df[fisher_pval < .05, gene_name_sig := " *"]
  event_df[, gene_name_sig := paste0(gene_name, gene_name_sig)]
  setorder(event_df, prop_seeding)

  gene_order <- unique(event_df$gene_name_sig)
  gene_order <- gene_order[gene_order != "background"]
  gene_order <- c(gene_order, "background")
  gene_face <- rep('plain', length(gene_order))
  gene_face[grepl("[ *]", gene_order)] <- 'bold'

  event_df[, gene_name_ordered := factor(gene_name_sig, levels = gene_order)]

  # Plot counts + proportions ordered
  # Proportions ordered
  label_col <- sym(label_col)
  gprop <- ggplot(event_df,
              aes(gene_name_ordered, prop, fill = event_class)) +
    geom_bar(stat = "identity", position = position_fill()) +
    geom_text(aes(label = !!label_col), position = position_fill(vjust = 0.5)) +
    scale_fill_manual(values = colour_vec, name = "Event class") +
    labs(title = plot_title, x = "Gene", y = "Proportion") +
    theme_bw() +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12),
          axis.text.y = element_text(face = gene_face),
          strip.text = element_text(size = 16),
          legend.position = "none") +
    coord_flip() +
    scale_y_continuous(breaks = c(0, .5, 1), labels = c(0, 0.5, 1)) 
  pdf(out_filename, width = plot_width, height = plot_height)
  print(gprop)
  dev.off()
}

### Function to plot a cross-genome plot but for one chromosome only, annotating gene:
plot_chromosome_gene_anno <- function(gene_list,
                                      chromosome,
                                      crossgenome_df,
                                      max_cn = 30,
                                      event_col,
                                      colour_vec,
                                      assets_dir,
                                      y_title = NULL,
                                      plot_title = NULL) {
  require(ggplot2)
  require(data.table)
  require(cowplot)
  require(scales)
  # melt table for plotting
  crossgenome_df <- melt_df_for_crossgenome_plot(df = crossgenome_df, assets_dir = assets_dir)
  setorder(crossgenome_df, chr, genome_position)

  classes_to_compare <- levels(crossgenome_df$event_class)
  print("Plotting cross-genome plot, comparing classes:")
  print(classes_to_compare)

  event_col <- sym(event_col)

  if (length(gene_list)>0){
    
    chrom_summary_df <- crossgenome_df[chr == chromosome]
    chrom_summary_df[, chr_midpoint_chrompos := chr_midpoint - shift]
    chrom_summary_df[, chromosome_position_mb := chromosome_position / 1e+06]
    
    
    g <- ggplot(chrom_summary_df, aes(chromosome_position_mb, !!event_col * 100)) +
      geom_line(aes(group = event_class, colour = event_class), alpha = .9, size = .4) +
      theme_classic() +
      scale_colour_manual(values = colour_vec, name = "Event class") +
      labs( x = paste0('Chr ', chromosome, ' (Mb)'), 
            y = y_title,
            title = plot_title
      ) +
      theme(
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 11),
        plot.title = element_text(size = 16),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "bottom"
    ) +
    guides(colour = guide_legend(nrow = 2))

    for (gene in gene_list){
      gene_start <- unique(chrom_summary_df[gene_name == gene & variable == "start", chromosome_position_mb])
      gene_end <- unique(chrom_summary_df[gene_name == gene & variable == "end", chromosome_position_mb])
      
      g <- g +
        annotate(geom = 'rect', xmin = gene_start, xmax = gene_end, ymin = 0, ymax = max_cn * 1/3 + 10, fill = 'grey50') + 
        annotate('text', label = gene, x = gene_start, y = max_cn * 1/3 + 12, colour = 'grey40')
    }
    
    return(g)  
  }
}



###############################################
### Functions to plot examples for figure 3 ###
###############################################
# plot_clones_alleles_gene <- function()



# ### 2. Copy number distance
  # # ALPACA copy number data:
  # # get copy number distance between every node on the tree:
  # # clone cn weighted by segment fraction of the genome:
  # alpaca_clone <- unique(alpaca_tum[, .(segment, clone, pred_CN_A, pred_CN_B)])
  # alpaca_clone[, pred_total_CN := pred_CN_A + pred_CN_B]
  # alpaca_clone[, chr := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][1]), segment)]
  # alpaca_clone[, start := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][2]), segment)]
  # alpaca_clone[, end := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][3]), segment)]
  # alpaca_clone[, seg_len := end - start]
  # alpaca_clone[, genome_len := sum(seg_len), by = "clone"]
  # alpaca_clone[, seg_fraction := seg_len / genome_len]
  # # Compute copy number distance between each pair of clones on the tree
  # cn_dist_tum <- pairwise_cn_distance_alpaca(alpaca_clone)
  # # process cn distance for every edge of the tree
  # tree_df[, c("l1_cn_branchlength_A", "l1_cn_branchlength_B", "l1_mean_cn_branchlength", "l2_cn_branchlength_A", "l2_cn_branchlength_B", "l2_mean_cn_branchlength") := rbindlist(lapply(clone, function(child) {
  #   parent <- tree_structure[tree_structure[, 2] == child, 1]
  #   cn_dist <- cn_dist_tum[(clone1 == parent & clone2 == child) | (clone1 == child & clone2 == parent),
  #     .(l1_cn_distance_A, l1_cn_distance_B, l1_mean_dist_allele, l2_cn_distance_A, l2_cn_distance_B, l2_mean_dist_allele)]
  #   if (nrow(cn_dist) == 0) cn_dist <- data.table(cn_distance_A = NA, cn_distance_B = NA, mean_dist_allele = NA)
  #   return(cn_dist)
  # }))]
  # # add cn distance for every node compared to MRCA
  # tree_df[, c("l1_cn_distance_MRCA_A", "l1_cn_distance_MRCA_B", "l1_mean_cn_distance_MRCA", "l2_cn_distance_MRCA_A", "l2_cn_distance_MRCA_B", "l2_mean_cn_distance_MRCA") := rbindlist(lapply(clone, function(child) {
  #   cn_dist <- cn_dist_tum[(clone1 == trunk & clone2 == child) | (clone1 == child & clone2 == trunk),
  #     .(l1_cn_distance_A, l1_cn_distance_B, l1_mean_dist_allele, l2_cn_distance_A, l2_cn_distance_B, l2_mean_dist_allele)]
  #   if(nrow(cn_dist) == 0) {
  #     cn_dist <- data.table(l1_cn_distance_A = NA, l1_cn_distance_B = NA, l1_mean_dist_allele = NA, l2_cn_distance_A = NA, l2_cn_distance_B = NA, l2_mean_dist_allele = NA)
  #   }
  #   if (child == trunk) {
  #     cn_dist <- data.table(l1_cn_distance_A = 0, l1_cn_distance_B = 0, l1_mean_dist_allele = 0, l2_cn_distance_A = 0, l2_cn_distance_B = 0, l2_mean_dist_allele = 0)
  #   }
  #   return(cn_dist)
  # }))]
  # # melt cn distance dataframes by allele
  # l1_cn_distance_df <- melt.data.table(tree_df[, .(clone, l1_cn_branchlength_A, l1_cn_branchlength_B, l1_mean_cn_branchlength)], measure.vars = c("l1_cn_branchlength_A", "l1_cn_branchlength_B"), variable.name = "allele", value.name = "l1_cn_branchlength")
  # l2_cn_distance_df <- melt.data.table(tree_df[, .(clone, l2_cn_branchlength_A, l2_cn_branchlength_B, l2_mean_cn_branchlength)], measure.vars = c("l2_cn_branchlength_A", "l2_cn_branchlength_B"), variable.name = "allele", value.name = "l2_cn_branchlength")
  # l1_cn_distance_MRCA_df <- melt.data.table(tree_df[, .(clone, l1_cn_distance_MRCA_A, l1_cn_distance_MRCA_B, l1_mean_cn_distance_MRCA)], measure.vars = c("l1_cn_distance_MRCA_A", "l1_cn_distance_MRCA_B"), variable.name = "allele", value.name = "l1_cn_distance_MRCA")
  # l2_cn_distance_MRCA_df <- melt.data.table(tree_df[, .(clone, l2_cn_distance_MRCA_A, l2_cn_distance_MRCA_B, l2_mean_cn_distance_MRCA)], measure.vars = c("l2_cn_distance_MRCA_A", "l2_cn_distance_MRCA_B"), variable.name = "allele", value.name = "l2_cn_distance_MRCA")
  # l1_cn_distance_df[, allele := gsub("l1_cn_branchlength_", "", allele)]
  # l2_cn_distance_df[, allele := gsub("l2_cn_branchlength_", "", allele)]
  # l1_cn_distance_MRCA_df[, allele := gsub("l1_cn_distance_MRCA_", "", allele)]
  # l2_cn_distance_MRCA_df[, allele := gsub("l2_cn_distance_MRCA_", "", allele)]
