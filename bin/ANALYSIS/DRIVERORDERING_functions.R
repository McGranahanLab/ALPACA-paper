#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

#####################################
### BRADLEY TERRY MODEL FUNCTIONS ###
#####################################

### Function to generate driver muttable
process_driver_muttable_for_btm <- function(muttable, gene_drivclass) {
  driver_muttab <- muttable[DriverMut == TRUE]
  driver_muttab[, gene_name := Hugo_Symbol]

  driver_muttab <- merge.data.table(driver_muttab, gene_drivclass[!is.na(mut_driver_class), .(gene_name, is_lung_mut_driver, mut_driver_class)], all = FALSE)
  
  # summarise all driver mutations present per clone
  driver_muttab[, mut_drivers_clone := paste(sort(unique(Hugo_Symbol)), collapse = ";"), by = .(tumour_id, PyCloneCluster_SC)]
  # driver_muttab_sum <- unique(driver_muttab[!is.na(PyCloneCluster_SC), .(tumour_id, clone = as.character(PyCloneCluster_SC), mut_drivers_clone)])

  driver_muttab[, n_tums_with_event := length(unique(tumour_id)), by = Hugo_Symbol]
  driver_muttab_small <- unique(driver_muttab[, .(tumour_id,
                                                  cluster_with_driver = PyCloneCluster_SC,
                                                  cluster_with_driver_clonality = PyCloneClonal_SC,
                                                  gene_name = Hugo_Symbol,
                                                  is_lung_mut_driver,
                                                  mut_driver_class,
                                                  n_tums_with_event,
                                                  mut_id)]
  )
  return(driver_muttab_small)
}

### Function to select events based on tumour frequency 
get_event_frequency_fromlist <- function(selected_events_list, min_n_tums_w_event) {
  all_events <- rbindlist(lapply(seq(selected_events_list), function(i) {
    event_list <- selected_events_list[[i]]
    name <- event_list$name
    df <- copy(event_list$df)
    full_df <- copy(event_list$full_df)
    n_events <- event_list$n_events
    cat('\nName: ', name)
    cat('\nNo events: ', n_events, '\n')
    df[, n_tums_with_event := length(unique(tumour_id)), by = .(event_id)]
    df_summary <- unique(df[, .(event_id, n_tums_with_event)])
    setorder(df_summary, -n_tums_with_event)
    top_events <- df_summary[n_tums_with_event >= min_n_tums_w_event][1:n_events, event_id] 
    full_df[, top_event := event_id %in% top_events]
    full_df[, event_name := name]
    full_df[, n_top_events := n_events]
    full_df <- merge.data.table(full_df, df_summary, by = "event_id", all = TRUE)
    setorder(full_df, -n_tums_with_event)
    return(full_df)
  }))
  all_events <- unique(all_events)
  return(all_events)
}


# Function to remove a vector of clones from a tree edge matrix
remove_nodes_from_tree <- function(edge_mat, nodes_to_remove){

  new_edge_mat <- copy(edge_mat)
  for (node in nodes_to_remove) {
    # print(node)
    parent_of_node <- new_edge_mat[clone == node, parent_clone]
    children_of_node <- new_edge_mat[parent_clone == node, clone]
    # remove edge connecting parent of node to node
    new_edge_mat <- new_edge_mat[!(parent_clone == parent_of_node & clone == node)]
    
    # remove edge connecting node to children of node
    new_edge_mat <- new_edge_mat[!(parent_clone == node & clone %in% children_of_node)]
    
    # add new edges connecting parent of node to children of node
    if (length(children_of_node) > 0) {
        new_edges <- data.table(parent_clone = parent_of_node,
                                clone = children_of_node)
        new_edge_mat <- rbind(new_edge_mat, new_edges)
    }
  }
  return(new_edge_mat)
}

### Function to get all ancestor - descendant relationships from a tree:
get_ancestor_descendant_relationships <- function(tree_graph) {
  all_clones <- unique(c(as.matrix(tree_graph)))
  all_relationships <- rbindlist(lapply(all_clones, function(c) {
    clone_ancestors <- get_cluster_ancestors(tree_graph, c, include_trunk = TRUE, include_cluster = FALSE)
    if (length(clone_ancestors) > 0) {
      return(data.table(expand.grid(clone_ancestors, c)))
    }
  }))
  colnames(all_relationships) <- c("ancestor", "descendant")
  return(all_relationships)
}


# Function to get all ancestor and descendant driver alterations for a tumour
# previously: get_driver_edges_tumour
get_driver_ancestor_descendant_pairs_tumour <- function(tum_events_df) {
  # ensure tum_events_df has the columns:
  # clone, parent_clone, driver_alterations_clone
  tum_events_df[, clone := as.character(clone)]
  tum_events_df[, parent_clone := as.character(parent_clone)]

  clones_without_drivers <- tum_events_df[is.na(driver_alterations_clone), clone]

  # reduce tree to only nodes with drivers
  new_edge_mat <- remove_nodes_from_tree(edge_mat = tum_events_df[, .(parent_clone, clone)],
                                          nodes_to_remove = clones_without_drivers
  )

  # Get every ancestor -> descendant clone pair
  anc_decs <- get_ancestor_descendant_relationships(new_edge_mat)

  descendant_drivers <- tum_events_df[clone %in% new_edge_mat$clone, .(descendant = clone, descendant_drivers = driver_alterations_clone)]
  descendant_drivers <- merge.data.table(descendant_drivers, anc_decs, by = "descendant", all = TRUE)
  ancestor_drivers <- tum_events_df[clone %in% new_edge_mat$parent_clone, .(ancestor = clone, ancestor_drivers = driver_alterations_clone)]

  # create driver edges
  driver_edges <- merge.data.table(descendant_drivers, ancestor_drivers, by = "ancestor", all = TRUE)
  # remove diploid -> trunk
  driver_edges <- driver_edges[ancestor != "diploid"]
  # remove cases where the ancestor == NA (only relevant when restricting to subclonal events only)
  driver_edges <- driver_edges[!is.na(ancestor_drivers)]
  
  # melt to long format, each row a new ancestor -> descendant driver alteration
  if (nrow(driver_edges) > 0) {
    driver_edges_long <- rbindlist(lapply(seq(nrow(driver_edges)), function(i) {
      ancestor_drivers <- strsplit(driver_edges[i, ancestor_drivers], ";")[[1]]
      descendant_drivers <- strsplit(driver_edges[i, descendant_drivers], ";")[[1]]
      edge_driver_order <- as.data.table(expand.grid(ancestor_drivers, descendant_drivers))
      colnames(edge_driver_order) <- c('ancestor_event', 'descendant_event')
      edge_driver_order <- unique(edge_driver_order)
      edge_driver_order <- cbind(driver_edges[i, .(ancestor, descendant)], edge_driver_order)
      return(edge_driver_order)
    }))
    driver_edges_long[, tumour_id := unique(tum_events_df$tumour_id)]
    driver_edges_long[, ancestor_event := as.character(ancestor_event)]
    driver_edges_long[, descendant_event := as.character(descendant_event)]
    return(driver_edges_long)
  } else {
    return(NULL)
  }
}

### TODO: ###
# Add function to get driver ancestor descendant pairs, taking into account:
# - truncal loh/wgd 
# - ties within clones
# get_driver_ancestor_descendant_pairs_tumour_new <- function(tum_distance_df, add_draws = TRUE, add_truncal_loh = TRUE)


### Function to get frequency matrix of ancestor descendant pairs:
get_pairwise_event_freq <- function(driver_binary_df) {
  driver_edge_data <- copy(driver_binary_df)
  all_unique_alterations <- sort(unique(c(driver_binary_df$ancestor, driver_binary_df$descendant)))

  driver_edge_data[, n_anc_dec := .N, by = c('ancestor', 'descendant')]
  n_anc_dec_df <- unique(driver_edge_data[, .(ancestor, descendant, n_anc_dec)])
  
  all_pairs <- data.table(expand.grid(all_unique_alterations, all_unique_alterations))
  colnames(all_pairs) <- c('ancestor', 'descendant')
  all_pairs_counts <- merge.data.table(all_pairs, n_anc_dec_df, by = c('ancestor', 'descendant'), all = TRUE)
  all_pairs_counts[is.na(n_anc_dec), n_anc_dec := 0]

  # now get number of co-occurrences:
  all_pairs_counts[, pair := paste(sort(c(as.character(ancestor), as.character(descendant))), collapse = ';'), by = seq(nrow(all_pairs_counts))]
  all_pairs_counts[, n_cooccurrences := sum(n_anc_dec), by = pair]

  # now get fraction of co-occurrences:
  all_pairs_counts[, frac_anc_dec := n_anc_dec / n_cooccurrences]
  return(all_pairs_counts)
}

# get frequency counts of pairwise ancestor-descendant pairs
get_driver_edge_counts <- function(driver_binary_df) {
  driver_edge_data <- copy(driver_binary_df)
  # Preprocess columns to get unique pairs of drivers (for each 'match played')
  driver_edge_data[, ancestor := as.character(ancestor)]
  driver_edge_data[, descendant := as.character(descendant)]
  driver_edge_data[, n_anc_dec := .N, by = c('ancestor', 'descendant')]

  all_unique_alterations <- unique(c(driver_edge_data$ancestor,
                                      driver_edge_data$descendant))
  
  print("Preprocessing input data")
  driver_edge_data[, pair := paste(sort(c(ancestor, descendant)), collapse = ';'), by = seq(nrow(driver_edge_data))]
  driver_edge_data[, driver1 := gsub(";.*", "", pair)]
  driver_edge_data[, driver2 := gsub(".*;", "", pair)]
  driver_edge_data[, which_event_ancestor := ifelse(driver1 == ancestor,
  "driver1_ancestor", "driver2_ancestor")
  ]
  
  # Count number of instances where ancestor -> descendant, to use as input to bradley terry model
  driver_edge_counts <- dcast(driver_edge_data,
                              pair + driver1 + driver2 ~ which_event_ancestor,
                              fun.aggregate = length)
  setorder(driver_edge_counts, pair)

  # remove cases where ancestor == descendant
  driver_edge_counts <- driver_edge_counts[driver1 != driver2]

  # make correct column names:
  driver_edge_counts_df <- driver_edge_counts[, .(driver1,
                                                  driver2,
                                                  pair,
                                                  win1 = driver1_ancestor,
                                                  win2 = driver2_ancestor)]
  driver_edge_counts_df[, driver1 := factor(driver1, levels = all_unique_alterations)]
  driver_edge_counts_df[, driver2 := factor(driver2, levels = all_unique_alterations)]   
  
  return(driver_edge_counts_df)
}

# Function to run Bradley Terry model
run_btm <- function(driver_binary_df, all_unique_alterations) {
  suppressPackageStartupMessages(require(BradleyTerry2))
  driver_edge_counts_df <- get_driver_edge_counts(driver_binary_df)

  # Run Bradley-Terry model (from Rpackage BradleyTerry2)
  print("Running Bradley Terry model")
  driver_order_model <- BTm(cbind(win1, win2),
                            driver1,
                            driver2,
                            formula = ~ driver,
                            id = "driver",
                            data = driver_edge_counts_df
  )
  
  print("Computing bias-reduced estimate")
  # Update by doing bias-reduced estimate
  driver_order_model_br <- update(driver_order_model, br = TRUE)
  bt_ordering <- as.data.frame(summary(driver_order_model_br)$coefficients)
  bt_ordering$alteration <- rownames(bt_ordering)
  bt_ordering <- as.data.table(bt_ordering)
  bt_ordering[, gene := gsub("driver", "", alteration)]
  return(bt_ordering)
}


### Function to get the lower triangle matrix to plot 
# Input data is taken from output of the function get_pairwise_event_freq()
get_lower_tri_pairwise <- function(top_events_pairs_df, top_events_prevalence_order) {
  top_events_pairs_ordered <- copy(top_events_pairs_df)
  top_events_pairs_ordered[, ancestor := factor(as.character(ancestor), levels = top_events_prevalence_order)]
  top_events_pairs_ordered[, descendant := factor(as.character(descendant), levels = top_events_prevalence_order)]
  pair_mat <- dcast(top_events_pairs_ordered, formula = ancestor ~ descendant, value.var = 'n_anc_dec')
  pair_mat <- as.data.frame(pair_mat)
  rownames(pair_mat) <- pair_mat[, 1]
  pair_mat[, 1] <- NULL
  pair_mat <- as.matrix(pair_mat)
  lowertri <- lower.tri(pair_mat, diag = TRUE)
  pair_mat[lowertri] <- NA
  lowertri_names <- as.data.frame(which(is.na(pair_mat), arr.ind = TRUE))
  pairs_to_keep <- data.table(row = rownames(pair_mat)[lowertri_names$row], col = colnames(pair_mat)[lowertri_names$col])
  pairs_to_keep[, pair := paste(c(as.character(row), as.character(col)), collapse = ';'), by = seq(nrow(pairs_to_keep))]
  return(pairs_to_keep)
}
