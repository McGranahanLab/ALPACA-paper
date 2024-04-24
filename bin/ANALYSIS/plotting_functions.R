### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(igraph))

# Function to plot a tree
plot_tree <- function(tree_edges,
                      vertex_size = 30,
                      vertex_label_cex = 1,
											vertex_frame_width = NULL,
                      colour_vector = NULL,
                      label_tree = TRUE,
                      plot_title = '',
                      title_size = 5) {
  require(igraph)
  tree_edges <- apply(tree_edges, c(1, 2), as.character)
  tree_edges <- as.data.frame(tree_edges)
  g <- graph.data.frame(tree_edges, directed = FALSE)
  v_index <- V(g)$name
  
  if (!is.null(colour_vector)) {
    vcol <- colour_vector[v_index]
  } else vcol <- "white"
	
	if (!is.null(vertex_frame_width)) {
		vbd <- vertex_frame_width[v_index]
	} else vbd <- 2


  if (length(vertex_size) > 1) {
    vsize <- vertex_size[v_index]
  } else {
    vsize <- vertex_size
  }
  

  trunk <- unique(tree_edges[, 1][!tree_edges[, 1] %in% tree_edges[, 2]])
  l <- layout_as_tree(g, root = trunk)
  
  if (label_tree) {
    vertex_label <- V(g)$name
  } else {
    vertex_label <- rep("", length(colour_vector))
  }


  # pdf(file.path(tum_out_dir, "treeplot_seeding.pdf"), width = 5, height = 5)
  plot.igraph(g
              , layout = l
              , vertex.color = vcol
              , vertex.size = vsize
              , vertex.frame.color = "black"
              , vertex.lwd = 5
              , edge.color = "black"
              , vertex.label = vertex_label
              , vertex.label.cex = vertex_label_cex
							, vertex.frame.width = 1
              , vertex.label.pos = 2
              , vertex.label.family = "Helvetica"
              , vertex.label.font = 2
              , vertex.label.color = "black"
              , main = plot_title
              , cex.main = title_size)
	# dev.off()
}



# Function to colour the tree using TRACERx colours,
# based on tree building output t$graph_pyclone$edgelength
colour_tree_tx <- function(node_vec, opacity = 255) {
  #opacity: transparency. Nicky's default is 153
  suppressPackageStartupMessages(require(RColorBrewer))
  ncols <- length(node_vec)
  max.cols <- 12
  cols     <- paste(brewer.pal(min(max.cols, ncols), name = "Paired"), sep = "")
  cols     <- rep(cols, ceiling(ncols / max.cols))[1:ncols]
  #Make colors with more than 12 clusters more transparent
  nopaq <- round(opacity / ceiling(ncols / 12))
  opacity_val <- c(
        sapply(1:ceiling(ncols / 12),
        function(x) rep(opacity - nopaq * (x - 1), 12))
        )[1:ncols]
  cols.opac <- paste(cols, as.hexmode(opacity_val), sep = "")
  names(cols.opac) <- node_vec
  return(cols.opac)
}

# Function to colour the tree with new RColorBrewer colours for each new clone
colour_tree_clonemap <- function(node_vec, opacity = 255) {
  suppressPackageStartupMessages(require(cloneMap))
  suppressPackageStartupMessages(require(RColorBrewer))
  brewerpals <- c("Paired", "Set3", "Accent", "Dark2",
                        "Pastel2", "Set1", "Set2", "Pastel1", "Spectral")

  ncols <- length(node_vec)
  chosenpals <- brewer.pal.info[brewerpals, ]
  chosenpals$cum_ncols <- sapply(
		seq(nrow(chosenpals)),
		function(i) {
						sum(chosenpals[1:i, "maxcolors"])
	})
  min_npals <- which(chosenpals$cum_ncols == min(chosenpals$cum_ncols[chosenpals$cum_ncols >= ncols]))
  col_pals <- chosenpals[1:min_npals, , drop = FALSE]
  final_pal_req <- rownames(col_pals)[nrow(col_pals)]
  if (nrow(col_pals) > 1) {
		n_finalpal <- ncols - col_pals[nrow(col_pals) - 1, "cum_ncols"]
  } else {
		n_finalpal <- ncols
  }

  cols <- unlist(lapply(rownames(col_pals), function(cur_pal) {
		if (cur_pal == final_pal_req) colour_vec <- brewer.pal(n_finalpal, cur_pal)
		else colour_vec <- brewer.pal(col_pals[cur_pal, "maxcolors"], cur_pal)
		return(colour_vec)
  }))
  names(cols) <- node_vec
  return(cols)
}

# Function to plot all output for one tumour
plot_tree_clonemaps <- function(tum_id,
                                tree_structure,
                                ccf_cluster_table,
                                node_size_vec,
                                colour_vec,
                                tree_out_dir,
                                vertex_size = 15,
                                label_tree = FALSE) {
  suppressPackageStartupMessages(require(cloneMap))
  vertex_size_vec <- 15 * log(node_size_vec)

  ### Plot tree individually
  print("Plotting tree:")
  pdf(file.path(tree_out_dir, "treeplot.pdf"), width = 10, height = 10)
  plot_tree(tree_edges = tree_structure,
          # vertex_size = vertex_size_vec,
          vertex_size = vertex_size,
          vertex_label_cex = 1.2,
          colour_vector = colour_vec,
          label_tree)
  dev.off()

  ### Plot each clone map individually
  print("Plotting clone maps:")
  samples <- colnames(ccf_cluster_table)
  for (r in samples) {
    print(r)
    reg_ccf <- as.data.frame(ccf_cluster_table[, r, drop = FALSE])
    reg_ccf$clones <- rownames(reg_ccf)
    colnames(reg_ccf) <- c("CCF", "clones")
    reg_ccf <- reg_ccf[, c("clones", "CCF")]
    cm_out <- paste0(tree_out_dir, "/clonemap_", r, ".pdf")
    if (!file.exists(cm_out)) {
      pdf(cm_out,
          width = 4,
          height = 4,
          bg = "transparent")
      cloneMap(tree.mat = tree_structure,
              CCF.data = reg_ccf,
              border.thickness = 2.2,
              clone.cols = colour_vec)
      dev.off()
    }
  }
  

  ### Plot tree and clonemaps arranged on a pdf
  print("Plotting tree and clone maps together:")
  tree_cm_out <- paste0(tree_out_dir, "/", tum_id, "_tree_clonemap.pdf")
  if (!file.exists(tree_cm_out)) {
    nregs <- length(samples)
    nrows <- 3
    ncols <- 1 + ceiling(nregs / 3)
    remainder <- 3 * (ceiling(nregs / 3) - nregs / 3)
    nzeros <- remainder %% 3
    plot_positions <- c(1, 1, 1, 2:(nregs + 1))
    pdf(tree_cm_out,
            width = 5 + 5 * (ceiling(nregs / 3) / 3),
            height = 5)
    layout_matrix <- matrix(0L, nrow = nrows, ncol = ncols)
    layout_matrix[seq(plot_positions)] <- plot_positions
    layout(mat = layout_matrix,
            heights = rep(3, nrows),
            widths = c(3, rep(1, (ncols - 1)))
    )
    # layout.show(n=(nregs+1))
    par(mar = c(0, 0, 0, 0))
    plot_tree(tree_edges = tree_structure,
            # vertex_size = vertex_size_vec,
            vertex_size = vertex_size,
            vertex_label_cex = 1.7,
            colour_vector = colour_vec,
            label_tree
    )
    for (r in samples) {
      print(r)
      par(mar = c(0, 0, 0, 0))
      reg_ccf <- as.data.frame(ccf_cluster_table[, r, drop = FALSE])
      reg_ccf$clones <- rownames(reg_ccf)
      colnames(reg_ccf) <- c("CCF", "clones")
      reg_ccf <- reg_ccf[, c("clones", "CCF")]
      cloneMap(tree.mat = tree_structure,
              CCF.data = reg_ccf,
              border.thickness = 2.2,
              clone.cols = colour_vec)
    }
    dev.off()
  }
}


# Function to plot seeding tree coloured by seeding clone paths
colour_seeding_tree <- function(tree_paths, tum_mets_cloneInfo) {
	
	# Get clone colouring order, by tree level
	tree_paths <- lapply(tree_paths, function(p) gsub("clone", "", p))
	tree_structure <- convert_tree_paths_to_edgemat(tree_paths)
	all_clones <- sort(unique(as.numeric(c(tree_structure))))

	col_df <- data.table(clone = all_clones, colour = '')
	col_df[, tree_level := sapply(clone, function(c){
		ancscs <- get_cluster_ancestors(tree_structure, c, include_trunk = TRUE, include_cluster = TRUE)
    if (length(ancscs) == 1) return(0)
    else return(length(ancscs) - 1)
	})]
	
	col_order <- merge.data.table(col_df, tum_mets_cloneInfo, by = "clone", all.x = TRUE, all.y = FALSE)
	setorder(col_order, tree_level)

	### Assign shared clones
	shared_col <- "grey"
	col_order[sharedClones == TRUE, colour := shared_col]

	### Assign MRCA
	mrca_col <- "#4F7B4F"
	col_order[clonalClones == TRUE, colour := mrca_col]
	
	### Assign seeding clones & their descendants
	seeding_clones_order <- col_order[seedingClones == TRUE, clone]
	seeding_colours <- data.table(seeding_clone = c(1:3), seeding_start = c("#855593", "#08519C", "#A50F15"))#,# seeding_end = c("#DAC0E1", "#F7FBFF", "#FFF5F0"))

	seeding_col_no <- 1
	for (s in seeding_clones_order) {
		# print(s)
		descendants <- get_cluster_descendants(tree_structure, s, include_clone = TRUE)
		descendants_order <- as.character(col_order[clone %in% descendants, clone])
		start_col <- seeding_colours[seeding_clone == seeding_col_no, seeding_start]
		descendants_colour_vec <- colorRampPalette(c(start_col, "white"))(length(descendants_order) + 1)[1:length(descendants_order)]
		names(descendants_colour_vec) <- descendants_order
		for (d in descendants_order) {
			col_order[clone == d, colour := descendants_colour_vec[as.character(d)]]
		}
		seeding_col_no <- seeding_col_no + 1
	}

	### Assign (remaining) non-seeding primary clones & their descendants	
	primary_clones_order <- col_order[primaryClones == TRUE & colour == '', clone]
	primary_col <- "#74C476"#), primary_end = c("#BCD5BC"))

	primary_colour_vec <- colorRampPalette(c(primary_col, "white"))(length(primary_clones_order) + 1)[1:length(primary_clones_order)]
	names(primary_colour_vec) <- as.character(primary_clones_order)
	for (p in primary_clones_order) {
		col_order[clone == p, colour := primary_colour_vec[as.character(p)]]
	}

	colour_vec <- col_order[, colour]
	names(colour_vec) <- col_order[, clone]
	return(colour_vec)
}


# Function to plot a tree with seeding colours
plot_seeding_tree <- function(tree_paths, tum_mets_cloneInfo, tum_out_dir) {
	tree_structure <- convert_tree_paths_to_edgemat(tree_paths)
	colour_vec <- colour_seeding_tree(tree_paths, tum_mets_cloneInfo)
	node_bd_vec <- sapply(names(colour_vec), function(x) {
		size <- ifelse(x %in% tum_mets_cloneInfo[seedingClones == TRUE, clone], 10, 2)
		return(size)
	})
	pdf(file.path(tum_out_dir, "treeplot_seeding.pdf"), width = 5, height = 5)
	plot_tree(tree_edges = tree_structure,
					# vertex_size = vertex_size_vec,
					vertex_size = 30,
					vertex_label_cex = 1.7,
					vertex_frame_width = node_bd_vec,
					colour_vector = colour_vec,
					label_tree = TRUE)
	dev.off()
  # Also save in figures/ directory
  tum_id <- basename(tum_out_dir)
  pdf(paste0('figures/Fig2_', tum_id, "_treeplot_seeding.pdf"), width = 5, height = 5)
	plot_tree(tree_edges = tree_structure,
					# vertex_size = vertex_size_vec,
					vertex_size = 30,
					vertex_label_cex = 1.7,
					vertex_frame_width = node_bd_vec,
					colour_vector = colour_vec,
					label_tree = TRUE)
	dev.off()
}

# Plot each clone map individually
plot_region_clonemaps <- function(tree_paths, cp_table, tum_out_dir, colour_vec) {
	tree_structure <- convert_tree_paths_to_edgemat(tree_paths)
	cp_table <- as.data.frame(cp_table)
	rownames(cp_table) <- gsub("clone", "", cp_table$clone)
	cp_table$clone <- NULL
	cp_matrix <- as.matrix(cp_table)
	ccf_cluster_table <- cp_to_ccf(tree_structure, cp_matrix)
	print("Plotting clone maps:")
	samples <- colnames(ccf_cluster_table)
	for (r in samples) {
		print(r)
		reg_ccf <- as.data.frame(ccf_cluster_table[, r, drop = FALSE])
		reg_ccf$clones <- rownames(reg_ccf)
		colnames(reg_ccf) <- c("CCF", "clones")
		reg_ccf <- reg_ccf[, c("clones", "CCF")]
		cm_out <- paste0(tum_out_dir, "/clonemap_", r, ".pdf")
		# if (!file.exists(cm_out)) {
		pdf(cm_out, width = 4, height = 4, bg = "transparent")
		cloneMap(tree.mat = tree_structure, CCF.data = reg_ccf, border.thickness = 2.2, clone.cols = colour_vec)
		dev.off()
		# }
	}
}


### Function to read in asas tables for example tumour cases in Figure 2
loadASAS <- function(asas_path, tree_region_IDs=NULL){
  # This needs to be run after loading tree
  asas_patient <- readRDS(asas_path)
  # rbind asas tables from each region
  asas_patient_dt <- rbindlist(lapply(names(asas_patient), function(region){
      r <- asas_patient[[region]]
      segs <- rbindlist(lapply(names(r), function(seg){
          seg_list <- r[[seg]]
          seg_df <- as.data.table(seg_list)
          seg_df[, seg_name := seg]
      }), fill = T)
      segs[, sample := region]
  }),
  fill = T)
  # filter for NAs
  asas_patient_dt <- asas_patient_dt[!is.na(ph_cpnA_vec)]
  asas_patient_dt[, region := gsub('-', '.', sample)]
  if (!is.null(tree_region_IDs)) asas_patient_dt <- asas_patient_dt[region %in% tree_region_IDs]
  return(asas_patient_dt)
}


### Functions
plot_samples_alleles_gene <- function(tum_id,
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
                                      plot_height) {

  cruk_id <- substr(tum_id, 1, 8)

  # Load fractional CN
  input_df <- fread(file.path(tum_out_dir, 'ALPACA_input_table.csv'))
  input_df[, startpos := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][2]), segment)]
  input_df[, endpos := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][3]), segment)]
  mphase_selected <- input_df[chr == chromosome & startpos > min_startpos & endpos < max_endpos]
  segs_to_plot <- unique(mphase_selected$segment)
  mphase_melt <- melt.data.table(mphase_selected, measure.vars = c("cpnA", "cpnB"), variable.name = "allele", value.name = "cpn")
  mphase_melt[, allele := gsub("cpn", "", allele)]
  mphase_melt[, max_cn_chr := max(cpn), by = .(sample, chr)]
  mphase_melt[, sample_label := sample_levels[sample]]
  mphase_melt[, sample_label := factor(sample_label, levels = sample_levels)]

  # Load gene info
  tum_seg_gene <- fread(tum_seg_gene_path)
  gene_df <- tum_seg_gene[gene_name %in% gene_list]
  gene_info <- merge.data.table(gene_df[, .(gene_name, segment, gene_start = query_start, gene_end = query_end)],
                                unique(mphase_melt[segment %in% gene_df$segment, .(segment, sample_label, max_cn_chr)]),
                                by = "segment",
                                all = TRUE,
                                allow.cartesian = TRUE)
  gene_info[, label_height := ceiling(max_cn_chr) + 1]
  gene_info[, rect_height := ceiling(max_cn_chr) + .5]
  gene_info[, gene_label := gene_name]
  gene_info[sample_label != sample_levels[1], gene_label := NA]


  # PLOT
  print("Plotting chromosome fractional CN")
  g <- ggplot(mphase_melt) 
  
  # Load SNPs
  if (file.exists(asas_path)) {
    asas_df <- fread(asas_path)
    asas_df[, segment := gsub(":", "_", seg_name)]
    asas_table_seg <- asas_df[segment %in% segs_to_plot]
    asas_table_seg[, sample_label := region]
    asas_table_seg[, sample_label := sample_levels[sample_label]]
    asas_table_seg[, sample_label := factor(sample_label, levels = sample_levels)]
    asas_table_melt <- melt.data.table(asas_table_seg, measure.vars = c("ph_cpnA_vec", "ph_cpnB_vec"), variable.name = "allele", value.name = "snp_cn")
    asas_table_melt[, allele := gsub("ph_cpn|_vec", "", allele)]
    # cap at 0
    asas_table_melt[snp_cn < 0, snp_cn := 0]
    # Cap max snp cpn per segment
    asas_table_melt <- merge.data.table(asas_table_melt, unique(gene_info[, .(sample_label, max_cn_chr)]), by = "sample_label", all = TRUE)
    asas_table_melt[, max_snp_cn := max_cn_chr + 2]
    asas_table_melt <- asas_table_melt[pos_vec < max_endpos & snp_cn < max_snp_cn]

    g <- g + geom_point(data = asas_table_melt, 
              aes(x = pos_vec, y = snp_cn, colour = allele), alpha = .3, size = .8)
  }
  g <- g +
    geom_rect(aes(xmin = gene_start, xmax = gene_end, ymin = 0, ymax = rect_height), fill = "black", colour = "black", data = gene_info) +
    geom_text_repel(aes(x = gene_start, y = label_height, label = gene_label), direction = "y", data = gene_info[sample_label == sample_levels[1]], size = 5) +
    geom_segment(aes(x = startpos, xend = endpos, y = cpn, yend = cpn, colour = allele), linewidth = 2.5) +
    facet_wrap(~sample_label, ncol = 1, strip.position = "right", scales = "free_y") +
    scale_colour_manual(values = c(A = "red", B = "blue")) +
    theme_cowplot(font_family = "Helvetica") +
    theme(strip.background = element_rect(fill = "white", colour = "black"),
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 17),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_blank(),
          plot.title = element_text(size = 20, face = "bold")) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks = xaxis_breaks, labels = xaxis_breaks) +
    labs(y = "Observed fractional CN", x = paste0("Chr ", chromosome," (Mb)"))
  pdf(file.path(tum_out_dir, paste0(tum_id, "_mphase_", chromosome, ".pdf")), width = plot_width, height = plot_height)
  print(g)
  dev.off()
  # Also save in figures/ directory:
  pdf(paste0("figures/Fig2_", tum_id, "_mphase_", chromosome, ".pdf"), width = plot_width, height = plot_height)
  print(g)
  dev.off()
}


plot_clones_alleles_gene <- function(tum_id,
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
                                      plot_height) {

  all_clones <- sort(as.numeric(unique(c(clone_levels))))
  # ALPACA output
  alpaca[, chr := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][1]), segment)]
  alpaca[, startpos := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][2]), segment)]
  alpaca[, endpos := mapply(function(x) as.numeric(strsplit(x, "_")[[1]][3]), segment)]
  cn_selected <- alpaca[chr == chromosome & startpos > min_startpos & endpos < max_endpos]
  segs_to_plot <- unique(alpaca$segment)
  cn_selected[, clone := gsub("clone", "", clone)]
  clone_chr_cn_selected <- cn_selected[clone %in% clone_levels]
  clone_chr <- melt.data.table(clone_chr_cn_selected, measure.vars = c("pred_CN_A", "pred_CN_B"), variable.name = "allele", value.name = "pred_CN")
  clone_chr[, allele := gsub("pred_CN_", "", allele)]
  clone_chr[, pred_CN := as.numeric(pred_CN)]
  clone_chr[allele == "A", pred_CN := as.numeric(pred_CN) + .1]
  clone_chr[allele == "B", pred_CN := as.numeric(pred_CN) - .1]
  clone_chr[pred_CN < 0, pred_CN := 0]
  clone_chr[, max_pred_CN := max(pred_CN), by = clone]
  clone_chr[, clone_fac := factor(clone, levels = clone_levels)]
  clone_colour <- colour_vec[as.character(clone_levels)]

  # Load gene info
  tum_seg_gene <- fread(tum_seg_gene_path)
  gene_df <- tum_seg_gene[gene_name %in% gene_list]
  gene_info_clone <- merge.data.table(gene_df[, .(gene_name, segment, gene_start = query_start, gene_end = query_end)], unique(clone_chr[, .(segment, clone_fac, max_pred_CN)]), by = "segment", all = FALSE)
  gene_info_clone[, gene_label := gene_name]
  gene_info_clone[clone_fac != clone_levels[1], gene_label := NA]
  
  print(clone_chr[segment %in% gene_info_clone$segment, .(segment, clone_fac, pred_CN)])
  
  # Plot
  g <- ggplot(clone_chr) +
    facet_wrap(~clone_fac, ncol = 1, strip.position = "right", scales = "free_y") +
    geom_rect(aes(xmin = gene_start, xmax = gene_end, ymin = 0, ymax = max_pred_CN + 2), fill = "black", colour = "black", data = gene_info_clone) +
    theme_cowplot() +
    geom_hline(yintercept = 0) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 17),
      plot.title = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(size = 13)
      # , colour = clone_colour
      # strip.background = element_blank(),
      # strip.text = element_blank()
    ) +
    scale_x_continuous(breaks = xaxis_breaks, labels = as.character(xaxis_breaks)) +
    scale_colour_manual(values = c(A = "red", B = "blue")) +
    geom_segment(aes(x = startpos, xend = endpos, y = pred_CN, yend = pred_CN, colour = allele), size = 1.4) +
    geom_text_repel(aes(x = gene_start, y = max_pred_CN, label = gene_label), direction = "y", data = gene_info_clone[clone_fac == clone_levels[1]], size = 5) +
    labs(y = NULL, x = paste0("Chr ", chromosome," (Mb)"))
    # , title = "Predicted clone CN")
    
  g <- ggplot_gtable(ggplot_build(g))
  stripr <- which(grepl('strip-r', g$layout$name))

  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- clone_colour[k]
    k <- k + 1
  }

  pdf(file.path(tum_out_dir, paste0(tum_id, "_alpaca_", chromosome, "_seeding.pdf")), width = plot_width, height = plot_height)
  grid::grid.draw(g)
  # print(g)
  dev.off()
  # Also save in figures. directory
  pdf(paste0("figures/Fig2_", tum_id, "_alpaca_", chromosome, "_seeding.pdf"), width = plot_width, height = plot_height)
  grid::grid.draw(g)
  # print(g)
  dev.off()
}
