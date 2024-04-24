# Functions for running survival analysis in R
### Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(corrplot))

#######################
### Preprocess data ###
#######################

process_edgeevents_for_survival_analysis <- function(distances_df, evo_df) {
  dist_tmp <- copy(distances_df)
  evo_df[, sc_GD := (First_GD == "Subclonal" | Second_GD == "Subclonal")]
  scna_ith <- evo_df[tumour_id %in% unique(distances_df$tumour_id), .(tumour_id, frac_abberant_genom_subcl, First_GD, Second_GD, sc_GD, genom_frac_event, purity_bootstrapped, ploidy_bootstrapped)]

  dist_tmp[, n_clones := length(unique(clone)), by = tumour_id]
  # Merge with scna_ith metrics from tx421:
  dist_tmp <- merge.data.table(dist_tmp, scna_ith, by = "tumour_id", all.x = TRUE, all.y = FALSE)

  ### Compute the mean & max pairwise distance between any clone and the MRCA
  dist_tmp[parent_clone != "diploid", sum_subclonal_gains := sum(total_interval_gains_bin), by = tumour_id]
  dist_tmp[parent_clone != "diploid", sum_subclonal_losses := sum(total_interval_losses_bin), by = tumour_id]
  dist_tmp[parent_clone != "diploid", slope_variance := var(slope_raw, na.rm = TRUE), by = tumour_id]
  dist_tmp[parent_clone != "diploid" & is.na(slope_variance), slope_variance := 0]
  dist_tmp[, total_tum_gains := sum(total_interval_gains_bin), by = tumour_id]
  dist_tmp[, total_tum_losses := sum(total_interval_losses_bin), by = tumour_id]
  dist_tmp[, frac_subclonal_scnas := sum_subclonal_cn_events_tumour / sum_cn_events_tumour]
  dist_tmp[, frac_subclonal_gains := sum_subclonal_gains / total_tum_gains]
  dist_tmp[, frac_subclonal_losses := sum_subclonal_losses / total_tum_losses]
  dist_tmp[, n_scnas_per_clone := sum_cn_events_tumour / n_clones]
  dist_tmp[, n_amplitude_scnas_per_clone := sum_cn_amplitude_events_tumour / n_clones]
  dist_tmp[, n_amplitude_gains_per_clone := sum_amplitude_gains_tumour / n_clones]
  dist_tmp[, n_amplitude_losses_per_clone := sum_amplitude_losses_tumour / n_clones]
  return(dist_tmp)
}

process_pairwisecndist_for_survival_analysis <- function(pairwise_dist, ln_prims) {
  ### Compute the mean & max pairwise distance between any two clones in the tumour
  # Remove the trunk
  pairwise_dist_sc_tmp <- pairwise_dist[(clone1 != "diploid") & (clone2 != "diploid")]
  # Remove clones that are exclusively in the LN:
  pairwise_dist_sc_tmp[, ID_1 := paste0(tumour_id, "_", clone1)]
  pairwise_dist_sc_tmp[, ID_2 := paste0(tumour_id, "_", clone2)]
  pairwise_dist_sc_tmp <- merge.data.table(pairwise_dist_sc_tmp, ln_prims[, .(ID_1 = ID, exclusive_LN_1 = exclusive_LN)], by = "ID_1", all.x = TRUE, all.y = FALSE)
  pairwise_dist_sc_tmp <- merge.data.table(pairwise_dist_sc_tmp, ln_prims[, .(ID_2 = ID, exclusive_LN_2 = exclusive_LN)], by = "ID_2", all.x = TRUE, all.y = FALSE)
  # pairwise_dist_sc <- copy(pairwise_dist_sc_tmp)
  pairwise_dist_sc <- pairwise_dist_sc_tmp[is.na(exclusive_LN_1) & is.na(exclusive_LN_2)]


  ### Get allele-specific distances:
  pairwise_dist_sc[, maxtum_l1_dist_pairwise_allelespecific := max(c(l1_cn_distance_A, l1_cn_distance_B)), by = .(tumour_id, seg_size)]
  pairwise_dist_sc[, maxtum_l1_w_dist_pairwise_allelespecific := max(c(l1_cn_distance_A_w, l1_cn_distance_B_w)), by = .(tumour_id, seg_size)]
  pairwise_dist_sc[, maxtum_l2_dist_pairwise_allelespecific := max(c(l2_cn_distance_A, l2_cn_distance_B)), by = .(tumour_id, seg_size)]
  pairwise_dist_sc[, maxtum_l2_w_dist_pairwise_allelespecific := max(c(l2_cn_distance_A_w, l2_cn_distance_B_w)), by = .(tumour_id, seg_size)]

  pairwise_dist_sc[, l2_dist_pairwise_sumallele := l2_cn_distance_A + l2_cn_distance_B]
  pairwise_dist_sc[, maxtum_l2_dist_pairwise_sumallele := max(l2_dist_pairwise_sumallele), by = .(tumour_id, seg_size)]
  pairwise_dist_sc[, maxtum_l2_dist_pairwise_AB := max(l2_cn_distance_AB), by = .(tumour_id, seg_size)]
  # pairwise_dist_sc[, pair_is_max := (maxtum_l2_dist_pairwise_allelespecific == l2_cn_distance_A) | (maxtum_l2_dist_pairwise_allelespecific == l2_cn_distance_B)] 
  # maxpairs <- pairwise_dist_sc[pair_is_max == TRUE & seg_size == "all"]
  # length(maxpairs[!is.na(exclusive_LN_1) | !is.na(exclusive_LN_2), unique(tumour_id)])
  # length(pairwise_dist_sc[!is.na(exclusive_LN_1) | !is.na(exclusive_LN_2), unique(tumour_id)])

  return(pairwise_dist_sc)
}

plot_cor_matrix <- function(tum_df_raw, out_path) {
	data_df <- as.data.frame(tum_df_raw)
	rownames(data_df) <- tum_df_raw$tumour_id
	data_df$tumour_id <- NULL
	res <- cor(data_df)

	pdf(out_path, width = 8, height = 8)
	corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
	dev.off()
}

preprocess_survival_data <- function(all_patient_df, tum_df, metrics_list) {
	# Process clinical data in same way as in TRACERx421
	clinical_data_df <- all_patient_df  %>%
    dplyr::select(cruk_id, tumour_id_per_patient, tx100, age, pathologyTNM, packyears, smoking_status_merged, 
                  histology_lesion1, adjuvant_treatment_YN, sex,
                  cens_os, os_time, cens_dfs, dfs_time, cens_dfs_any_event, dfs_time_any_event, cens_lung_event, lung_event_time, 
                  Recurrence_time_use, Relapse_cat_new
    ) %>%
    filter(!cruk_id %in% c("CRUK0682")) %>%  # no lesion 1 sampled : CRUK0682
    mutate(tumour_id_per_patient = case_when(cruk_id == "CRUK0704" ~  "CRUK0704_Tumour1.2.3",  # these are collision tumours
                                             cruk_id == "CRUK0881" ~  "CRUK0881_Tumour1.2",
                                             TRUE ~ tumour_id_per_patient)) %>%
    mutate(cens_dfs = case_when(cruk_id %in% c("CRUK0512", "CRUK0373","CRUK0428","CRUK0511") ~ 0, 
                                TRUE ~ as.numeric(cens_dfs)),
           dfs_time = case_when(cruk_id %in% c("CRUK0512", "CRUK0373","CRUK0428","CRUK0511") ~ dfs_time_any_event, 
                                TRUE ~ dfs_time )) %>%
    
		mutate(Relapse_cat_new = case_when(cruk_id %in% c("CRUK0173","CRUK0249", "CRUK0781","CRUK0129") ~ "Unknown Site", # "CRUK0173" : relapsed before death but no CT reports available. "CRUK0249", "CRUK0781","CRUK0129": lung cancer death as first event 
                                       TRUE ~ as.character(Relapse_cat_new) )) %>%
    mutate(Relapse_cat_new = factor(Relapse_cat_new, levels = c("No rec", "Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site"))) %>%
    mutate(is.rec = case_when(Relapse_cat_new %in% c("Intrathoracic", "Intra & Extra", "Extrathoracic", "Unknown Site") ~ TRUE, 
                              TRUE ~ FALSE)) %>%
    mutate(time.to.rec = case_when(is.rec == TRUE ~ Recurrence_time_use)) %>%
    mutate(Relapse.Site = case_when(Relapse_cat_new %in% c("Intrathoracic") ~ "Intrathoracic",
                                    Relapse_cat_new %in% c("Intra & Extra", "Extrathoracic") ~ "Extrathoracic")) %>%
    mutate(Relapse.Site = factor(Relapse.Site , levels = c("Intrathoracic", "Extrathoracic"))) %>%
    
    mutate(TNM_stage = case_when(pathologyTNM %in% c("IA","IB") ~ "I",
                                 pathologyTNM %in% c("IIA","IIB") ~ "II",
                                 pathologyTNM %in% c("IIIA","IIIB") ~ "III")) %>%
    mutate(pTNM_stage = case_when(pathologyTNM %in% c("IIIA","IIIB") ~ "III",
                                  TRUE ~ pathologyTNM)) %>%
    mutate(adjuvant_tx = factor(adjuvant_treatment_YN, levels = c("No adjuvant", "Adjuvant"))) %>%
    mutate(packyears = packyears,
           smoking_history = factor(smoking_status_merged, levels = c("Never Smoked","Ex-Smoker",  "Smoker" ))) %>%
    mutate(histology = case_when(histology_lesion1 == "Invasive adenocarcinoma" ~ "LUAD",
                                 histology_lesion1 == "Squamous cell carcinoma" ~ "LUSC",
                                 TRUE ~ "Other")) %>%
    mutate(histology_LUAD = case_when(histology_lesion1 == "Invasive adenocarcinoma" ~ "LUAD",
                                      TRUE ~ "Non-LUAD")) %>%
    mutate(histology = factor(histology, levels = c("LUAD", "LUSC", "Other")),
           histology_LUAD = factor(histology_LUAD, levels = c("LUAD","Non-LUAD")),
           adjuvant_tx = factor(adjuvant_treatment_YN, levels = c("No adjuvant", "Adjuvant")),
           smoking_history = factor(smoking_status_merged, levels = c("Never Smoked","Ex-Smoker",  "Smoker" )))  %>%
    dplyr::select(-pathologyTNM, -histology_lesion1, -adjuvant_treatment_YN, -smoking_status_merged )
  

	### Deal with multi-tumour cases: take the tumour with the highest stage (column tumour_id_per_patient from all_patient_df):
	## SCNA-ITH  
	cols_to_keep <- c("tumour_id", metrics_list)
	metrics_df <- tum_df[, ..cols_to_keep]
	metrics_df[, cruk_id := gsub("-Cluster.*", "", tumour_id)]
	metrics_df[, tumour_id_per_patient := gsub("-Cluster", "_Tumour", tumour_id)]
	metrics_df[, tumour_id := NULL]

	# Deal with collision tumours:  CRUK0704, CRUK0881 (CRUK0039 also had collision tumour, though only 1 of them was analysed in primary cohort )
	# use max SCNA-ITH (frac_abberant_genom_subcl) to represent the patient's SCNA-ITH, based on NEJM2017 that high SCNA-ITH is associated with poor prognosis
	# Collision tumours
	evo_metrics_collision <- metrics_df[cruk_id %in% c("CRUK0704", "CRUK0881")]
	for (m in metrics_list) {
		evo_metrics_collision[, eval(m) := max(get(m), na.rm = TRUE), by = cruk_id]
	}
	evo_metrics_collision <- evo_metrics_collision %>%
    mutate(tumour_id_per_patient = case_when(grepl("CRUK0704", tumour_id_per_patient) ~ "CRUK0704_Tumour1.2.3", 
                                 grepl("CRUK0881", tumour_id_per_patient) ~ "CRUK0881_Tumour1.2", 
                                 TRUE ~ tumour_id_per_patient))
	evo_metrics_collision <- unique(evo_metrics_collision)
	
	# Non-collision tumours
	evo_metrics_noncollision <- metrics_df[!cruk_id %in% c("CRUK0704", "CRUK0881")]
	metrics_df <- rbindlist(list(evo_metrics_noncollision, evo_metrics_collision), use.names = TRUE)

	# This will be the analytical cohort  n=387 (CRUK0682 excluded from survival analysis, and some tumours didnt get through ALPACA)
	surv_df_surv_pre <- merge.data.table(clinical_data_df, metrics_df, by = c("cruk_id", "tumour_id_per_patient"), all = FALSE)
	surv_df_surv_pre <- unique(surv_df_surv_pre)
	## Note, CRUK0721_Tumour1 is a single region tumour but is the highest assigned stage tumour for that case
	## So this patient gets removed


	### Get median/tertiles cutoff per analytical cohort for each metric
	for (cur_metric in metrics_list) {
		metric_distribution <- surv_df_surv_pre[, as.numeric(get(cur_metric))]
		
		# Compute median / tertiles
		median_metric <- median(metric_distribution)
  		tertiles_metric <- quantile(metric_distribution, probs = c(1/3, 2/3))
  
		# Categorise each variable
		cur_group_median <- paste0(cur_metric, "_cat")
		surv_df_surv_pre[, eval(cur_group_median) := ifelse(get(cur_metric) >= median_metric, "high", "low")]

		cur_group_tert <- paste0(cur_metric, "_tertiles")
		surv_df_surv_pre[, eval(cur_group_tert) := ifelse(get(cur_metric) >= as.numeric(tertiles_metric[2]), "high",
																								ifelse((get(cur_metric) >= as.numeric(tertiles_metric[1])) & (get(cur_metric) < as.numeric(tertiles_metric[2])), "mid",
																									ifelse(get(cur_metric) < as.numeric(tertiles_metric[1]), "low", NA)))]
		# cat('\n\n\n')
		# print(cur_metric)
		# print(sd(metric_distribution))

		# Save the standard deviation of each metric
		cur_group_sd <- paste0(cur_metric, "_sd")
		surv_df_surv_pre[, eval(cur_group_sd) := sd(metric_distribution)]

		# Rescale continuous variables by dividing by S.D.
		# to get the increase of their Hazard ratio per one s.d. 
		metric_by_sd <- paste0(cur_metric, "_by_sd")
		surv_df_surv_pre[, eval(metric_by_sd) := get(cur_metric) / get(cur_group_sd)]
	}
	
	# Rescale clinical variables age and packyears to get 
	# the increase of their Hazard ratio per every 10 years
	in_data_surv <- surv_df_surv_pre %>%
		mutate(`Age (+10 years)` = age/10,  # to get OR per 10 increase
				`PackYears (+10)` = packyears/10)
	return(in_data_surv)
}

#########################
### survival analysis ###  uni-Cox  & KM plot  for 2 groups
#########################

plot_km_bymedian <- function(in_data_surv, group_vec, out_pdf, title_start = 'DFS', time_type = "dfs_time", status_type = "cens_dfs") {
	
	surv_df <- as.data.frame(in_data_surv)
	for (n in 1:length(group_vec)) {
		out_pdf_group <- gsub('.pdf', paste0(group_vec[n], '.pdf'), out_pdf)
		
		print(group_vec[n])
		xlim_value <- 2000
		xbreak_value <- 500
		
		surv_df$group <- surv_df[[group_vec[n]]]
		surv_df$time <- surv_df[[time_type]]
		surv_df$status <- surv_df[[status_type]]  
		plot_title <- gsub("_cat|_tertiles", "", group_vec[n])

		surv_df <- surv_df %>%
			dplyr::filter(!is.na(group), !is.na(status)) %>%
			mutate(group = factor(group, levels = c("low", "high")))
	
		
		fit <- survfit(Surv(time, status) ~ group, data = surv_df)
		cox_fit <- coxph( Surv(time, status) ~ group, data = surv_df )
		CI_tab <- summary(cox_fit)$conf.int
		round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)

		label <- paste0('Hazard ratio:\n',
										'High vs Low - ', round_cust(CI_tab[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab[, 3][1]), '-', round_cust(CI_tab[, 4][1]),
										# '\nPolycl. & Polyphyl. vs Monocl. - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
										'\n (P=', signif(summary(cox_fit)$coefficients[,5][1],2),')')
		labs_vec <-  c('Low', 'High')
	
		
		colours <- brewer.pal(3, 'Set1')[1:2]
		
	
		g <- ggsurvplot(
			data = surv_df,
			fit = fit,
			pval = FALSE, #TRUE,
			pval.method = FALSE, #TRUE,
			#test.for.trend = TRUE, 
			risk.table = TRUE,
			#risk.table.height = 0.16,
			legend.title = '',
			font.legend = 16,
			palette = rev(colours),
			xlab = "Days", 
			ylab = title_start,
			xlim = c(0,xlim_value), 
			font.tickslab = c(18),
			font.x = c(20),
			font.y = c(20),
			break.x.by = xbreak_value,
			legend.labs = labs_vec,
			title = plot_title) 
		
		g$plot <- g$plot +
			annotate("text", hjust = 0,
								x = 10, y = 0.15, # x and y coordinates of the text
								label = label, size = 5) +
			theme(plot.title = element_text(size = 18))
		g$table <- g$table + 
			labs(title = '', y = '', x = '') 
		
		pdf(out_pdf_group, width = 5, height = 7, onefile=FALSE)
		print(g)
		dev.off()
	}

}

plot_km_bytertiles <- function(in_data_surv, group_vec, out_pdf, title_start = 'DFS', time_type = "dfs_time", status_type = "cens_dfs") {
	
	surv_df <- as.data.frame(in_data_surv)
	
	for (n in 1:length(group_vec)) {
		out_pdf_group <- gsub('.pdf', paste0(group_vec[n], '.pdf'), out_pdf)

		print(group_vec[n])
		xlim_value <- 2000
		xbreak_value <- 500
		
		surv_df$group <- surv_df[[group_vec[n]]]
		surv_df$time <- surv_df[[time_type]]
		surv_df$status <- surv_df[[status_type]]  		
		plot_title <- gsub("_cat|_tertiles", "", group_vec[n])
		
		surv_df <- surv_df %>%
        dplyr::filter(!is.na(group), !is.na(status)) %>%
        mutate(group = factor(group, levels = c("low", "mid", "high")))%>%
        mutate(group_2 = case_when(group == "high" ~ "high",
                                   TRUE ~ "all other")) %>%
        mutate(group_2 = factor(group_2, levels = c("all other", "high"))) %>%
				mutate(group_3 = case_when(group == "low" ~ "low",
                                   TRUE ~ "all other")) %>%
        mutate(group_3 = factor(group_3, levels = c("all other", "low")))
		
		fit <- survfit(Surv(time, status) ~ group, data = surv_df)
		cox_fit <- coxph( Surv(time, status) ~ group, data = surv_df )
		CI_tab <- summary(cox_fit)$conf.int
		round_cust <- function(x, digits = 2) format(round(x, digits), nsmall = digits)
      
		fit <- survfit(Surv(time, status) ~ group, data = surv_df)
		cox_fit <- coxph( Surv(time, status) ~ group, 
											data = surv_df )
		CI_tab <- summary(cox_fit)$conf.int
		
		cox_fit2 <- coxph( Surv(time, status) ~ group_2, 
												data = surv_df )
		CI_tab2 <- summary(cox_fit2)$conf.int
		
		cox_fit3 <- coxph( Surv(time, status) ~ group_3, 
												data = surv_df )
		CI_tab3 <- summary(cox_fit3)$conf.int

		label <- paste0('Hazard ratio:\n',
										'High vs Low - ', round_cust(CI_tab[, 1][2]), ' ; 95% CI: ',round_cust(CI_tab[, 3][2]), '-', round_cust(CI_tab[, 4][2]),
										' (P=', signif(summary(cox_fit)$coefficients[,5][2],2),')',
										'\nHigh vs all other - ', round_cust(CI_tab2[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab2[, 3][1]), '-', round_cust(CI_tab2[, 4][1]),
										' (P=', signif(summary(cox_fit2)$coefficients[,5][1],2),')',
										'\nLow vs all other - ', round_cust(CI_tab3[, 1][1]), ' ; 95% CI: ',round_cust(CI_tab3[, 3][1]), '-', round_cust(CI_tab3[, 4][1]),
										' (P=', signif(summary(cox_fit3)$coefficients[,5][1],2),')'
		)
		labs_vec <-  c('Low', 'Mid', 'High')
		
		colours <- brewer.pal(3, 'Set1')[c(1,3,2)]
		
		
		g <- ggsurvplot(
			data = surv_df,
			fit = fit,
			pval = FALSE, #TRUE,
			pval.method = FALSE, #TRUE,
			#test.for.trend = TRUE, 
			risk.table = TRUE,
			#risk.table.height = 0.16,
			legend.title = '',
			font.legend = 16,
			palette = rev(colours),
			xlab = "Days", 
			ylab = title_start,
			xlim = c(0, xlim_value), 
			font.tickslab = c(18),
			font.x = c(20),
			font.y = c(20),
			break.x.by = xbreak_value,
			legend.labs = labs_vec,
			title = plot_title)
		
		g$plot <- g$plot +
			annotate("text", hjust = 0,
								x = 10, y = 0.15, # x and y coordinates of the text
								label = label, size = 4) +
			theme(plot.title = element_text(size = 18))
		g$table <- g$table + 
			labs(title = '', y = '', x = '')  #theme(plot.title = element_blank())

		pdf(out_pdf_group, width = 5, height = 7, onefile=FALSE)
		print(g)
		dev.off()
	}
}


run_multivariate_analysis <- function(in_data_surv,
										variables_cox,
										out_pdf,
										title_start = 'DFS',
										time_type = "dfs_time",
										status_type = "cens_dfs") {

	surv_df <- as.data.frame(in_data_surv)

	# as numerical
	for (var in c(variables_cox)) {
		surv_df <- surv_df[!is.na(surv_df[[var]]),]
	}

	surv_df$time <- surv_df[[time_type]]
	surv_df$status <- surv_df[[status_type]]  

	f <- as.formula(paste("Surv(time, status)", paste(paste0("`", variables_cox, "`"), collapse = " + "), sep = " ~ "))
	cox_out <- coxph(f, data = surv_df)

	pdf(out_pdf, width = 14, height = 13, onefile=FALSE)
	g <- survminer::ggforest(model = cox_out,
									data = surv_df,
									main = paste0("Hazard ratio: ", title_start),
									fontsize = 1.8,
									refLabel = "reference") + 
			theme(panel.grid.major = element_line(colour = "black"))
	print(g)
	dev.off()
}
