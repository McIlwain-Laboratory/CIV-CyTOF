#Sasha's code for CIV analysis
#Oct2023
#After you run Arcsinh_IH_with_sasha_mods.R script, use this script to produce signaling markers csv with p values
library(ggplot2)
require(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(stats)

# This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis is taking place
setwd("Replace me")


#------------------- Do not edit past this point ---------------------------
selected_samp_type <- "PB"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/signaling_markers")
mark_file_path <- file.path(getwd(), "output/signaling_markers")

#-------------------- Task 3 and 4: Cell Signaling Levels p_val table --------------------

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#make empty lists to build into
unique_stims <- unique(all_mark_data$stimulation)
unique_cell_pops <- unique(all_mark_data$population)
unique_markers <- unique(all_mark_data$reagent)
resulting_list <- list()
overall_result_list <- list() # store overall anova p vals
#make sure to bind above below to create finalized list

#nested loop for making plots showing all populations x Days x Treatments
for (s in unique_stims) {
  for (c in unique_cell_pops) {
    for (m in unique_markers) {
        #create subset of data
    subset_data <- all_mark_data %>% filter(stimulation == s, population == c, reagent == m)
    
    #run anova
    anova_result <- aov(feature ~ group, data = subset_data)
    
    #pull p-value from anova
    overall_p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
    overall_df <- data.frame(cell_pop = c, stim_cond = s, marker = m, anova_marker_pval_unadj = overall_p_value)
    overall_result_list[[length(overall_result_list) + 1]] <- overall_df
    
    #do post-hoc ttest
    t_test <- pairwise.t.test(subset_data$feature, subset_data$group, p.adjust.method = "none")
    p_values_matrix <- as.data.frame(t_test$p.value)
    p_values_matrix_long <- p_values_matrix %>%
      rownames_to_column(var = "RowName") %>%
      pivot_longer(
        cols = -RowName, 
        names_to = "ColumnName", 
        values_to = "p_value") %>%
      unite("group_pair", RowName, ColumnName, sep = "-") %>%
      filter(!is.na(p_value))
    p_values_matrix_long$cell_pop <- c
    p_values_matrix_long$stim_cond <- s
    p_values_matrix_long$marker <- m
    
    #combine ttest information into dataframe and add to list
    df <- p_values_matrix_long
    names(df)[2] <- "ttest_groupmarker_pval_unadj"
    resulting_list[[length(resulting_list) + 1]] <- df
      }}}

marker_ttest_p_values_df <- bind_rows(resulting_list)
marker_anova_p_values_df <- bind_rows(overall_result_list)

#Need to correct for multiple hypothesis with BH p.adj(dataframe$pvalue, method = "BH")
marker_ttest_p_values_df$ttest_groupmarker_pval_BH <- p.adjust(marker_ttest_p_values_df$ttest_groupmarker_pval_unadj, method = "BH")
marker_anova_p_values_df$anova_marker_pval_BH <- p.adjust(marker_anova_p_values_df$anova_marker_pval_unadj, method = "BH")

#general file path
mark_file_path <- file.path(getwd(), "output/signaling_markers")

#file path for anova values
anova_file_name <- paste("sig_markers_anova_p_vals_", selected_samp_type, ".csv",sep="")
full_anova_file_path <- file.path(mark_file_path, anova_file_name)
marker_anova_p_values_df$tissue_type <- selected_samp_type
write.csv(marker_anova_p_values_df, full_anova_file_path, row.names = FALSE)

#file path for ttest values
ttest_file_name <- paste("sig_markers_ttest_p_vals_", selected_samp_type, ".csv",sep="")
full_ttest_file_path <- file.path(mark_file_path, ttest_file_name)
marker_ttest_p_values_df$tissue_type <- selected_samp_type
write.csv(marker_ttest_p_values_df, full_ttest_file_path, row.names = FALSE)

######--------------------------------------------------------------------------

# code chunk below is same as above but for LN instead of PB
selected_samp_type <- "LN"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/signaling_markers")
mark_file_path <- file.path(getwd(), "output/signaling_markers")

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#make empty lists to build into
unique_stims <- unique(all_mark_data$stimulation)
unique_cell_pops <- unique(all_mark_data$population)
unique_markers <- unique(all_mark_data$reagent)
resulting_list <- list()
overall_result_list <- list() # store overall anova p vals
#make sure to bind above below to create finalized list

#nested loop for making plots showing all populations x Days x Treatments
for (s in unique_stims) {
  for (c in unique_cell_pops) {
    for (m in unique_markers) {
        #create subset of data
        subset_data <- all_mark_data %>% filter(stimulation == s, population == c, reagent == m)
        
        #run anova
        anova_result <- aov(feature ~ group, data = subset_data)
        
        #pull p-value from anova
        overall_p_value <- summary(anova_result)[[1]]$`Pr(>F)`[1]
        overall_df <- data.frame(cell_pop = c, stim_cond = s, marker = m, anova_marker_pval_unadj = overall_p_value)
        overall_result_list[[length(overall_result_list) + 1]] <- overall_df
        
        #do post-hoc ttest
        t_test <- pairwise.t.test(subset_data$feature, subset_data$group, p.adjust.method = "none")
        p_values_matrix <- as.data.frame(t_test$p.value)
        p_values_matrix_long <- p_values_matrix %>%
          rownames_to_column(var = "RowName") %>%
          pivot_longer(
            cols = -RowName, 
            names_to = "ColumnName", 
            values_to = "p_value") %>%
          unite("group_pair", RowName, ColumnName, sep = "-") %>%
          filter(!is.na(p_value))
        p_values_matrix_long$cell_pop <- c
        p_values_matrix_long$stim_cond <- s
        p_values_matrix_long$marker <- m
        
        #combine ttest information into dataframe and add to list
        df <- p_values_matrix_long
        names(df)[2] <- "ttest_groupmarker_pval_unadj"
        resulting_list[[length(resulting_list) + 1]] <- df
      }}}

marker_ttest_p_values_df <- bind_rows(resulting_list)
marker_anova_p_values_df <- bind_rows(overall_result_list)

#Need to correct for multiple hypothesis with BH p.adj(dataframe$pvalue, method = "BH")
marker_ttest_p_values_df$ttest_groupmarker_pval_BH <- p.adjust(marker_ttest_p_values_df$ttest_groupmarker_pval_unadj, method = "BH")
marker_anova_p_values_df$anova_marker_pval_BH <- p.adjust(marker_anova_p_values_df$anova_marker_pval_unadj, method = "BH")

#general file path
mark_file_path <- file.path(getwd(), "output/signaling_markers")

#file path for anova values
anova_file_name <- paste("sig_markers_anova_p_vals_", selected_samp_type, ".csv",sep="")
full_anova_file_path <- file.path(mark_file_path, anova_file_name)
marker_anova_p_values_df$tissue_type <- selected_samp_type
write.csv(marker_anova_p_values_df, full_anova_file_path, row.names = FALSE)

#file path for ttest values
ttest_file_name <- paste("sig_markers_ttest_p_vals_", selected_samp_type, ".csv",sep="")
full_ttest_file_path <- file.path(mark_file_path, ttest_file_name)
marker_ttest_p_values_df$tissue_type <- selected_samp_type
write.csv(marker_ttest_p_values_df, full_ttest_file_path, row.names = FALSE)

#-------------------------------------------------------------------------------
#Combine the files into a single file per statistical test

#general file path
mark_file_path <- file.path(getwd(), "output/signaling_markers")

# combine both tissues into one file for ANOVA
file_path_anova_LN_mark <- file.path(mark_file_path, "sig_markers_anova_p_vals_LN.csv")
file_path_anova_PB_mark <- file.path(mark_file_path, "sig_markers_anova_p_vals_PB.csv")
sigdata_anova_LN <- read.csv(file_path_anova_LN_mark)
sigdata_anova_PB <- read.csv(file_path_anova_PB_mark)
sigcombined_anova_data <- rbind(sigdata_anova_LN, sigdata_anova_PB)
sigcombined_anova_file_path <- file.path(mark_file_path, "sig_markers_anova_p_vals_both.csv")
write.csv(sigcombined_anova_data, sigcombined_anova_file_path, row.names = FALSE)

# combine both tissues into one file for ttest
file_path_ttest_LN_mark <- file.path(mark_file_path, "sig_markers_ttest_p_vals_LN.csv")
file_path_ttest_PB_mark <- file.path(mark_file_path, "sig_markers_ttest_p_vals_PB.csv")
sigdata_ttest_LN <- read.csv(file_path_ttest_LN_mark)
sigdata_ttest_PB <- read.csv(file_path_ttest_PB_mark)
sigcombined_ttest_data <- rbind(sigdata_ttest_LN, sigdata_ttest_PB)
sigcombined_ttest_file_path <- file.path(mark_file_path, "sig_markers_ttest_p_vals_both.csv")
write.csv(sigcombined_ttest_data, sigcombined_ttest_file_path, row.names = FALSE)