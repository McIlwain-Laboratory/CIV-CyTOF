# Sasha's code for CIV analysis
# Jun2024
# After you run Arcsinh_IH_with_sasha_mods.R script, use this script to produce cell abundance csv with p values

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

#for PB

# script below creates an output folder, and reads in the relevant preprocessed csv based on selected sample type
selected_samp_type <- "PB"
dir.create("output")
dir.create("output/cell_abundance")
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)

#-------------------- Task 1: Cell Abundance P value table --------------------

all_freq_data <- preprocessed_data %>% filter(reagent == "frequency")
all_freq_data <- all_freq_data %>% filter(stimulation == "U")

# d7 only

d7_freq_data <- all_freq_data %>% filter(time == "d7")
result_list <- list() # store group pair p vals
overall_result_list <- list() # store overall p vals
unique_populations <- unique(d7_freq_data$population)
unique_stimulations <- unique(d7_freq_data$stimulation) #should only be U

for (pop in unique_populations) {
  for (stim in unique_stimulations) {
    subset_data <- d7_freq_data %>% filter(population == pop, stimulation == stim)
    
    # ANOVA for all three groups
    overall_anova_result <- aov(feature ~ group, data = subset_data)
    overall_p_value <- summary(overall_anova_result)[[1]]$`Pr(>F)`[1]
    overall_df <- data.frame(cell_pop = pop, stim_cond = stim, anova_freq_pval_unadj = overall_p_value)
    overall_df$time = "d7"
    overall_result_list[[length(overall_result_list) + 1]] <- overall_df
    
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
    p_values_matrix_long$cell_pop <- pop
    p_values_matrix_long$stim_cond <- stim

    df <- p_values_matrix_long
    names(df)[2] <- "ttest_groupfreq_pval_unadj"
    
    
    df$time = "d7"
    result_list[[length(result_list) + 1]] <- df

  }
}

# Apply Bonferroni correction to overall p-values
d7_pairwise_p_values_df <- bind_rows(result_list)
d7_overall_p_values_df <- bind_rows(overall_result_list)
d7_overall_p_values_df$groups <- "Control Protein mRNA"

# Apply BH correction to all pairwise p-values
d7_pairwise_p_values_df$ttest_groupfreq_pval_BH <- p.adjust(d7_pairwise_p_values_df$ttest_groupfreq_pval_unadj, method = "BH")

#Apply BH correction to all of the overall ANOVA values
d7_overall_p_values_df$anova_freq_pval_BH <- p.adjust(d7_overall_p_values_df$anova_freq_pval_unadj, method = "BH")

# merge dataframes for the diff timepoints
ca_file_path <- file.path(getwd(), "output/cell_abundance")
file_name <- paste("freq_pairwise_p_vals_", selected_samp_type, ".csv", sep = "")
full_file_path <- file.path(ca_file_path, file_name)
all_pairwise_p_vals <- d7_pairwise_p_values_df
all_pairwise_p_vals$tissue_type = selected_samp_type
write.csv(all_pairwise_p_vals, full_file_path, row.names = FALSE)

file_name <- paste("freq_overall_p_vals_", selected_samp_type, ".csv", sep = "")
full_file_path <- file.path(ca_file_path, file_name)
all_overall_p_vals <- d7_overall_p_values_df
all_overall_p_vals$tissue_type = selected_samp_type
write.csv(all_overall_p_vals, full_file_path, row.names = FALSE)

#---------------------------------------------------------------------------

#same codechunk as above but for LN
# script below creates an output folder, and reads in the relevant preprocessed csv based on selected sample type
selected_samp_type <- "LN"
dir.create("output")
dir.create("output/cell_abundance")
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)

#-------------------- Task 1: Cell Abundance P value table (LN) --------------------

all_freq_data <- preprocessed_data %>% filter(reagent == "frequency")
all_freq_data <- all_freq_data %>% filter(stimulation == "U")

# d7 only

d7_freq_data <- all_freq_data %>% filter(time == "d7")
result_list <- list() # store group pair p vals
overall_result_list <- list() # store overall p vals
unique_populations <- unique(d7_freq_data$population)
unique_stimulations <- unique(d7_freq_data$stimulation)

for (pop in unique_populations) {
  for (stim in unique_stimulations) {
    subset_data <- d7_freq_data %>% filter(population == pop, stimulation == stim)
    
    # ANOVA for all three groups
    overall_anova_result <- aov(feature ~ group, data = subset_data)
    overall_p_value <- summary(overall_anova_result)[[1]]$`Pr(>F)`[1]
    overall_df <- data.frame(cell_pop = pop, stim_cond = stim, anova_freq_pval_unadj = overall_p_value)
    overall_df$time = "d7"
    overall_result_list[[length(overall_result_list) + 1]] <- overall_df
    
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
    p_values_matrix_long$cell_pop <- pop
    p_values_matrix_long$stim_cond <- stim
    
    df <- p_values_matrix_long
    names(df)[2] <- "ttest_groupfreq_pval_unadj"
    
    
    df$time = "d7"
    result_list[[length(result_list) + 1]] <- df
  }
}

# Apply Bonferroni correction to overall p-values
d7_pairwise_p_values_df <- bind_rows(result_list)
d7_overall_p_values_df <- bind_rows(overall_result_list)
d7_overall_p_values_df$groups <- "Control Protein mRNA"

# Apply BH correction to all pairwise p-values
d7_pairwise_p_values_df$ttest_groupfreq_pval_BH <- p.adjust(d7_pairwise_p_values_df$ttest_groupfreq_pval_unadj, method = "BH")

#Apply BH correction to all of the overall ANOVA values
d7_overall_p_values_df$anova_freq_pval_BH <- p.adjust(d7_overall_p_values_df$anova_freq_pval_unadj, method = "BH")

# merge dataframes for the diff timepoints
ca_file_path <- file.path(getwd(), "output/cell_abundance")
file_name <- paste("freq_pairwise_p_vals_", selected_samp_type, ".csv", sep = "")
full_file_path <- file.path(ca_file_path, file_name)
all_pairwise_p_vals <- d7_pairwise_p_values_df
all_pairwise_p_vals$tissue_type = selected_samp_type
write.csv(all_pairwise_p_vals, full_file_path, row.names = FALSE)

file_name <- paste("freq_overall_p_vals_", selected_samp_type, ".csv", sep = "")
full_file_path <- file.path(ca_file_path, file_name)
all_overall_p_vals <- d7_overall_p_values_df
all_overall_p_vals$tissue_type = selected_samp_type
write.csv(all_overall_p_vals, full_file_path, row.names = FALSE)

#------for Fig 2 supplement: "statistics for frequency comparisons"------
#code below combines LN and PB csv files into one
file_path_LN <- file.path(ca_file_path, "freq_pairwise_p_vals_LN.csv")
file_path_PB <- file.path(ca_file_path, "freq_pairwise_p_vals_PB.csv")
data_LN <- read.csv(file_path_LN)
data_PB <- read.csv(file_path_PB)
combined_data <- rbind(data_LN, data_PB)
combined_file_path <- file.path(ca_file_path, "freq_stats_pairwise_p_vals.csv")
write.csv(combined_data, combined_file_path, row.names = FALSE)

file_path_LN <- file.path(ca_file_path, "freq_overall_p_vals_LN.csv")
file_path_PB <- file.path(ca_file_path, "freq_overall_p_vals_PB.csv")
data_LN <- read.csv(file_path_LN)
data_PB <- read.csv(file_path_PB)
combined_data <- rbind(data_LN, data_PB)
combined_file_path <- file.path(ca_file_path, "freq_stats_overall_p_vals.csv")
write.csv(combined_data, combined_file_path, row.names = FALSE)

