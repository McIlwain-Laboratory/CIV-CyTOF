#Sasha's code for CIV analysis
#Apr2024
#AFTER RUNNING VOLCANO PLOTS SCRIPT, Run this for producing box plots for the top hits
library(ggplot2)
require(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(stats)
library("ggpubr")
library(ggrepel)

# This code is built to be mostly plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis will take place
setwd("Replace me")

stat_test_selection <- "pairwise" #choose either "overall" (overall ANOVA) or "pairwise" (pairwise comparison)

#------------------- Do not edit past this point (unless adjusting parameters)---------------------------

#=========== Preparing the dataframes ===========

#reads in the stats datafile depending on which stat test was selected above
current_directory <- getwd()
volc_data_filename<- paste0("mark_",stat_test_selection,"_pvals_both.csv")
volc_stat_data_path <- file.path(current_directory, "output", "signaling_markers", "volcano_plot",volc_data_filename)
volc_stat_data <- read.csv(volc_stat_data_path)

#add new col containing all meta information combined
volc_stat_data$pop_mark_tiss_stim <- paste(volc_stat_data$cell_pop,volc_stat_data$marker,volc_stat_data$tissue_type,volc_stat_data$stim_cond,sep="_")

#rename specific p_val column to general term
volc_stat_data <- volc_stat_data %>%
  { if("overall_p_val" %in% colnames(.)) rename(., p_val = overall_p_val) else . } %>%
  { if("pairwise_p_val" %in% colnames(.)) rename(., p_val = pairwise_p_val) else . }

#reads in preprocessed data which contains values for individual monkey ID's
preprocessed_volc_data_filename<- "volc_combined_preprocessed.csv"
volc_data_path <- file.path(current_directory, "preprocessed data for volc plot",preprocessed_volc_data_filename)
volc_data <- read.csv(volc_data_path)

#add new col containing all meta information combined
volc_data$pop_mark_tiss_stim <- paste(volc_data$population,volc_data$reagent,volc_data$tissue_type,volc_data$stimulation,sep="_")

unique_stims <- unique(volc_stat_data$stim_cond)
unique_tiss <- unique(volc_stat_data$tissue_type)

# initialize list to hold data frames for the 4 stim+tiss combos
data_frames <- list()

# Loop over combinations of stim and tissue types
for (s in unique_stims) {  #for R and P
  for (t in unique_tiss) {  #for LN and PB
    stat_subset <- volc_stat_data %>%  filter(stim_cond == s, tissue_type == t) %>%arrange(p_val) %>% slice_head(n = 10) # select the top 10 rows with the smallest 'p_val' values
    sig_ref_list <- unique(stat_subset$pop_mark_tiss_stim) #make reference list to identify significant pop_mark's
    
    for (q in sig_ref_list) {  # cycle through significant pop_marks
      sampleID_subset <- volc_data %>% filter(pop_mark_tiss_stim == q) #subset for only the top hits for this stim+tissue combo
      
      df_key <- paste0(t, "_", s) #dynamically construct dataframe name
      
      # append data to the corresponding data frame in the list
      if (!is.null(data_frames[[df_key]])) {
        data_frames[[df_key]] <- bind_rows(data_frames[[df_key]], sampleID_subset)
      } else {
        data_frames[[df_key]] <- sampleID_subset
      }
    }
  }
}

PB_P <- data_frames[["PB_P"]]
PB_R <- data_frames[["PB_R"]]
LN_P <- data_frames[["LN_P"]]
LN_R <- data_frames[["LN_R"]]

# Identify unique keys in both datasets
unique_keys_volc_data <- unique(volc_data$pop_mark_tiss_stim)
unique_keys_stat_subset <- unique(stat_subset$pop_mark_tiss_stim)

# Check if all keys in the subset are in the main data
not_in_stat_subset <- setdiff(unique_keys_volc_data, unique_keys_stat_subset)
not_in_volc_data <- setdiff(unique_keys_stat_subset, unique_keys_volc_data)

# Print out mismatches
print(not_in_stat_subset)
print(not_in_volc_data)
#=========== Making the boxplots ===========

#setting directory
current_directory <- getwd()
dir.create("output/signaling_markers/boxplots")

#----PB_P-----
dir.create("output/signaling_markers/boxplots/PB_P")
boxplot_folder_path <- file.path(current_directory, "output", "signaling_markers", "boxplots","PB_P")

pop_mark <- unique(PB_P$pop_mark_tiss_stim)

#for adding p_val info
p_vals_for_pmts <- volc_stat_data %>% 
  select(p_val, pop_mark_tiss_stim)

for (pm in pop_mark) {
  boxplotsubset <- PB_P %>% filter(pop_mark_tiss_stim == pm)
  
  # Retrieve the p_value for the current pop_mark_tiss_stim
  current_p_val <- filter(p_vals_for_pmts, pop_mark_tiss_stim == pm)$p_val
  
  # Create the ggplot object
  b <- ggplot(boxplotsubset, aes(x = group, y = feature, color = group)) +
    geom_boxplot() +
    geom_point(position = position_dodge(0.75), size = 2, alpha = 0.6) +  # Adjusted for correct positioning
    scale_color_brewer(palette = "Set1") +  # Apply a color palette for different groups
    labs(title = pm, x = "Vaccine Group", y = "Signaling Intensity (Normalized)") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.background = element_rect(fill = "white", colour = "grey80"),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    annotate("text", x = Inf, y = Inf, label = paste("p =", format(current_p_val, digits = 3)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black")  # Annotate with p-value
  
  # Generate filename and full path for saving
  boxplot_filename <- paste0("boxplot_", pm, ".png")
  full_boxplot_path <- file.path(boxplot_folder_path, boxplot_filename)
  
  # Save the plot with specified dimensions
  ggsave(full_boxplot_path, plot = b, width = 10, height = 8, bg = "white")  
}


#------PB_R-----
dir.create("output/signaling_markers/boxplots/PB_R")
boxplot_folder_path <- file.path(current_directory, "output", "signaling_markers", "boxplots","PB_R")

pop_mark <- unique(PB_R$pop_mark_tiss_stim)

#for adding p_val info
p_vals_for_pmts <- volc_stat_data %>% 
  select(p_val, pop_mark_tiss_stim)

for (pm in pop_mark) {
  boxplotsubset <- PB_R %>% filter(pop_mark_tiss_stim == pm)
  
  # Retrieve the p_value for the current pop_mark_tiss_stim
  current_p_val <- filter(p_vals_for_pmts, pop_mark_tiss_stim == pm)$p_val
  
  # Create the ggplot object
  b <- ggplot(boxplotsubset, aes(x = group, y = feature, color = group)) +
    geom_boxplot() +
    geom_point(position = position_dodge(0.75), size = 2, alpha = 0.6) +  # Adjusted for correct positioning
    scale_color_brewer(palette = "Set1") +  # Apply a color palette for different groups
    labs(title = pm, x = "Vaccine Group", y = "Signaling Intensity (Normalized)") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.background = element_rect(fill = "white", colour = "grey80"),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    annotate("text", x = Inf, y = Inf, label = paste("p =", format(current_p_val, digits = 3)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black")  # Annotate with p-value
  
  # Generate filename and full path for saving
  boxplot_filename <- paste0("boxplot_", pm, ".png")
  full_boxplot_path <- file.path(boxplot_folder_path, boxplot_filename)
  
  # Save the plot with specified dimensions
  ggsave(full_boxplot_path, plot = b, width = 10, height = 8, bg = "white")  
}



#-----LN_P-----
dir.create("output/signaling_markers/boxplots/LN_P")
boxplot_folder_path <- file.path(current_directory, "output", "signaling_markers", "boxplots","LN_P")

pop_mark <- unique(LN_P$pop_mark_tiss_stim)

#for adding p_val info
p_vals_for_pmts <- volc_stat_data %>% 
  select(p_val, pop_mark_tiss_stim)

for (pm in pop_mark) {
  boxplotsubset <- LN_P %>% filter(pop_mark_tiss_stim == pm)
  
  # Retrieve the p_value for the current pop_mark_tiss_stim
  current_p_val <- filter(p_vals_for_pmts, pop_mark_tiss_stim == pm)$p_val
  
  # Create the ggplot object
  b <- ggplot(boxplotsubset, aes(x = group, y = feature, color = group)) +
    geom_boxplot() +
    geom_point(position = position_dodge(0.75), size = 2, alpha = 0.6) +  # Adjusted for correct positioning
    scale_color_brewer(palette = "Set1") +  # Apply a color palette for different groups
    labs(title = pm, x = "Vaccine Group", y = "Signaling Intensity (Normalized)") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.background = element_rect(fill = "white", colour = "grey80"),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    annotate("text", x = Inf, y = Inf, label = paste("p =", format(current_p_val, digits = 3)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black")  # Annotate with p-value
  
  # Generate filename and full path for saving
  boxplot_filename <- paste0("boxplot_", pm, ".png")
  full_boxplot_path <- file.path(boxplot_folder_path, boxplot_filename)
  
  # Save the plot with specified dimensions
  ggsave(full_boxplot_path, plot = b, width = 10, height = 8, bg = "white")  
}


#LN_R
dir.create("output/signaling_markers/boxplots/LN_R")
boxplot_folder_path <- file.path(current_directory, "output", "signaling_markers", "boxplots","LN_R")

pop_mark <- unique(LN_R$pop_mark_tiss_stim)

#for adding p_val info
p_vals_for_pmts <- volc_stat_data %>% 
  select(p_val, pop_mark_tiss_stim)

for (pm in pop_mark) {
  boxplotsubset <- LN_R %>% filter(pop_mark_tiss_stim == pm)
  
  # Retrieve the p_value for the current pop_mark_tiss_stim
  current_p_val <- filter(p_vals_for_pmts, pop_mark_tiss_stim == pm)$p_val
  
  # Create the ggplot object
  b <- ggplot(boxplotsubset, aes(x = group, y = feature, color = group)) +
    geom_boxplot() +
    geom_point(position = position_dodge(0.75), size = 2, alpha = 0.6) +  # Adjusted for correct positioning
    scale_color_brewer(palette = "Set1") +  # Apply a color palette for different groups
    labs(title = pm, x = "Vaccine Group", y = "Signaling Intensity (Normalized)") +
    theme_minimal(base_size = 14) +  
    theme(
      panel.background = element_rect(fill = "white", colour = "grey80"),
      plot.background = element_rect(fill = "white", colour = NA)
    ) +
    annotate("text", x = Inf, y = Inf, label = paste("p =", format(current_p_val, digits = 3)), 
             hjust = 1.1, vjust = 2, size = 5, color = "black")  # Annotate with p-value
  
  # Generate filename and full path for saving
  boxplot_filename <- paste0("boxplot_", pm, ".png")
  full_boxplot_path <- file.path(boxplot_folder_path, boxplot_filename)
  
  # Save the plot with specified dimensions
  ggsave(full_boxplot_path, plot = b, width = 10, height = 8, bg = "white")  
}




