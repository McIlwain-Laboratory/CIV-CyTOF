#Sasha's code for CIV analysis
#Oct2024
#Run script to produce correlation plots for top freq differences vs metadata 
#Fig 5a in manuscript
library(ggplot2)
require(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(stats)
library(readr)

# This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis is taking place
setwd("Replace me")

# download metadata to that folder and paste the file path to it
metadata<- read.csv("Replace me/CIV01 Lung Metadata RNA and Pathology.csv") 

# paste file path to frequency stats p value table produced from script CIV_abund_p_vals_V3.R
pvals <- read.csv("Replace me")

#choose tissue to perform analysis: type either "PB" or "LN"
selected_samp_type <- "LN"

#------------------- Do not edit past this point ---------------------------

#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/correlation_plots")
correl_file_path <- file.path(getwd(), "output/correlation_plots")

#get freq data
freq_diff <- preprocessed_data %>% filter(reagent == "frequency",time == "d7",stimulation == "U") 
pvals <- pvals %>% filter(Tissue == selected_samp_type) #changed for Cameron's stats file

#sort the p-values in ascending order and select the top 10 hits
top_features <- pvals %>%
  arrange(padj_BH) %>%   #changed for Cameron's stats file
  slice(1:13) %>% #this used to be top 10, but decided to look at all 13 pops
  select(cell_pop, stim_cond,padj_BH) %>%
  distinct()

#-------------------- Pathology Score vs freq correlation plots --------------------

metadata_pathscore <- metadata %>% filter(Measure == "Pathology Score")

#merge metadata and frequency differences data
merged_data <- metadata_pathscore %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(freq_diff, by = "sampleID")
merged_data <- merged_data %>% select(-X,-Group,-reagent,-time) #clean up 

#filter merged data for top features
merged_top <- merged_data %>%
  filter(population %in% top_features$cell_pop & stimulation %in% top_features$stim_cond)
# ensure the Value column is numeric
merged_top$Value <- as.numeric(as.character(merged_top$Value))

#where plots will populate
dir.create("output/correlation_plots/freq_pathology_score")
correl_fp_pathscore <- file.path(getwd(), "output/correlation_plots/freq_pathology_score")

unique_cellpops <- unique(merged_top$population)
unique_stims <- unique(merged_top$stimulation)
color_scale <- scale_color_manual(values = c("Control" = "#F8766D", "Protein" = "#619CFF", "mRNA" = "#00BA38"))

# empty list to store p-values for correction
all_p_values <- c()
# empty data frame to store p and r values
correlation_results <- data.frame(
  cell_pop = character(),
  stim_cond = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

# collect all the p-values PRIOR to mult hypoth adjustment
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    # calculate the correlation and store the p-value
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    all_p_values <- c(all_p_values, cor_test$p.value)  # store p-values
  }
}

#apply Benjamini-Hochberg (BH) correction to all p-values
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

#go back thru the loop again to create plots using the adjusted p-values
index <- 1  # counter to keep track of adjusted p-values
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calculate the correlation coefficient (R value)
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    r_value <- cor_test$estimate
    
    #get the adjusted p-value
    p_value_adj <- adjusted_p_values[index]
    index <- index + 1
    
    # append results to df with r and p vals
    correlation_results <- rbind(correlation_results, data.frame(
      cell_pop = cellpop,
      stim_cond = stimcond,
      r_value = r_value,
      adjusted_p_value = p_value_adj
    ))
    
    # create plot
    p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
      geom_point(size = 3) +  
      geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
      labs(title = paste(cellpop, "_", selected_samp_type, sep = ""),
           x = "Pathology Score",
           y = paste(cellpop, "_", selected_samp_type, "_frequency", sep = "")) +
      color_scale +  
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey95"),
            plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous() +
      #annotate with R value and adjusted p-value
      annotate("text", x = mean(range(subset_data$Value, na.rm = TRUE)), 
               y = Inf, label = paste("R =", round(r_value, 2), "\nAdj p val =", format.pval(p_value_adj, digits = 2)), 
               vjust = 1.5, hjust = 0.5, size = 4.5)
    
    # save plot
    file_name <- paste("path_", selected_samp_type, "_", cellpop, ".png", sep = "") 
    correl_fp_pathscore_full <- file.path(correl_fp_pathscore, file_name)
    ggsave(correl_fp_pathscore_full, plot = p)
  }
}

#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "pathology_score"
correlation_results$tissue <- selected_samp_type
correlation_results1 <- correlation_results #used for combining csvs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/freq_pathology_score/",selected_samp_type,"_path_score_freq_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)


#-------------------- CoV RNA vs freq correlation plots --------------------

#note to self: for d7 N Gene log10 copies per 30 mg there is one n.t. value. Dave said to just get rid of this value
#note to self: -1 is the lower limit of detection 

metadata_covrna <- metadata %>% filter(Measure == "d7 N Gene log10 copies per 30 mg")
metadata_covrna <- metadata_covrna %>% filter(Value != "n.t.") #get rid of n.t. value

#merge metadata and frequency differences data
merged_data <- metadata_covrna %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(freq_diff, by = "sampleID")
merged_data <- merged_data %>% select(-X,-Group,-reagent,-time) #clean up 

#filter merged data for top features
merged_top <- merged_data %>%
  filter(population %in% top_features$cell_pop & stimulation %in% top_features$stim_cond)
# ensure the Value column is numeric
merged_top$Value <- as.numeric(as.character(merged_top$Value))

#where plots will populate
dir.create("output/correlation_plots/freq_cov_rna")
correl_fp_covrna <- file.path(getwd(), "output/correlation_plots/freq_cov_rna")

unique_cellpops <- unique(merged_top$population)
unique_stims <- unique(merged_top$stimulation)
color_scale <- scale_color_manual(values = c("Control" = "#F8766D", "Protein" = "#619CFF", "mRNA" = "#00BA38"))

#empty list to store p-values for all tests
all_p_values <- c()

# empty data frame to store p and r values
correlation_results <- data.frame(
  cell_pop = character(),
  stim_cond = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

#collect all the p-values prior to adjusting
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    # Ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calc and store the p-value
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    all_p_values <- c(all_p_values, cor_test$p.value)  # Store p-values
  }
}

#apply the Benjamini-Hochberg (BH) correction to all p-values
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

#go thru  loop again to create plots using the adjusted p-values
index <- 1  # counter to keep track of adjusted p-values
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calc correlation coefficient (R value)
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    r_value <- cor_test$estimate
    
    #pull the adjusted p-value
    p_value_adj <- adjusted_p_values[index]
    index <- index + 1
    
    # append results to df with r and p vals
    correlation_results <- rbind(correlation_results, data.frame(
      cell_pop = cellpop,
      stim_cond = stimcond,
      r_value = r_value,
      adjusted_p_value = p_value_adj
    ))
    
    #plot
    p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
      geom_point(size = 3) +  
      geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
      labs(title = paste(cellpop, "_", selected_samp_type, sep = ""),
           x = "d7 N Gene log10 copies per 30 mg",
           y = paste(cellpop, "_", selected_samp_type, "_frequency", sep = "")) +
      color_scale +  
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey95"),
            plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous() +
      #annotate plot with R value and adjusted p-value
      annotate("text", x = mean(range(subset_data$Value, na.rm = TRUE)), 
               y = Inf, label = paste("R =", round(r_value, 2), "\nAdj p val =", format.pval(p_value_adj, digits = 2)), 
               vjust = 1.5, hjust = 0.5, size = 4.5)
    
    #save plot
    file_name <- paste("covrna_", selected_samp_type, "_", cellpop, ".png", sep = "") 
    correl_fp_covrna_full <- file.path(correl_fp_covrna, file_name)
    ggsave(correl_fp_covrna_full, plot = p)
  }
}

#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "cov_RNA"
correlation_results$tissue <- selected_samp_type
correlation_results2 <- correlation_results #used for combining csvs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/freq_cov_rna/",selected_samp_type,"_cov_rna_freq_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)


#-------------------- Radiography Score d0-d7 vs freq correlation plots --------------------

metadata_radiogrd0tod7 <- metadata %>% filter(Measure == "Radiograph Score Day 0 to Day 7")

#merge metadata and frequency differences data
merged_data <- metadata_radiogrd0tod7 %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(freq_diff, by = "sampleID")
merged_data <- merged_data %>% select(-X,-Group,-reagent,-time) #clean up 

#filter merged data for top features
merged_top <- merged_data %>%
  filter(population %in% top_features$cell_pop & stimulation %in% top_features$stim_cond)
# ensure the Value column is numeric
merged_top$Value <- as.numeric(as.character(merged_top$Value))

#where plots will populate
dir.create("output/correlation_plots/freq_radiogrd0tod7")
correl_fp_rdgrph07 <- file.path(getwd(), "output/correlation_plots/freq_radiogrd0tod7")

unique_cellpops <- unique(merged_top$population)
unique_stims <- unique(merged_top$stimulation)
color_scale <- scale_color_manual(values = c("Control" = "#F8766D", "Protein" = "#619CFF", "mRNA" = "#00BA38"))

#empty list to store p-values for all tests
all_p_values <- c()

# empty data frame to store p and r values
correlation_results <- data.frame(
  cell_pop = character(),
  stim_cond = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

#collect all p-values before adjusting
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calc correlation and store p-value
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    all_p_values <- c(all_p_values, cor_test$p.value)  # Store p-values
  }
}

#apply the Benjamini-Hochberg (BH) correction to all p-values
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

#go thru loop again to create plots using the adjusted p-values
index <- 1  #counter to keep track of adjusted p-values
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calc the correlation coefficient (R value)
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    r_value <- cor_test$estimate
    
    #pull adjusted p-value for the current test
    p_value_adj <- adjusted_p_values[index]
    index <- index + 1
    
    # append results to df with r and p vals
    correlation_results <- rbind(correlation_results, data.frame(
      cell_pop = cellpop,
      stim_cond = stimcond,
      r_value = r_value,
      adjusted_p_value = p_value_adj
    ))
    
    #plot
    p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
      geom_point(size = 3) +  
      geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
      labs(title = paste(cellpop, "_", selected_samp_type, sep = ""),
           x = "Radiograph Score Day 0 to Day 7",
           y = paste(cellpop, "_", selected_samp_type, "_frequency", sep = "")) +
      color_scale +  
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey95"),
            plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous() +
      #annotate plot with R value and adjusted p-value
      annotate("text", x = mean(range(subset_data$Value, na.rm = TRUE)), 
               y = Inf, label = paste("R =", round(r_value, 2), "\nAdj p val =", format.pval(p_value_adj, digits = 2)), 
               vjust = 1.5, hjust = 0.5, size = 4.5)
    
    #save plot
    file_name <- paste("rdgrphd0d7_", selected_samp_type, "_", cellpop, ".png", sep = "") 
    correl_fp_rdgrph07_full <- file.path(correl_fp_rdgrph07, file_name)
    ggsave(correl_fp_rdgrph07_full, plot = p)
  }
}

#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "radiograph_d0tod7"
correlation_results$tissue <- selected_samp_type
correlation_results3 <- correlation_results #used for combining csvs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/freq_radiogrd0tod7/",selected_samp_type,"_radiogrd0tod7_freq_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)


#-------------------- Radiography Score d7 vs freq correlation plots --------------------

metadata_radiogrd7 <- metadata %>% filter(Measure == "Radiograph Score Day 7")

#merge metadata and frequency differences data
merged_data <- metadata_radiogrd7 %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(freq_diff, by = "sampleID")
merged_data <- merged_data %>% select(-X,-Group,-reagent,-time) #clean up 

#filter merged data for top features
merged_top <- merged_data %>%
  filter(population %in% top_features$cell_pop & stimulation %in% top_features$stim_cond)
# ensure the Value column is numeric
merged_top$Value <- as.numeric(as.character(merged_top$Value))

#where plots will populate
dir.create("output/correlation_plots/freq_radiogrd7")
correl_fp_rdgrph7 <- file.path(getwd(), "output/correlation_plots/freq_radiogrd7")

unique_cellpops <- unique(merged_top$population)
unique_stims <- unique(merged_top$stimulation)
color_scale <- scale_color_manual(values = c("Control" = "#F8766D", "Protein" = "#619CFF", "mRNA" = "#00BA38"))

#empty list to store p-values for all tests
all_p_values <- c()

# empty data frame to store p and r values
correlation_results <- data.frame(
  cell_pop = character(),
  stim_cond = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

#collect all p-values before adjusting
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calc correlation and store the p-value
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    all_p_values <- c(all_p_values, cor_test$p.value)  # Store p-values
  }
}

#apply the Benjamini-Hochberg (BH) correction to all p-values
adjusted_p_values <- p.adjust(all_p_values, method = "BH")

# go thru loop again to create plots using the adjusted p-values
index <- 1  # counter to keep track of adjusted p-values
for (cellpop in unique_cellpops) {
  for (stimcond in unique_stims) {
    subset_data <- subset(merged_top, population == cellpop & stimulation == stimcond)
    
    #ensure the subset data Value and feature columns are numeric
    subset_data$Value <- as.numeric(as.character(subset_data$Value))
    subset_data$feature <- as.numeric(as.character(subset_data$feature))
    
    #calculate correlation coefficient (R value)
    cor_test <- cor.test(subset_data$Value, subset_data$feature)
    r_value <- cor_test$estimate
    
    #get adjusted p-value for the current test
    p_value_adj <- adjusted_p_values[index]
    index <- index + 1
    
    # append results to df with r and p vals
    correlation_results <- rbind(correlation_results, data.frame(
      cell_pop = cellpop,
      stim_cond = stimcond,
      r_value = r_value,
      adjusted_p_value = p_value_adj
    ))
    
    #plot
    p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
      geom_point(size = 3) +  
      geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
      labs(title = paste(cellpop, "_", selected_samp_type, sep = ""),
           x = "Radiograph Score Day 7",
           y = paste(cellpop, "_", selected_samp_type, "_frequency", sep = "")) +
      color_scale +  
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_line(color = "grey95"),
            plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous() +
      #annotate plot with R value and adjusted p-value
      annotate("text", x = mean(range(subset_data$Value, na.rm = TRUE)), 
               y = Inf, label = paste("R =", round(r_value, 2), "\nAdj p val =", format.pval(p_value_adj, digits = 2)), 
               vjust = 1.5, hjust = 0.5, size = 4.5)
    
    #save
    file_name <- paste("rdgrphd7_", selected_samp_type, "_", cellpop, ".png", sep = "") 
    correl_fp_rdgrph7_full <- file.path(correl_fp_rdgrph7, file_name)
    ggsave(correl_fp_rdgrph7_full, plot = p)
  }
}

#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "radiograph_d7"
correlation_results$tissue <- selected_samp_type
correlation_results4 <- correlation_results #used for combining csvs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/freq_radiogrd7/",selected_samp_type,"_radiograph_d7_freq_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)

#combine all tests into one csv
merged_csv_all_test <- rbind(correlation_results1, correlation_results2, correlation_results3, correlation_results4)
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/",selected_samp_type,"_freq_corr_results_ALL.csv"))
write.csv(merged_csv_all_test, file = output_csv_path, row.names = FALSE)

