#Sasha's code for CIV analysis
#Nov2024
#Run this for producing frequency/abundance graphs comparing between tissues (PB and LN)
library(ggplot2)
require(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(stats)
library("ggpubr")
library(ggrepel)

# This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis is taking place
setwd("Replace me")

#------------------- Do not edit past this point ---------------------------

#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
PB_file_path <- file.path(getwd(), preproc_folder_name, paste("PB_preprocessed.csv", sep = ""))
LN_file_path <- file.path(getwd(), preproc_folder_name, paste("LN_preprocessed.csv", sep = ""))
preprocessed_PB_data <- read.csv(PB_file_path)
preprocessed_LN_data <- read.csv(LN_file_path)
dir.create("output")
dir.create("output/comparing_tissues")
dir.create("output/comparing_tissues/frequency")
dir.create("output/comparing_tissues/frequency/log_scaled_plots")
dir.create("output/comparing_tissues/frequency/unscaled_plots")
freq_tissue_file_path <- file.path(getwd(), "output/comparing_tissues/frequency")
freq_tissue_file_path_log <- file.path(getwd(), "output/comparing_tissues/frequency/log_scaled_plots")
freq_tissue_file_path_nolog <- file.path(getwd(), "output/comparing_tissues/frequency/unscaled_plots")

#-------------------- Task 6: Cell Abundances Tissue Comparison for PB vs LN--------------------

#snippet below reads in frequency values then merges PB and LN into one data frame
#read in frequency values only
PB_freq_data <- preprocessed_PB_data %>% filter(reagent == "frequency")
LN_freq_data <- preprocessed_LN_data %>% filter(reagent == "frequency")

PB_freq_data <- PB_freq_data %>% filter(stimulation == "U", time == "d7")
LN_freq_data <- LN_freq_data %>% filter(stimulation == "U", time == "d7")

#rename columns
colnames(PB_freq_data)[8] ="feature.PB"
colnames(LN_freq_data)[8] ="feature.LN"

#drop the useless reagent column
PB_freq_data <- subset(PB_freq_data, select = -c(reagent))
LN_freq_data <- subset(LN_freq_data, select = -c(reagent))

#basically Vlookup by the following parameters and add PB feature values to matching parameters 
freqdata_by_col <- merge(LN_freq_data, PB_freq_data, by = c("population", "group","sampleID","time","stimulation"), all.x = TRUE)
#get rid of useless row# columns to make final freq dataframe to work with
freqdata_by_col <- freqdata_by_col %>% select(-X.x,-X.y,-time)

#find averages and add ANOVA statistics
#calculate averages for each combination of population, stimulation, and group
freq_average_data <- freqdata_by_col %>%
  group_by(population, stimulation, group) %>%
  summarise(mean_PB = mean(feature.PB), 
            mean_LN = mean(feature.LN), 
            .groups = 'drop')

#empty data frame to store p-values
freq_p_values <- data.frame(population = character(), stimulation = character(), 
                            p_value_PB = numeric(), p_value_LN = numeric())

#nested loop for ANOVA tests
for (pop in unique(freqdata_by_col$population)) {
  for (stim in unique(freqdata_by_col$stimulation)) {
    subset_data <- freqdata_by_col[freqdata_by_col$population == pop & freqdata_by_col$stimulation == stim, ]
    
    if (nrow(subset_data) > 1) {
      p_value_PB <- summary(aov(feature.PB ~ group, data = subset_data))[[1]]$'Pr(>F)'[1]
      p_value_LN <- summary(aov(feature.LN ~ group, data = subset_data))[[1]]$'Pr(>F)'[1]
    } else {
      p_value_PB <- NA
      p_value_LN <- NA
    }
    
    freq_p_values <- rbind(freq_p_values, data.frame(population = pop, stimulation = stim, 
                                                     p_value_PB = p_value_PB, p_value_LN = p_value_LN))
  }
}

# Apply Benjamini-Hochberg correction to p-values
freq_p_values$p_value_PB <- p.adjust(freq_p_values$p_value_PB, method = "BH")
freq_p_values$p_value_LN <- p.adjust(freq_p_values$p_value_LN, method = "BH")

#merge the average data with corrected p-values
freqdata_by_col_final <- left_join(freq_average_data, freq_p_values, by = c("population", "stimulation"))

#add column for whether significance in EITHER PB or LN
freqdata_by_col_final$significance <- ifelse(freqdata_by_col_final$p_value_PB < 0.05 | freqdata_by_col_final$p_value_LN < 0.05, "yes", "no")

#Hardcoding in significance for pDCs. They are significant via a t-test between mRNA and control in only PB 
freqdata_by_col_final[37, "significance"] <- "yes"
freqdata_by_col_final[38, "significance"] <- "yes"
freqdata_by_col_final[39, "significance"] <- "yes"

#output csv with averages and p values
write.csv(freqdata_by_col_final, file.path(freq_tissue_file_path, "LN_and_PB_freq_avgs_with_p_val.csv"), row.names = FALSE)

#unique_timept<- unique(freqdata_by_col_averages$time) #there is actually only d7 in the dataset because LN was only collected day7
unique_stim <- unique(freqdata_by_col_final$stimulation)

#defien color scheme for cell populations:
# Define HEX color scheme for cell populations
cell_pop_colorscheme <- c(
  "Bcell" = "#f86d87",
  "Bcell IgM+" = "#e18a01",
  "CD4T" = "#be9c00",
  "CD4T HLADR+Ki67+" = "#8cab00",
  "CD8T" = "#24b700",
  "CD8T HLADR+Ki67+" = "#1bbe70",
  "cMC" = "#27c1ac",
  "intMC" = "#31bbda",
  "mDC" = "#38adfc",
  "ncMC" = "#8b93ff",
  "NK CD56+CD16-" = "#d575fe",
  "NK CD56loCD16+" = "#f962dd",
  "pDC" = "#fc64ac"
)


#--------------------------------
#This section does some simple manipulation of freqdata_by_col_final to run a pearsons correlation on 
#the mean data for a given cell population frequency between lymphnode and PB.
correlation_LNvsPB <- freqdata_by_col_final %>% select(population, stimulation, group, mean_PB, mean_LN)

group_vector <- unique(correlation_LNvsPB$group)

correlation_LNvsPB_df <- data.frame(group = character(), r_value = numeric(), p_value = numeric())

for(grouping in group_vector){
  #Create subset df to filter by grouping
  subset_df <- correlation_LNvsPB %>% filter(group == grouping)
  
  #run pearsons
  result <- cor.test(subset_df$mean_PB, subset_df$mean_LN, method = "pearson")
  
  #make results into dataframe and add to overall dataframe.
  result_df <-data.frame(group = grouping, r_value = result$estimate, p_value = result$p.value)
  correlation_LNvsPB_df <- rbind(correlation_LNvsPB_df, result_df)
}

#do the same for the entire set of everything
overall_results_LNvsPB <- cor.test(correlation_LNvsPB$mean_PB, correlation_LNvsPB$mean_LN, method = "pearson")

correlation_LNvsPB_df <- add_row(correlation_LNvsPB_df, group = "all", 
                                 r_value = overall_results_LNvsPB$estimate, 
                                 p_value = overall_results_LNvsPB$p.value)

#These values are used in the manuscript while describing figure 3b.


#------abundance plots (log scaled) V1 ------
#plot avg cell pop abundance values. Connect into a triangle if significant for either t-test or ANOVA (all ANOVAs are included in code and t-test for pDCs is hardcoded in)

subset_data <- freqdata_by_col_final %>% filter(stimulation == "U")
  
  # Use shapes that support fill and color (border)
shape_mapping <- c("Control" = 21, "Protein" = 22, "mRNA" = 24)
  
  # Filter for statistically significant Control group data points
significant_control_data <- subset_data %>% filter(group == "Control", significance == "yes")
  
  p1 <- ggplot() +
    geom_point(data = subset_data, aes(x = mean_PB, y = mean_LN, fill = population, shape = group), color = "black", size = 3, stroke = 0.5) +
    scale_shape_manual(values = shape_mapping) +
    scale_fill_manual(values = cell_pop_colorscheme) +
    geom_polygon(data = subset_data %>% filter(significance == "yes"), aes(x = mean_PB, y = mean_LN, group = population, fill = population), alpha = 0.4) +
    geom_text_repel(data = significant_control_data, aes(x = mean_PB, y = mean_LN, label = population), size = 3, color = "black") + # Label for significant Control group points
    scale_x_log10() +
    scale_y_log10() +
    ggtitle(paste("Frequency by population", stim, "(log10, p<0.05 shaded, BH corrected)", sep = " ")) +
    theme_minimal() +
    labs(x = "Average Freq PB (log10)", y = "Average Freq LN (log10)", color = "Cell Populations Color Key", fill = "Stat Signif Cell Populations", shape = "Group") +
    guides(fill = guide_legend(title = "Cell Populations Color Key", ncol = 2)) +
    theme(legend.position = "right")
  
  ggsave(paste("freq_triangle_plot_", stim, "_log10", ".svg", sep = ""), 
         path = freq_tissue_file_path_log, 
         plot = p1, 
         width = 10, height = 8, units = "in", bg = "white") 


#------abundance plots (log scaled) V2 ------
#abundance plots (log scaled)
for (stim in unique_stim) {
  subset_data <- freqdata_by_col_final %>% filter(stimulation == stim)
  
  # Use shapes that support fill and color (border)
  shape_mapping <- c("Control" = 21, "Protein" = 22, "mRNA" = 24)
  
  # Filter for statistically significant Control group data points
  significant_control_data <- subset_data %>% filter(group == "Control", significance == "yes")
  
  p1 <- ggplot() +
    geom_point(data = subset_data, aes(x = mean_PB, y = mean_LN, fill = population, shape = group), color = "black", size = 3, stroke = 0.5) +
    scale_shape_manual(values = shape_mapping) +
    scale_fill_manual(values = cell_pop_colorscheme) +
    geom_polygon(data = subset_data %>% filter(significance == "yes"), aes(x = mean_PB, y = mean_LN, group = population, fill = population), alpha = 0.4) +
    geom_text_repel(data = significant_control_data, aes(x = mean_PB, y = mean_LN, label = population), size = 3, color = "black") +
    scale_x_log10() +
    scale_y_log10() +
    ggtitle(paste("Frequency by population", stim, "(log10, p<0.05 shaded, BH corrected)", sep = " ")) +
    theme_minimal() +
    labs(x = "Average Freq PB (log10)", y = "Average Freq LN (log10)", color = "Cell Populations Color Key", fill = "Stat Signif Cell Populations", shape = "Group") +
    guides(fill = guide_legend(title = "Cell Populations Color Key", override.aes = list(color = NA, shape = 22, alpha = 1))) + # Override alpha for legend
    theme(legend.position = "right")
  
  ggsave(paste("freq_triangle_plot_", stim, "_log10", ".svg", sep = ""), 
         path = freq_tissue_file_path_log, 
         plot = p1, 
         width = 10, height = 8, units = "in", bg = "white") 
}

#abundance plots (NOT log scaled)
#for a given stim condition, plot avg cell pop abundance values. Connect into a triangle if significant
for (stim in unique_stim) {
  subset_data <- freqdata_by_col_final %>% filter(stimulation == stim)
  
  # Use shapes that support fill and color (border)
  shape_mapping <- c("Control" = 21, "Protein" = 22, "mRNA" = 24)
  
  # Filter for statistically significant Control group data points
  significant_control_data <- subset_data %>% filter(group == "Control", significance == "yes")
  
  p1 <- ggplot() +
    geom_point(data = subset_data, aes(x = mean_PB, y = mean_LN, fill = population, shape = group), color = "black", size = 3, stroke = 0.5) +
    scale_shape_manual(values = shape_mapping) +
    scale_fill_manual(values = cell_pop_colorscheme) +
    geom_polygon(data = subset_data %>% filter(significance == "yes"), aes(x = mean_PB, y = mean_LN, group = population, fill = population), alpha = 0.4) +
    geom_text_repel(data = significant_control_data, aes(x = mean_PB, y = mean_LN, label = population), size = 3, color = "black") + # Label for significant Control group points
    xlim(0, 50) + # Keeping the custom x-limits
    ggtitle(paste("Frequency by population", stim, "(p<0.05 shaded, BH corrected)", sep = " ")) +
    theme_minimal() +
    labs(x = "Average Freq PB", y = "Average Freq LN", color = "Cell Populations Color Key", fill = "Stat Signif Cell Populations", shape = "Group") +
    guides(fill = guide_legend(title = "Cell Populations Color Key", override.aes = list(color = NA, shape = 22, alpha = 1))) +
    theme(legend.position = "right")
  
  ggsave(paste("freq_triangle_plot_", stim, "_nolog", ".svg", sep = ""), 
         path = freq_tissue_file_path_nolog, 
         plot = p1, 
         width = 10, height = 8, units = "in", bg = "white") 
}

