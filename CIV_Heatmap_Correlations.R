#Cameron Meikle
#December 31st, 2024
#This script creates the heatmaps for the correlation plots. Very simply the correlation
#heatmaps will have three columns, one for each of the 

#Library in packages
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggpattern)


#Set working directory
setwd("Replace me")

working_directory <- getwd()

#Read in files for frequency, filter, and combine
LN_Freq <- read.csv(paste0(working_directory, "/output/correlation_plots/LN_freq_corr_results_ALL.csv"))
PB_Freq <- read.csv(paste0(working_directory, "/output/correlation_plots/PB_freq_corr_results_ALL.csv"))

combined_Freq_data <- rbind(LN_Freq, PB_Freq) %>% filter(metadata_test != "radiograph_d0tod7") %>%
  mutate(significance = ifelse(adjusted_p_value < 0.05, "significant(p<0.05)", "p>0.05"))

#Create basic ggplot heatmap that is facet_wrapped by tissue type
freq_plot <- ggplot(combined_Freq_data, aes(y = cell_pop, x = metadata_test, fill = r_value, pattern = significance)) + 
  geom_tile_pattern(color = "black", lwd = 0.5, linetype = 1,
                    pattern_alpha = 0.7, pattern_size = 0.05, 
                    pattern_spacing = 0.15, pattern_colour = "black",
                    pattern_fill = "grey") +
  scale_pattern_discrete(choices = c("crosshatch", "none")) +
  facet_wrap(~ tissue) +
  scale_fill_gradient2(low = "#0C6291", mid = "#FBFEF9", high ="#A63446", limits = c(-1, 1), midpoint = 0 ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave("output/correlation_plots/freq_corr_heatmap.pdf", plot = freq_plot, dpi = 300, height = 8, width = 6)

#Read in files for markers, filter, and combine
#RNA_Cov_marker
PB_R_cov_rna_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_cov_rna/PB_R_cov_rna_mark_corr_results.csv"))
PB_R_cov_rna_mark$stim <- "R"

PB_P_cov_rna_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_cov_rna/PB_P_cov_rna_mark_corr_results.csv"))
PB_P_cov_rna_mark$stim <- "P"

LN_R_cov_rna_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_cov_rna/LN_R_cov_rna_mark_corr_results.csv"))
LN_R_cov_rna_mark$stim <- "R"

LN_P_cov_rna_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_cov_rna/LN_P_cov_rna_mark_corr_results.csv"))
LN_P_cov_rna_mark$stim <- "P"

#Pathology Score
PB_R_path_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_pathology_score/PB_R_path_score_mark_corr_results.csv"))
PB_R_path_mark$stim <- "R"

PB_P_path_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_pathology_score/PB_P_path_score_mark_corr_results.csv"))
PB_P_path_mark$stim <- "P"

LN_R_path_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_pathology_score/LN_R_path_score_mark_corr_results.csv"))
LN_R_path_mark$stim <- "R"

LN_P_path_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_pathology_score/LN_P_path_score_mark_corr_results.csv"))
LN_P_path_mark$stim <- "P"

#Radiography day 7 score
PB_R_radiod7_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_radiogrd7/PB_R_radiogrd7_mark_corr_results.csv"))
PB_R_radiod7_mark$stim <- "R"

PB_P_radiod7_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_radiogrd7/PB_P_radiogrd7_mark_corr_results.csv"))
PB_P_radiod7_mark$stim <- "P"

LN_R_radiod7_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_radiogrd7/LN_R_radiogrd7_mark_corr_results.csv"))
LN_R_radiod7_mark$stim <- "R"

LN_P_radiod7_mark <- read.csv(paste0(working_directory, "/output/correlation_plots/mark_radiogrd7/LN_P_radiogrd7_mark_corr_results.csv"))
LN_P_radiod7_mark$stim <- "P"


#Combine all of the dataframes into marker dataframe and rename metadata
combined_Marker_data <- rbind(PB_P_cov_rna_mark, PB_R_cov_rna_mark, LN_P_cov_rna_mark, LN_R_cov_rna_mark, 
                             PB_P_path_mark, PB_R_path_mark, LN_P_path_mark, LN_R_path_mark, 
                             PB_P_radiod7_mark, PB_R_radiod7_mark, LN_P_radiod7_mark, LN_R_radiod7_mark) %>%
  mutate(metadata_test = recode(metadata_test, "d7 N Gene log10 copies per 30 mg" = "cov_RNA", 
                                "Radiograph Score Day 7" = "Radiograph_d7")) %>%
  mutate(significance = ifelse(adjusted_p_value < 0.05, "significant(p<0.05)", "p>0.05"))


#Create vectors for for loops
stim_vector <- unique(combined_Marker_data$stim)
tissue_vector <- unique(combined_Marker_data$tissue)

#Create a list to store plots in for further work
plot_list <- list()

#Nested For loop structure for running 
for(tissue_type in tissue_vector){
  for(stim_type in stim_vector){
  subset_df <- combined_Marker_data %>% filter(tissue == tissue_type & stim == stim_type)
  
  plot <- ggplot(subset_df, aes(y = pop_mark, x = metadata_test, fill = r_value, pattern = significance)) +
    geom_tile_pattern(color = "black", lwd = 0.5, linetype = 1,
                      pattern_alpha = 0.7, pattern_size = 0.05, 
                      pattern_spacing = 0.15, pattern_colour = "black",
                      pattern_fill = "grey") +
    scale_pattern_manual(values = c("p>0.05" = "crosshatch", "significant(p<0.05)" = "none")) +
    scale_fill_gradient2(low = "#0C6291", mid = "#FBFEF9", high ="#A63446", limits = c(-1, 1), midpoint = 0 ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  print(plot)
  plot_list <- append(plot_list, list(plot))
  
  }
}

#Combine listed plots with patchwork
combined_plots <- wrap_plots(plot_list, ncol = 2)
print(combined_plots)

ggsave(filename = "output/correlation_plots/marker_corr_heatmap.pdf", plot = combined_plots, width = 10, height = 6, dpi = 300)
