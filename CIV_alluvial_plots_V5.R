#Sasha's code for CIV analysis
#Nov2024
#After you run Arcsinh_IH_with_sasha_mods.R script, run this for alluvial plots for Fig3a
library(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(readr)

# This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis will take place
setwd("/Users/cameronmeikle/Desktop/UNR_RA_position/Projects/CIV/CIV R Analysis 26July24")

#------------------- Do not edit past this point ----

#script below creates an output folder, and reads in the relevant preprocessed csv
dir.create("output")
dir.create("output/cell_abundance")
dir.create("output/cell_abundance/alluvial_plots")
preproc_folder_name <- "preprocessed data"
comb_file_path <- file.path(getwd(), preproc_folder_name, "combined_preprocessed.csv")
comb_cleaned_data <- read.csv(comb_file_path)
comb_cleaned_data <- comb_cleaned_data %>% filter(reagent == "frequency")
comb_cleaned_data <- comb_cleaned_data %>% select(-reagent, -X)
comb_cleaned_data <- comb_cleaned_data %>% filter(stimulation == "U") #filter for only stim P
file_name <- "comb_data_for_alluvial.csv"
combi_file_path <- file.path(getwd(), preproc_folder_name, file_name)
write.csv(comb_cleaned_data, combi_file_path, row.names = FALSE)

#group data by population, vaccine group, and tiss type
mean_data <- comb_cleaned_data %>%
  group_by(population, group, tissue_type) %>%
  summarise(mean_feature = mean(feature, na.rm = TRUE)) %>%
  ungroup()

#Correct the values for B Cells, CD4T cells, and CD8T cells due to nesting of 
#B cell IgM+, CD4T HLADR+Ki67+, and CD8T HLADR+Ki67+ cells

mean_data_pivoted <- mean_data %>% pivot_wider(names_from = population, values_from = mean_feature) %>%
  mutate(Bcell = Bcell - `Bcell IgM+`) %>% mutate(CD4T = CD4T - `CD4T HLADR+Ki67+`) %>%
  mutate(CD8T = CD8T - `CD8T HLADR+Ki67+`) 

#Get everything scaled to 100 percent per tissue_type, vaccine group box, not using but can be used
#scale_to_100 <- function(row) {
#  row_sum <- sum(row)
#  if (row_sum == 0) {
#    return(row)  # If the row sum is 0, return it unchanged
#  }
#  return(row * (100 / row_sum))
#}

#apply to only numeric columns
#numeric_cols <- sapply(mean_data_pivoted, is.numeric)

#mean_data_pivoted[numeric_cols] <- t(apply(mean_data_pivoted[numeric_cols], 1, scale_to_100))


#Change back into mean_data setup

mean_data <- mean_data_pivoted %>% 
  pivot_longer(cols = Bcell:pDC, names_to = "population", values_to = "mean_feature")

#Just to check the total percentage of frequency that is actually plotted. Should be 100% per group
df_correction_hypothesis <- mean_data %>% group_by(group, tissue_type) %>% summarise(total_freq = sum(mean_feature, na.rm = TRUE))


#write grouped dataframe to CSV
mean_data_file_name <- "mean_data_for_alluvial_stim_P.csv"
meandata_file_path <- file.path(getwd(),"output/cell_abundance/alluvial_plots", mean_data_file_name)
write.csv(mean_data, meandata_file_path, row.names = FALSE)

#alluvial color scheme
alluvial_color_scheme <- c(
  "LN" = "#ba02cf",
  "PB" = "#d6c95c",
  "Control" = "#FF4040",
  "mRNA" = "#00BA38",
  "Protein" = "#619CFF",
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

#make plot
alluv_plot <- ggplot(data = mean_data,
             aes(axis1 = tissue_type, axis2 = group, axis3 = population, y = mean_feature)) +
  scale_x_discrete(limits = c("Tissue Type", "Group", "Population"), expand = c(.2, .05)) +
  xlab("") +
  geom_flow(aes(fill = tissue_type), stat = "alluvium", lode.guidance = "rightleft") +
  geom_stratum(aes(fill = tissue_type), show.legend = FALSE) +  
  geom_stratum(aes(fill = group), show.legend = FALSE) +        
  geom_stratum(aes(fill = population), show.legend = FALSE) +   
  scale_fill_manual(values = alluvial_color_scheme) +  
  geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum))) +
  theme_minimal(base_size = 12) +
  ggtitle("Distribution of cells between tissue type, vaccine group, and cell populations") +
  ylab("Frequency") +
  guides(fill = "none")

file_name <- "alluvial_stim_U.svg" 
alluv_plot_path <- file.path(getwd(),"output/cell_abundance/alluvial_plots", file_name)
ggsave(alluv_plot_path, plot = alluv_plot)

