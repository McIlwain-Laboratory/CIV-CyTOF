#Cameron's Code for CIV Analysis
#Jan9th, 2025
#This code is mean to reproduce Fig2a and Fig2c in the CIV analysis
#These plots are both heatmaps of populations veresus vaccine group.
#The population values are log2/control
#This is meant to make sure that the coloring is correct more than anything else
#It should faithfully recreate these figures from CellEngine


#Upload the relevant packages
library(dplyr)
library(tidyverse)
library(patchwork)

# set the working directory to the folder where this analysis will take place
setwd("Replace me")


#Create directory space
dir.create("output")
dir.create("output/cell_abundance")
dir.create("output/cell_abundance/log2heatmaps")

#Load the dataset
preproc_folder_name <- "preprocessed data"
comb_file_path <- file.path(getwd(), preproc_folder_name, "combined_preprocessed.csv")
comb_cleaned_data <- read.csv(comb_file_path)
comb_cleaned_data <- comb_cleaned_data %>% filter(reagent == "frequency")
comb_cleaned_data <- comb_cleaned_data %>% select(-reagent, -X)
comb_cleaned_data <- comb_cleaned_data %>% filter(stimulation == "P") #filter for only stim P

#Calculate the mean_data
mean_data <- comb_cleaned_data %>%
  group_by(population, group, tissue_type) %>%
  summarise(mean_feature = mean(feature, na.rm = TRUE)) %>%
  ungroup()

#Pivot data and calculate Log2 ratio (Log2(x/control))
mean_data_pivoted <- mean_data %>% 
  pivot_wider(names_from = group, values_from = mean_feature) %>%
  mutate(Control_log2_control = log2(Control/Control)) %>%
  mutate(mRNA_log2_control = log2(mRNA/Control)) %>%
  mutate(Protein_log2_control = log2(Protein/Control)) %>% select(-mRNA, -Protein, -Control)

#Undo pivot
mean_data_final  <- mean_data_pivoted %>%
  pivot_longer(cols = c(Control_log2_control, Protein_log2_control, mRNA_log2_control),
               names_to = "group", values_to = "mean_feature")


#Make a heatmap
unique_tissue_vector <- unique(mean_data_final$tissue_type)


plot_list <- list()

for(tissue in unique_tissue_vector) {
  sub_df <- mean_data_final %>% filter(tissue_type == tissue)
  
  heatmap <- ggplot(sub_df, aes(y = population, x = group, fill = mean_feature )) +
    geom_tile(color = "black", linewidth = 0.5) +
    scale_fill_gradient2(low = "#0C6291", mid = "#FBFEF9", high ="#A63446", limits = c(-2.1, 2.1), midpoint = 0 ) +
    theme_minimal() +
    scale_x_discrete(name = tissue, labels = 
                       c("Control_log2_control" = "Control", 
                         "Protein_log2_control" = "Protein", 
                         "mRNA_log2_control" = "mRNA")
                     ) +
    scale_y_discrete(limits = c(
      "NK CD56loCD16+", "NK CD56+CD16-","cMC", "intMC", "ncMC", 
      "mDC", "pDC", "Bcell IgM+", "Bcell",
      "CD8T HLADR+Ki67+", "CD4T HLADR+Ki67+", "CD8T", "CD4T"
    ))
  
  plot_list <- append(plot_list, list(heatmap))
  
}

combined_plots <- wrap_plots(plot_list, ncol = 2)
print(combined_plots)

ggsave(filename = "output/cell_abundance/log2heatmaps/Fig2AandC.svg", plot = combined_plots, width = 10, height = 8)
