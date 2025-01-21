#Sasha's code for CIV analysis
#Dec2024
#Run script to produce correlation plots for top marker differences vs metadata 
#Fig 5b in manuscript
library(ggplot2)
require(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(stats)
library(readr)

#This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis is taking place
setwd("/Users/cameronmeikle/Desktop/UNR_RA_Position/Projects/CIV/CIV R Analysis 26July24")

# download metadata to that folder and paste the file path to it
metadata<- read.csv("/Users/cameronmeikle/Desktop/UNR_RA_Position/Projects/CIV/CIV R Analysis 26July24/metadata/CIV01 Lung Metadata RNA and Pathology.csv") 

#choose tissue to perform analysis: type either "PB" or "LN"
selected_samp_type <- "LN"

#choose stim to perform analysis: type either "P" or "R"
selected_stim <- "R"

#------------------- Do not edit past this point ---------------------------

#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/correlation_plots")
correl_file_path <- file.path(getwd(), "output/correlation_plots")

#get mark data
mark_diff <- preprocessed_data %>% filter(reagent != "frequency",time == "d7",stimulation == selected_stim) 
mark_diff <- mark_diff %>% mutate(pop_mark = paste(population, reagent, sep = "_"))
mark_diff <- mark_diff %>% select(-`X`)

#We are now doing top 10 markers that are shown in figure 4d
P_PB <- c("Bcell IgM+_IkBa",
          "Bcell_IkBa",
          "cMC_S6",
          "pDC_Erk1_2",
          "ncMC_Erk1_2",
          "NK CD56+CD16-_TBK1",
          "NK CD56+CD16-_P38",
  "intMC_PLCg2",
  "CD4T HLADR+Ki67+_MAPKAPK2",
  "NK CD56+CD16-_CREB")

P_LN <- c("Bcell IgM+_IkBa",
          "Bcell_IkBa",
          "cMC_S6",
          "pDC_Erk1_2",
          "ncMC_Erk1_2",
          "NK CD56+CD16-_TBK1",
          "NK CD56+CD16-_P38",
          "intMC_PLCg2",
          "CD4T HLADR+Ki67+_MAPKAPK2",
          "NK CD56+CD16-_CREB")

R_PB <- c("Bcell IgM+_IkBa",
          "Bcell_IkBa",
          "cMC_S6",
          "pDC_Erk1_2",
          "ncMC_Erk1_2",
          "NK CD56+CD16-_TBK1",
          "NK CD56+CD16-_P38",
          "intMC_PLCg2",
          "CD4T HLADR+Ki67+_MAPKAPK2",
          "NK CD56+CD16-_CREB")

R_LN <- c("Bcell IgM+_IkBa",
          "Bcell_IkBa",
          "cMC_S6",
          "pDC_Erk1_2",
          "ncMC_Erk1_2",
          "NK CD56+CD16-_TBK1",
          "NK CD56+CD16-_P38",
          "intMC_PLCg2",
          "CD4T HLADR+Ki67+_MAPKAPK2",
          "NK CD56+CD16-_CREB")

# create 4 separate dfs filtered for only rows where top three pop_marks are present
P_PB_mark_df <- mark_diff %>%
  filter(pop_mark %in% P_PB)

P_LN_mark_df <- mark_diff %>%
  filter(pop_mark %in% P_LN)

R_PB_mark_df <- mark_diff %>%
  filter(pop_mark %in% R_PB)

R_LN_mark_df <- mark_diff %>%
  filter(pop_mark %in% R_LN)

#selects a dataframe from above depending on input variables
df_name <- paste0(selected_stim,"_",selected_samp_type,"_mark_df")
df <- get(df_name)

#color scale for plots
color_scale <- scale_color_manual(values = c("Control" = "#FF4040", "Protein" = "#619CFF", "mRNA" = "#00BA38"))

#size of points for plots
point_size <- 3

#-------------------- Pathology Score vs marker correlation plots --------------------

metadata_pathscore <- metadata %>% filter(Measure == "Pathology Score")

#merge metadata and marker differences data
merged_data <- metadata_pathscore %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(df, by = "sampleID")
merged_data <- merged_data %>% select(-Group,-time) #clean up 

#where plots will populate
dir.create("output/correlation_plots/mark_pathology_score")
correl_mark_pathscore <- file.path(getwd(), "output/correlation_plots/mark_pathology_score")

unique_pop_marks <- unique(df$pop_mark)

#empty df to store correlation results
correlation_results <- data.frame(
  pop_mark = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

#empty df to store all p-values with given pop_mark
all_p_values <- data.frame(
  pop_mark = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

#loop thru each unique pop_mark to collect all the p-values PRIOR to mult hypoth adjustment
for (i in unique_pop_marks) {
  subset_data <- subset(merged_data, pop_mark == i)
  
  #ensure numeric
  subset_data$Value <- as.numeric(as.character(subset_data$Value))
  subset_data$feature <- as.numeric(as.character(subset_data$feature))
  
  #calc correlation and store p-val
  cor_test <- cor.test(subset_data$Value, subset_data$feature)
  
  #append givenpop_mark and p-value to the data frame
  all_p_values <- rbind(all_p_values, data.frame(pop_mark = i, p_value = cor_test$p.value))
}

print(all_p_values)

#apply Benjamini-Hochberg (BH) correction to all p-values
all_p_values$adjusted_p_value <- p.adjust(all_p_values$p_value, method = "BH")

print(all_p_values)

#create plots using the adjusted p-values
index <- 1
for (pop_mark in unique_pop_marks) {
  subset_data <- merged_data[merged_data$pop_mark == pop_mark, ]
  
  #ensure numeric
  subset_data$Value <- as.numeric(as.character(subset_data$Value))
  subset_data$feature <- as.numeric(as.character(subset_data$feature))
  
  #calculate corr coefficient R value
  cor_test <- cor.test(subset_data$Value, subset_data$feature)
  r_value <- cor_test$estimate
  
  #pull adjusted p-value
  p_value_adj <- all_p_values$adjusted_p_value[all_p_values$pop_mark == pop_mark]
  
  #append results to correl results df
  correlation_results <- rbind(correlation_results, data.frame(
    pop_mark = pop_mark,
    r_value = r_value,
    adjusted_p_value = p_value_adj
  ))
  
  #plot
  p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
    geom_point(size = point_size) +  
    geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
    labs(title = paste(pop_mark, "_", selected_samp_type,"_",selected_stim, sep = ""),
         x = "Pathology Score",
         y = paste(pop_mark, "_", selected_samp_type,"_",selected_stim, "_marker", sep = "")) +
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
  
  #save plot
  file_name <- paste("mark_pathscor_", selected_samp_type, "_", selected_stim,"_",pop_mark, ".svg", sep = "") 
  correl_mark_full <- file.path(correl_mark_pathscore, file_name)
  ggsave(correl_mark_full, plot = p)
}

#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "pathology_score"
correlation_results$tissue <- selected_samp_type
correlation_results1 <- correlation_results  # used for combining CSVs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/mark_pathology_score/", selected_samp_type, "_",selected_stim, "_path_score_mark_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)





#-------------------- Cov RNA vs marker correlation plots --------------------

metadata_covrna <- metadata %>% filter(Measure == "d7 N Gene log10 copies per 30 mg")

#merge metadata and marker differences data
merged_data <- metadata_covrna %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(df, by = "sampleID")
merged_data <- merged_data %>% select(-Group,-time) #clean up 

#where plots will populate
dir.create("output/correlation_plots/mark_cov_rna")
mark_cov_rna <- file.path(getwd(), "output/correlation_plots/mark_cov_rna")

unique_pop_marks <- unique(df$pop_mark)

#empty df to store correlation results
correlation_results <- data.frame(
  pop_mark = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

#empty df to store all p-values with given pop_mark
all_p_values <- data.frame(
  pop_mark = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

#loop thru each unique pop_mark to collect all the p-values PRIOR to mult hypoth adjustment
for (i in unique_pop_marks) {
  subset_data <- subset(merged_data, pop_mark == i)
  
  #ensure numeric
  subset_data$Value <- as.numeric(as.character(subset_data$Value))
  subset_data$feature <- as.numeric(as.character(subset_data$feature))
  
  #calc correlation and store p-val
  cor_test <- cor.test(subset_data$Value, subset_data$feature)
  
  #append givenpop_mark and p-value to the data frame
  all_p_values <- rbind(all_p_values, data.frame(pop_mark = i, p_value = cor_test$p.value))
}

print(all_p_values)

#apply Benjamini-Hochberg (BH) correction to all p-values
all_p_values$adjusted_p_value <- p.adjust(all_p_values$p_value, method = "BH")

print(all_p_values)

#create plots using the adjusted p-values
index <- 1
for (pop_mark in unique_pop_marks) {
  subset_data <- merged_data[merged_data$pop_mark == pop_mark, ]
  
  #ensure numeric
  subset_data$Value <- as.numeric(as.character(subset_data$Value))
  subset_data$feature <- as.numeric(as.character(subset_data$feature))
  
  #calculate corr coefficient R value
  cor_test <- cor.test(subset_data$Value, subset_data$feature)
  r_value <- cor_test$estimate
  
  #pull adjusted p-value
  p_value_adj <- all_p_values$adjusted_p_value[all_p_values$pop_mark == pop_mark]
  
  #append results to correl results df
  correlation_results <- rbind(correlation_results, data.frame(
    pop_mark = pop_mark,
    r_value = r_value,
    adjusted_p_value = p_value_adj
  ))
  
  #plot
  p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
    geom_point(size = point_size) +  
    geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
    labs(title = paste(pop_mark, "_", selected_samp_type,"_",selected_stim, sep = ""),
         x = "d7 N Gene log10 copies per 30 mg",
         y = paste(pop_mark, "_", selected_samp_type,"_",selected_stim, "_marker", sep = "")) +
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
  
  #save plot
  file_name <- paste("mark_cov_rna_", selected_samp_type, "_", selected_stim,"_",pop_mark, ".svg", sep = "") 
  correl_mark_full <- file.path(mark_cov_rna, file_name)
  ggsave(correl_mark_full, plot = p)
}

#did not edit bellow
#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "d7 N Gene log10 copies per 30 mg"
correlation_results$tissue <- selected_samp_type
correlation_results1 <- correlation_results  # used for combining CSVs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/mark_cov_rna/", selected_samp_type, "_",selected_stim, "_cov_rna_mark_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)




#-------------------- Radiograph Score Day 7 vs marker correlation plots --------------------

metadata_radiogrd7 <- metadata %>% filter(Measure == "Radiograph Score Day 7")

#merge metadata and marker differences data
merged_data <- metadata_radiogrd7 %>%
  rename(sampleID = `Animal.ID`) %>%
  left_join(df, by = "sampleID")
merged_data <- merged_data %>% select(-Group,-time) #clean up 

#where plots will populate
dir.create("output/correlation_plots/mark_radiogrd7")
mark_radiogrd7 <- file.path(getwd(), "output/correlation_plots/mark_radiogrd7")

unique_pop_marks <- unique(df$pop_mark)

#empty df to store correlation results
correlation_results <- data.frame(
  pop_mark = character(),
  r_value = numeric(),
  adjusted_p_value = numeric(),
  stringsAsFactors = FALSE
)

#empty df to store all p-values with given pop_mark
all_p_values <- data.frame(
  pop_mark = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

#loop thru each unique pop_mark to collect all the p-values PRIOR to mult hypoth adjustment
for (i in unique_pop_marks) {
  subset_data <- subset(merged_data, pop_mark == i)
  
  #ensure numeric
  subset_data$Value <- as.numeric(as.character(subset_data$Value))
  subset_data$feature <- as.numeric(as.character(subset_data$feature))
  
  #calc correlation and store p-val
  cor_test <- cor.test(subset_data$Value, subset_data$feature)
  
  #append givenpop_mark and p-value to the data frame
  all_p_values <- rbind(all_p_values, data.frame(pop_mark = i, p_value = cor_test$p.value))
}

print(all_p_values)

#apply Benjamini-Hochberg (BH) correction to all p-values
all_p_values$adjusted_p_value <- p.adjust(all_p_values$p_value, method = "BH")

print(all_p_values)

#create plots using the adjusted p-values
index <- 1
for (pop_mark in unique_pop_marks) {
  subset_data <- merged_data[merged_data$pop_mark == pop_mark, ]
  
  #ensure numeric
  subset_data$Value <- as.numeric(as.character(subset_data$Value))
  subset_data$feature <- as.numeric(as.character(subset_data$feature))
  
  #calculate corr coefficient R value
  cor_test <- cor.test(subset_data$Value, subset_data$feature)
  r_value <- cor_test$estimate
  
  #pull adjusted p-value
  p_value_adj <- all_p_values$adjusted_p_value[all_p_values$pop_mark == pop_mark]
  
  #append results to correl results df
  correlation_results <- rbind(correlation_results, data.frame(
    pop_mark = pop_mark,
    r_value = r_value,
    adjusted_p_value = p_value_adj
  ))
  
  #plot
  p <- ggplot(subset_data, aes(x = Value, y = feature, color = group)) +
    geom_point(size = point_size) +  
    geom_smooth(method = "lm", se = TRUE, color = "darkgrey", fill = "grey", size = 0.5) +
    labs(title = paste(pop_mark, "_", selected_samp_type,"_",selected_stim, sep = ""),
         x = "Radiograph Score Day 7",
         y = paste(pop_mark, "_", selected_samp_type,"_",selected_stim, "_marker", sep = "")) +
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
  
  #save plot
  file_name <- paste("mark_radiogrd7_", selected_samp_type, "_", selected_stim,"_",pop_mark, ".svg", sep = "") 
  correl_mark_full <- file.path(mark_radiogrd7, file_name)
  ggsave(correl_mark_full, plot = p)
}

#save R values and adjusted p-values to a CSV file
correlation_results$metadata_test <- "Radiograph Score Day 7"
correlation_results$tissue <- selected_samp_type
correlation_results1 <- correlation_results  # used for combining CSVs at the end
output_csv_path <- file.path(getwd(), paste0("output/correlation_plots/mark_radiogrd7/", selected_samp_type, "_",selected_stim, "_radiogrd7_mark_corr_results.csv"))
write.csv(correlation_results, file = output_csv_path, row.names = FALSE)



