#Sasha's code for CIV analysis
#Mar2024
#Run this for producing volcano plots for signaling responses
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

# download raw cell engine data to that folder and paste the file path to it
cell_engine_data <- read.csv("Replace me")

#------------------- Do not edit past this point (unless adjusting parameters)---------------------------

#areas for potential adjustments are marked with comment

#filtering CellEngine data for only d7 and only certain cell pops
cell_engine_data <- cell_engine_data %>% filter(Day == "d7") #added on 24May2024 because the stats scripts don't calculate stats properly if other timepts remain included
cell_engine_data <- cell_engine_data %>% filter(population %in% c("Bcell", 
                                                                  "Bcell IgM+", 
                                                                  "CD4T", 
                                                                  "CD4T HLADR+Ki67+", 
                                                                  "CD8T", 
                                                                  "CD8T HLADR+Ki67+", 
                                                                  "cMC", 
                                                                  "intMC", 
                                                                  "mDC", 
                                                                  "ncMC", 
                                                                  "NK CD56+CD16-", 
                                                                  "NK CD56loCD16+", 
                                                                  "pDC"))


#========================= Data Prep for Fig 4a Volcano plots ===================
selected_samp_type <- "PB"

#Custom function to filter CE data
filter_my_data <- function(unfiltered_data,stringpath) {
  unfiltered_data$value[is.na(unfiltered_data$value)] <- 0 #all NaN's in the value column are replaced with zeros
  unfiltered_data <- unfiltered_data %>% filter(SampleType == selected_samp_type)
  unfiltered_data <- unfiltered_data %>% filter(!is.na(VaccineGroup))
  median_df <- unfiltered_data %>% filter(statistic == "median")
  freq_df <- unfiltered_data %>% filter(statistic == "percent")
  filtered_median_df <- median_df %>% select(population, reagent, group = VaccineGroup, sampleID = SampleID, stimulation = Stimulation, time = Day, median = value)
  filtered_freq_df <- freq_df %>% select(population, reagent, group = VaccineGroup, sampleID = SampleID, stimulation = Stimulation, time = Day, frequency = value)
  mediansave <- paste0(stringpath,"/volc_CIVfiltered_median_df.csv")
  freqsave <- paste0(stringpath,"/volc_CIVfiltered_freq_df.csv")
  write.csv(filtered_median_df, mediansave,row.names = FALSE)
  write.csv(filtered_freq_df, freqsave,row.names = FALSE) 
}

#this is where the filtered and preprocessed data will go
dir.create("preprocessed data for volc plot")
current_directory <- getwd()
filtered_input_data <- file.path(current_directory, "preprocessed data for volc plot")

# execute custom function using variables defined above
filter_my_data(cell_engine_data,filtered_input_data)

folder_name <- "preprocessed data for volc plot"
median_file_name <- "volc_CIVfiltered_median_df.csv"
funct_path <- file.path(current_directory, folder_name, median_file_name)

freq_file_name <- "volc_CIVfiltered_freq_df.csv"
freq_path <- file.path(current_directory, folder_name, freq_file_name)

out_path <- file.path(current_directory, folder_name)

### Features
MyData = read.csv(funct_path) 
final_data = c()
fil = setdiff(colnames(MyData), c("stimulation", "median"))
  
  stims = unique(MyData$stimulation) 
  
  MyData_timepoint = MyData %>% 
    mutate(feature = asinh(median/5)) %>%  #arc sign transformation (if you change this, remember to also change for LN)
    #mutate(feature = (median)) %>% 
    select(-one_of("median")) %>%
    pivot_wider(names_from = stimulation, values_from = feature)
  
  for (stim in stims){
    if (stim!="U"){
      MyData_timepoint[stim] = MyData_timepoint[stim] - MyData_timepoint["U"] #(if you change this, remember to also change for LN)
    }
  }
  
  MyData_fin = MyData_timepoint %>% 
    pivot_longer(-all_of(fil),
                 names_to="stimulation",values_to="feature")
  
  final_data = rbind(final_data, MyData_fin)
  
freq = read.csv(freq_path) %>% rename(feature = frequency) %>% mutate(reagent = "frequency")
final_data = rbind(final_data, freq)
final_data$tissue_type <- selected_samp_type
final_data <- final_data %>% filter(time == "d7")
csv_file_path <- file.path(out_path, paste(selected_samp_type, "volc_preprocessed.csv", sep = "_"))
write.csv(final_data, csv_file_path)

#repeat of everything above but for LN tissue type 

selected_samp_type <- "LN"

#Custom function to filter CE data
filter_my_data <- function(unfiltered_data,stringpath) {
  unfiltered_data$value[is.na(unfiltered_data$value)] <- 0 #all NaN's in the value column are replaced with zeros
  unfiltered_data <- unfiltered_data %>% filter(SampleType == selected_samp_type)
  unfiltered_data <- unfiltered_data %>% filter(!is.na(VaccineGroup))
  median_df <- unfiltered_data %>% filter(statistic == "median")
  freq_df <- unfiltered_data %>% filter(statistic == "percent")
  filtered_median_df <- median_df %>% select(population, reagent, group = VaccineGroup, sampleID = SampleID, stimulation = Stimulation, time = Day, median = value)
  filtered_freq_df <- freq_df %>% select(population, reagent, group = VaccineGroup, sampleID = SampleID, stimulation = Stimulation, time = Day, frequency = value)
  mediansave <- paste0(stringpath,"/volc_CIVfiltered_median_df.csv")
  freqsave <- paste0(stringpath,"/volc_CIVfiltered_freq_df.csv")
  write.csv(filtered_median_df, mediansave,row.names = FALSE)
  write.csv(filtered_freq_df, freqsave,row.names = FALSE) 
}

#this is where the filtered and preprocessed data will go
dir.create("preprocessed data for volc plot")
current_directory <- getwd()
filtered_input_data <- file.path(current_directory, "preprocessed data for volc plot")

# execute custom function using variables defined above
filter_my_data(cell_engine_data,filtered_input_data)

folder_name <- "preprocessed data for volc plot"
median_file_name <- "volc_CIVfiltered_median_df.csv"
funct_path <- file.path(current_directory, folder_name, median_file_name)

freq_file_name <- "volc_CIVfiltered_freq_df.csv"
freq_path <- file.path(current_directory, folder_name, freq_file_name)

out_path <- file.path(current_directory, folder_name)

### Features
MyData = read.csv(funct_path) 
final_data = c()
fil = setdiff(colnames(MyData), c("stimulation", "median"))
for(timepoint in unique(MyData$time)){
  
  MyData_timepoint = MyData %>% 
    filter(time == timepoint) 
  
  stims = unique(MyData_timepoint$stimulation) 
  
  MyData_timepoint = MyData_timepoint %>% 
    mutate(feature = asinh(median/5)) %>%  #arc sign transformation (if you change this, remember to also change for PB)
    #mutate(feature = (median)) %>% 
    select(-one_of("median")) %>%
    pivot_wider(names_from = stimulation, values_from = feature)
  
  for (stim in stims){
    if (stim!="U"){
      MyData_timepoint[stim] = MyData_timepoint[stim] - MyData_timepoint["U"] #(if you change this, remember to also change for PB)
    }
  }
  
  MyData_fin = MyData_timepoint %>% 
    pivot_longer(-all_of(fil),
                 names_to="stimulation",values_to="feature")
  
  final_data = rbind(final_data, MyData_fin)
  
}
freq = read.csv(freq_path) %>% rename(feature = frequency) %>% mutate(reagent = "frequency")
final_data = rbind(final_data, freq)
final_data$tissue_type <- selected_samp_type
final_data <- final_data %>% filter(time == "d7")
csv_file_path <- file.path(out_path, paste(selected_samp_type, "volc_preprocessed.csv", sep = "_"))
write.csv(final_data, csv_file_path)


#Dave asked to make a file that combines the PB and LN preprocessed files
preproc_folder_name <- "preprocessed data for volc plot"
LN_file_path <- file.path(getwd(), preproc_folder_name, paste("LN", "_volc_preprocessed.csv", sep = ""))
PB_file_path <- file.path(getwd(), preproc_folder_name, paste("PB", "_volc_preprocessed.csv", sep = ""))
LN_preprocessed_data <- read.csv(LN_file_path)
PB_preprocessed_data <- read.csv(PB_file_path)
comb_preprocessed_data <- rbind(PB_preprocessed_data,LN_preprocessed_data)
comb_preprocessed_data <- comb_preprocessed_data[, -c(1)] #delete unnec col
comb_csv_file_path <- file.path(out_path, "volc_combined_preprocessed.csv")
write.csv(comb_preprocessed_data, comb_csv_file_path)



#========================= Fig 4a Volcano plots ======================
#-------------------- PB Overall ANOVA and pairwise comp stat tests--------------------

selected_samp_type <- "PB"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data for volc plot"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_volc_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/signaling_markers")
dir.create("output/signaling_markers/volcano_plot")
volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#changing NA's to 0 and Inf to the max value
all_mark_data$feature[is.na(all_mark_data$feature)] <- 0
max_value <- max(all_mark_data$feature[is.finite(all_mark_data$feature)], na.rm = TRUE)
all_mark_data$feature[is.infinite(all_mark_data$feature)] <- max_value

#make empty lists to build into
unique_stims <- unique(all_mark_data$stimulation)
unique_cell_pops <- unique(all_mark_data$population)
unique_markers <- unique(all_mark_data$reagent)

#list to store overall ANOVA results
overall_result_list <- list()
#list to store pairwise comparison results
resulting_list <- list()


#For running the t-test and creating a list of differential expression
all_mark_data_ttest <- all_mark_data %>% pivot_wider(id_cols = c("population", "reagent", "group", "sampleID") , names_from = stimulation, values_from = feature)
all_mark_data_ttest <- all_mark_data_ttest %>% rename(cell_pop = population, marker = reagent)

#Reverse the subtraction for the P and R stimulation to run the t-test values

all_mark_data_ttest$P <- all_mark_data_ttest$P + all_mark_data_ttest$U
all_mark_data_ttest$R <- all_mark_data_ttest$R + all_mark_data_ttest$U

ttestvalues_set <- data.frame(cell_pop = factor(), marker = factor(), 
                              P = numeric(), R = numeric())

for (c in unique_cell_pops) {
  for (m in unique_markers) {
    
    subset_ttest <- all_mark_data_ttest %>% filter(cell_pop == c, marker == m)
    #calculate t-test for any given set of data
    
    
    P <- t.test(subset_ttest$P, subset_ttest$U)
    R <- t.test(subset_ttest$R, subset_ttest$U)
    
    dynamic_df <- data.frame(cell_pop = c, marker = m, 
                             P = P$p.value, 
                             R = R$p.value )
    
    ttestvalues_set <- rbind(dynamic_df, ttestvalues_set)
    
  }
}

#Correct for multiple hypotheses
ttestvalues_set <- ttestvalues_set %>% pivot_longer(cols= c("P","R"), names_to = "stim_cond", values_to = 'ttest_stim_pval_unadj')

ttestvalues_set$ttest_stim_pval_BH <- p.adjust(ttestvalues_set$ttest_stim_pval_unadj, method = "BH")

ttestvalues_set_PB <- ttestvalues_set

#For the mean feature for any given set of pairs
all_mark_data <- all_mark_data %>% filter(stimulation != "U")

unique_stims <- unique_stims[unique_stims != "U"]

#calculates mean feature and then adds it to an overall list
for (s in unique_stims) {
    for (c in unique_cell_pops) {
      for (m in unique_markers) {
        subset_data <- all_mark_data %>% filter(stimulation == s, population == c, reagent == m)

        #calculate mean of 'feature' for the current subset
        mean_feature <- mean(subset_data$feature, na.rm = TRUE)  # all vaccine groups pooled

        #include mean 'feature' value in the overall_df
        overall_df <- data.frame(cell_pop = c, stim_cond = s, marker = m, mean_feature = mean_feature)
        overall_result_list[[length(overall_result_list) + 1]] <- overall_df
      }
    }
  }


mark_mean_feature <- bind_rows(overall_result_list)

#Make naming conventions the same between mark_overall_pvals and ttestvalues_set_PB and then merge

mark_overall_info <- merge(mark_mean_feature, ttestvalues_set_PB, by = c("cell_pop", "stim_cond", "marker"))



volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")
file_name <- paste("mark_overall_info_", selected_samp_type, ".csv",sep="")
full_file_path <- file.path(volc_file_path, file_name)
mark_overall_info$tissue_type <- selected_samp_type
write.csv(mark_overall_info, full_file_path, row.names = FALSE)

#-------------------- LN Overall ANOVA and pairwise comp stat tests--------------------
#same code chunk as above just with LN
selected_samp_type <- "LN"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data for volc plot"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_volc_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/signaling_markers")
dir.create("output/signaling_markers/volcano_plot")
volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#changing NA's to 0 and Inf to the max value
all_mark_data$feature[is.na(all_mark_data$feature)] <- 0
max_value <- max(all_mark_data$feature[is.finite(all_mark_data$feature)], na.rm = TRUE)
all_mark_data$feature[is.infinite(all_mark_data$feature)] <- max_value

#vectors of stim, markers, and cell-pops
unique_stims <- unique(all_mark_data$stimulation)
unique_cell_pops <- unique(all_mark_data$population)
unique_markers <- unique(all_mark_data$reagent)

#list to store overall ANOVA results
overall_result_list <- list()
#list to store pairwise comparison results
resulting_list <- list()


#For running the t-test and creating a list of differential expression
all_mark_data_ttest <- all_mark_data %>% pivot_wider(id_cols = c("population", "reagent", "group", "sampleID") , names_from = stimulation, values_from = feature)
all_mark_data_ttest <- all_mark_data_ttest %>% rename(cell_pop = population, marker = reagent)

#Reverse the subtraction for the P and R stimulation to run the t-test values

all_mark_data_ttest$P <- all_mark_data_ttest$P + all_mark_data_ttest$U
all_mark_data_ttest$R <- all_mark_data_ttest$R + all_mark_data_ttest$U

ttestvalues_set <- data.frame(cell_pop = factor(), marker = factor(), 
                              P = numeric(), R = numeric())

for (c in unique_cell_pops) {
  for (m in unique_markers) {
    
    subset_ttest <- all_mark_data_ttest %>% filter(cell_pop == c, marker == m)
    #calculate t-test for any given set of data

    P <- t.test(subset_ttest$P, subset_ttest$U)
    R <- t.test(subset_ttest$R, subset_ttest$U)
    
    dynamic_df <- data.frame(cell_pop = c, marker = m, 
                             P = P$p.value, 
                             R = R$p.value )
    
    ttestvalues_set <- rbind(dynamic_df, ttestvalues_set)
    
  }
}

#Correct for multiple hypotheses
ttestvalues_set <- ttestvalues_set %>% pivot_longer(cols= c("P","R"), names_to = "stim_cond", values_to = 'ttest_stim_pval_unadj')

ttestvalues_set$ttest_stim_pval_BH <- p.adjust(ttestvalues_set$ttest_stim_pval_unadj, method = "BH")

ttestvalues_set_LN <- ttestvalues_set

#For the mean feature for any given set of pairs
all_mark_data <- all_mark_data %>% filter(stimulation != "U")

unique_stims <- unique_stims[unique_stims != "U"]

#calculates mean feature and then adds it to an overall list
for (s in unique_stims) {
  for (c in unique_cell_pops) {
    for (m in unique_markers) {
      subset_data <- all_mark_data %>% filter(stimulation == s, population == c, reagent == m)
      
      #calculate mean of 'feature' for the current subset
      mean_feature <- mean(subset_data$feature, na.rm = TRUE)  # all vaccine groups pooled
      
      #include mean 'feature' value in the overall_df
      overall_df <- data.frame(cell_pop = c, stim_cond = s, marker = m, mean_feature = mean_feature)
      overall_result_list[[length(overall_result_list) + 1]] <- overall_df
    }
  }
}


mark_mean_feature <- bind_rows(overall_result_list)

#Make naming conventions the same between mark_overall_pvals and ttestvalues_set_PB and then merge

mark_overall_info <- merge(mark_mean_feature, ttestvalues_set_LN, by = c("cell_pop", "stim_cond", "marker"))



volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")
file_name <- paste("mark_overall_info_", selected_samp_type, ".csv",sep="")
full_file_path <- file.path(volc_file_path, file_name)
mark_overall_info$tissue_type <- selected_samp_type
write.csv(mark_overall_info, full_file_path, row.names = FALSE)

#-------------------- Overall ANOVA and pairwise comp stat tests combine both tissues--------------------

volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")
file_path_LN_mark <- file.path(volc_file_path, "mark_overall_info_LN.csv")
file_path_PB_mark <- file.path(volc_file_path, "mark_overall_info_PB.csv")
sigdata_LN <- read.csv(file_path_LN_mark)
sigdata_PB <- read.csv(file_path_PB_mark)
sigcombined_data <- rbind(sigdata_LN, sigdata_PB)
sigcombined_file_path <- file.path(volc_file_path, "mark_overall_info_both.csv")
write.csv(sigcombined_data, sigcombined_file_path, row.names = FALSE)

names(sigdata_LN)
names(sigdata_PB)
#now that the data prepping for the volcano plots is complete, lets make the plot

#-----------Making volcano plots-------

#load in the overall pval data (both tissues)
volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")
volc_data_path <- paste0(volc_file_path, "/mark_overall_info_both.csv", sep = "")
volc_data <- read.csv(volc_data_path)

# create a new column that combines 'cell_pop' and 'marker'
volc_data$label <- paste(volc_data$cell_pop, volc_data$marker, sep = "_")

# add a new column to classify points based on p-value and mean_feature
volc_data$class <- with(volc_data, ifelse(ttest_stim_pval_unadj < 0.05 & #adjust p value significance threshold here
                                            abs(mean_feature) > 0.1, #adjust mean feature significance threshold here
                                          ifelse(mean_feature > 0, "Upregulated", "Downregulated"), 
                                          "Not significant"))

unique_tissue_types <- unique(volc_data$tissue_type)
unique_stim_conds <- unique(volc_data$stim_cond)

#Note: if you are attempting to normalize with division instead of subtraction (Stim/unstim), you need to add the following line of code into the plot
#   subset_data$mean_feature_adj <- subset_data$mean_feature + 0.000001  # Add a small constant to mean_feature to avoid log(0), adjust this constant as needed

#Add in correct labels for ttest from tissue type

#plot

for (tissue in unique_tissue_types) {
  for (stim in unique_stim_conds) {
    
    #stim = "R"
    #tissue = "PB"
    
    subset_data <- subset(volc_data, tissue_type == tissue & stim_cond == stim)
    
    # adjust y-values exceeding the limit
    subset_data$adjusted_p_val <- with(subset_data, ifelse(-log10(ttest_stim_pval_unadj) > 20, 20, -log10(ttest_stim_pval_unadj))) # y axis max can be adjusted here (also adjust it within the gg plot below)
    
    # Adjust x-values exceeding the limits
    subset_data$adjusted_mean_feature <- with(subset_data, pmin(pmax(mean_feature, -2.25), 2.25)) # x axis min/max can be adjusted here  (also adjust it within the gg plot below)
    
    p <- ggplot(subset_data, aes(x = adjusted_mean_feature, y = adjusted_p_val, color = class)) +
      geom_jitter(size=1.8, width = 0.05, height = 0) +  # added jitter
      geom_text_repel(aes(label = ifelse(class != "Not significant", as.character(label), "")),
                      size = 1.5, 
                      box.padding = 0.35, 
                      point.padding = 0.5,
                      max.overlaps = Inf,
                      show.legend = FALSE,
                      segment.size = 0.2,  
                      segment.color = "grey70", 
                      segment.alpha = 0.5) +  
      scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "grey")) +
      ggtitle(paste("Volcano Plot:", tissue, stim)) +
      xlab("Mean Feature (arcsign transformed)") +  
      ylab("-log10(p-value)") +
      theme_classic() +  
      theme(legend.title = element_blank()) +
      scale_y_continuous(limits = c(min(subset_data$adjusted_p_val), max(subset_data$adjusted_p_val))) +  # y axis max can be adjusted here (also adjust it within subset data above)
      scale_x_continuous(limits = c(-2.25, 2.25))  # x axis min/max can be adjusted here (also adjust it within subset data above)
    
    plot(p)
    print(p)
    filename <- file.path(volc_file_path, paste0("volcano_", tissue, "_", stim, ".svg"))
    ggsave(filename, plot = p, width = 8, height = 7)
  }
}

