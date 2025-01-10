#Sasha's code for CIV analysis
#Oct2023

library(tidyverse)
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)

# This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis will take place
setwd("Replace me")

# download raw cell engine data to that folder and paste the file path to it
cell_engine_data <- read.csv("Replace me")


#------------------- Do not edit past this point ---------------------------

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
  mediansave <- paste0(stringpath,"/CIVfiltered_median_df.csv")
  freqsave <- paste0(stringpath,"/CIVfiltered_freq_df.csv")
  write.csv(filtered_median_df, mediansave,row.names = FALSE)
  write.csv(filtered_freq_df, freqsave,row.names = FALSE) 
  }

#this is where the filtered and preprocessed data will go
dir.create("preprocessed data")
current_directory <- getwd()
filtered_input_data <- file.path(current_directory, "preprocessed data")

# execute custom function using variables defined above
filter_my_data(cell_engine_data,filtered_input_data)

folder_name <- "preprocessed data"
median_file_name <- "CIVfiltered_median_df.csv"
funct_path <- file.path(current_directory, folder_name, median_file_name)

freq_file_name <- "CIVfiltered_freq_df.csv"
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
    mutate(feature = asinh(median/5)) %>% 
    select(-one_of("median")) %>%
    pivot_wider(names_from = stimulation, values_from = feature)
  
  for (stim in stims){
    if (stim!="U"){
      MyData_timepoint[stim] = MyData_timepoint[stim] - MyData_timepoint["U"]
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
csv_file_path <- file.path(out_path, paste(selected_samp_type, "preprocessed.csv", sep = "_"))
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
  mediansave <- paste0(stringpath,"/CIVfiltered_median_df (not yet adjusted).csv") #as in, not yet stim-unstim or arcsign transformed
  freqsave <- paste0(stringpath,"/CIVfiltered_freq_df.csv")
  write.csv(filtered_median_df, mediansave,row.names = FALSE)
  write.csv(filtered_freq_df, freqsave,row.names = FALSE) 
}

#this is where the filtered and preprocessed data will go
dir.create("preprocessed data")
current_directory <- getwd()
filtered_input_data <- file.path(current_directory, "preprocessed data")

# execute custom function using variables defined above
filter_my_data(cell_engine_data,filtered_input_data)

folder_name <- "preprocessed data"
median_file_name <- "CIVfiltered_median_df (not yet adjusted).csv"
funct_path <- file.path(current_directory, folder_name, median_file_name)

freq_file_name <- "CIVfiltered_freq_df.csv"
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
    mutate(feature = asinh(median/5)) %>% 
    select(-one_of("median")) %>%
    pivot_wider(names_from = stimulation, values_from = feature)
  
  for (stim in stims){
    if (stim!="U"){
      MyData_timepoint[stim] = MyData_timepoint[stim] - MyData_timepoint["U"]
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
csv_file_path <- file.path(out_path, paste(selected_samp_type, "preprocessed.csv", sep = "_"))
write.csv(final_data, csv_file_path)

#Dave asked to make a file that combines the PB and LN preprocessed files
preproc_folder_name <- "preprocessed data"
LN_file_path <- file.path(getwd(), preproc_folder_name, paste("LN", "_preprocessed.csv", sep = ""))
PB_file_path <- file.path(getwd(), preproc_folder_name, paste("PB", "_preprocessed.csv", sep = ""))
LN_preprocessed_data <- read.csv(LN_file_path)
PB_preprocessed_data <- read.csv(PB_file_path)
comb_preprocessed_data <- rbind(PB_preprocessed_data,LN_preprocessed_data)
comb_preprocessed_data <- comb_preprocessed_data[, -c(1)] #delete unnec col
comb_csv_file_path <- file.path(out_path, "combined_preprocessed.csv")
write.csv(comb_preprocessed_data, comb_csv_file_path)

