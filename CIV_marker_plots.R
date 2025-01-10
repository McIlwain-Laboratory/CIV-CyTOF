#Sasha's code for CIV analysis
#Oct2023
#After you run marker p-value script, run this for producing sig marker graphs
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

#PB sample type
selected_samp_type <- "PB"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/signaling_markers")
mark_file_path <- file.path(getwd(), "output/signaling_markers")

#-------------------- Task 5: Signaling Marker plots for PB--------------------

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#make empty lists to build into
mark_unique_timepoints <- unique(all_mark_data$time)
mark_unqiue_stims <- unique(all_mark_data$stimulation)
mark_unique_markers <- unique(all_mark_data$reagent)


#nested loop for making plots showing all populations x Days x Treatments
for (t in mark_unique_timepoints) {
  for (s in mark_unqiue_stims) {
    for (m in mark_unique_markers) {
    subset_data <- all_mark_data %>% filter(time == t, stimulation == s, reagent == m)
    graphtitle <- paste(m," ",selected_samp_type," ",t," ",s,sep="")
    ggplot(subset_data) +
      aes(x = feature, y = population, fill = group) +
      geom_boxplot() +
      scale_fill_hue(direction = 1) +
      labs(x = "Signaling Channel Level", y = "Cell Population", title = graphtitle) +
      theme_minimal()
    
    ggfilename<- paste("marker_levels_",m,"_",selected_samp_type,"_",t,"_",s,".png",sep="")
    ggsave(ggfilename,path = mark_file_path,bg="white")
  }
}}

#LN sample type
selected_samp_type <- "LN"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/signaling_markers")
mark_file_path <- file.path(getwd(), "output/signaling_markers")

#-------------------- Task 5: Signaling Marker plots for LN --------------------

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#make empty lists to build into
mark_unique_timepoints <- unique(all_mark_data$time)
mark_unqiue_stims <- unique(all_mark_data$stimulation)
mark_unique_markers <- unique(all_mark_data$reagent)


#nested loop for making plots showing all populations x Days x Treatments
for (t in mark_unique_timepoints) {
  for (s in mark_unqiue_stims) {
    for (m in mark_unique_markers) {
      subset_data <- all_mark_data %>% filter(time == t, stimulation == s, reagent == m)
      graphtitle <- paste(m," ",selected_samp_type," ",t," ",s,sep="")
      ggplot(subset_data) +
        aes(x = feature, y = population, fill = group) +
        geom_boxplot() +
        scale_fill_hue(direction = 1) +
        labs(x = "Signaling Channel Level", y = "Cell Population", title = graphtitle) +
        theme_minimal()
      
      ggfilename<- paste("marker_levels_",m,"_",selected_samp_type,"_",t,"_",s,".png",sep="")
      ggsave(ggfilename,path = mark_file_path,bg="white")
    }
  }}

