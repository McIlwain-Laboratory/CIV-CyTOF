#Sasha's code for CIV analysis
#Oct2023
#After you run cell abundance p-value script, run this for producing cell abundance graphs
#Fig 2 in manuscript
library(ggplot2)
require(tidyverse) #load packages
library(tidyselect)
library(reshape2)
library(tidyr) 
library(dplyr)
library(stats)
library(RColorBrewer)

# This code is built to be plug-and-play. Please edit the parameters below and run.

# set the working directory to the folder where this analysis is taking place
setwd("Replace me")


#------------------- Do not edit past this point ---------------------------

#start with PB
selected_samp_type <- "PB"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/cell_abundance/abund_plots")
caplots_file_path <- file.path(getwd(), "output/cell_abundance/abund_plots")

#-------------------- Task 2: Cell Abundance plots --------------------

all_freq_data <- preprocessed_data %>% filter(reagent == "frequency")
all_freq_data <- all_freq_data %>% filter(stimulation == "U", time == "d7")

#add a column that groups them by larger cell population/type (ex: CD4 and CD8 are group to "T_cell")
all_freq_data <- all_freq_data %>%
  mutate(grouped_pops = case_when(
    grepl("Bcell", population, ignore.case = TRUE) ~ "B cells and plasma cells",
    grepl("^DC$|mDC|pDC", population, ignore.case = TRUE) ~ "DCs and pDCs",
    grepl("CD4|CD8|NKT", population, ignore.case = TRUE) & !grepl("Ki67", population, ignore.case = TRUE) ~ "T cells",
    grepl("Ki67", population, ignore.case = TRUE) ~ "Activated T cells",
    grepl("cMC|intMC|ncMC", population, ignore.case = TRUE) ~ "Monocytes",
    TRUE ~ "Other"  # for types not explicitly matched
  ))

#make empty lists to build into
unique_grped_cell_pop <- unique(all_freq_data$grouped_pops)
unique_timepoints <- unique(all_freq_data$time)
unqiue_stims <- unique(all_freq_data$stimulation) #should only be one stimulation: U

#nested loop for making plots showing all populations x Days x Treatments
for (g in unique_grped_cell_pop) {
for (t in unique_timepoints) {
  for (s in unqiue_stims) {
    subset_data <- all_freq_data %>% filter(grouped_pops == g, time == t, stimulation == s)
    graphtitle <- paste(selected_samp_type," Cell Abundance for ",g," ",t," ",s,sep="")
    ggplot(subset_data) +
      aes(x = feature, y = population, fill = group) +
      geom_boxplot() +
      scale_x_continuous(limits =c(0, NA)) +
      scale_fill_manual(values = c("Control" = "#ff4040",
                                   "mRNA" = "#00ba38",
                                   "Protein" = "#619cff")) +
      labs(x = "% of CD45+CD66- Cells", y = "Cell Population", title = graphtitle) +
      theme_minimal()
    
    ggfilename<- paste("abund_",g,"_",selected_samp_type,"_",t,"_",s,".svg",sep="")
    ggsave(ggfilename,path = caplots_file_path,bg="white")
    }
}
}


#same code chunk as abpve but for LN
selected_samp_type <- "LN"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
dir.create("output/cell_abundance/abund_plots")
caplots_file_path <- file.path(getwd(), "output/cell_abundance/abund_plots")

all_freq_data <- preprocessed_data %>% filter(reagent == "frequency")
all_freq_data <- all_freq_data %>% filter(stimulation == "U") 

#add a column that groups them by larger cell population/type (ex: CD4 and CD8 are group to "T_cell")
all_freq_data <- all_freq_data %>%
  mutate(grouped_pops = case_when(
    grepl("Bcell", population, ignore.case = TRUE) ~ "B cells and plasma cells",
    grepl("^DC$|mDC|pDC", population, ignore.case = TRUE) ~ "DCs and pDCs",
    grepl("CD4|CD8|NKT", population, ignore.case = TRUE) & !grepl("Ki67", population, ignore.case = TRUE) ~ "T cells",
    grepl("Ki67", population, ignore.case = TRUE) ~ "Activated T cells",
    grepl("cMC|intMC|ncMC", population, ignore.case = TRUE) ~ "Monocytes",
    TRUE ~ "Other"  # for types not explicitly matched (Ask Dave how to categorize)
  ))

#make empty lists to build into
unique_grped_cell_pop <- unique(all_freq_data$grouped_pops)
unique_timepoints <- unique(all_freq_data$time)
unqiue_stims <- unique(all_freq_data$stimulation) #should only be one stimulation: U

#nested loop for making plots showing all populations x Days x Treatments
for (g in unique_grped_cell_pop) {
  for (t in unique_timepoints) {
    for (s in unqiue_stims) {
      subset_data <- all_freq_data %>% filter(grouped_pops == g, time == t, stimulation == s)
      graphtitle <- paste(selected_samp_type," Cell Abundance for ",g," ",t," ",s,sep="")
      ggplot(subset_data) +
        aes(x = feature, y = population, fill = group) +
        geom_boxplot() +
        scale_x_continuous(limits = c(0,NA)) +
        scale_fill_manual(values = c("Control" = "#ff4040",
                                     "mRNA" = "#00ba38",
                                     "Protein" = "#619cff")) +
        labs(x = "% of CD45+CD66- Cells", y = "Cell Population", title = graphtitle) +
        theme_minimal()
      
      ggfilename<- paste("abund_",g,"_",selected_samp_type,"_",t,"_",s,".svg",sep="")
      ggsave(ggfilename,path = caplots_file_path,bg="white")
    }
  }
}


#------------Fig 2 supplement--------
#Table of cell population frequencies for all animals and tissues

#Table cell population frequencies for all animals and tissues
#start with PB
selected_samp_type <- "PB"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
caplots_file_path <- file.path(getwd(), "output/cell_abundance")

preprocessed_data <- preprocessed_data %>% select(-X)
all_freq_data_PB <- preprocessed_data %>% filter(reagent == "frequency")
all_freq_data_PB <- all_freq_data_PB %>% filter(stimulation == "U")

selected_samp_type <- "LN"
#script below reads in the relevant preprocessed csv based on selected sample type
preproc_folder_name <- "preprocessed data"
file_path <- file.path(getwd(), preproc_folder_name, paste(selected_samp_type, "_preprocessed.csv", sep = ""))
preprocessed_data <- read.csv(file_path)
caplots_file_path <- file.path(getwd(), "output/cell_abundance")

preprocessed_data <- preprocessed_data %>% select(-X)
all_freq_data_LN <- preprocessed_data %>% filter(reagent == "frequency")
all_freq_data_LN <- all_freq_data_LN %>% filter(stimulation == "U")

all_freq_data_bothtiss <- rbind(all_freq_data_PB, all_freq_data_LN)
all_freq_data_bothtiss <- all_freq_data_bothtiss %>% select(-reagent) #get rid of reagent column
all_freq_data_bothtiss <- all_freq_data_bothtiss %>% rename(frequency = feature) #rename feature column to frequency

freqtable_file_path <- file.path(getwd(), "output/cell_abundance/freq_all_animals_and_tissues.csv")
write.csv(all_freq_data_bothtiss, file = freqtable_file_path, row.names = FALSE)

