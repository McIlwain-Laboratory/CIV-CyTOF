#Cameron's Code for the CIV analysis
#Oct 2024
#This script takes all of the dataframes from the CIV_abund_p_vals_V5.R and
#CIV_marker_p_vals_V2.R and combines them into two distinct dataframes. The first
#is all of the data for frequency and the second is all the data for the median
#marker value.
#These two datasets will be added to a supplement

#This code is built to change only the working directory.

library(tidyverse)
library(dplyr)

setwd("Replace me")

#-------------------------------------------------------------------------------

#Download all of the data using the working directory and file tree format
working_directory <- getwd()

###FREQUENCY###

#Download ANOVA-frequency file
anova_freq_pvals_filename <- file.path(paste(working_directory, "/output/cell_abundance/freq_stats_overall_p_vals.csv", sep = ""))
anova_freq_pvals <- read.csv(anova_freq_pvals_filename)

#Download ttest-frequency file
ttest_freq_pvals_filename <- file.path(paste(working_directory, "/output/cell_abundance/freq_stats_pairwise_p_vals.csv", sep = ""))
ttest_freq_pvals <- read.csv(ttest_freq_pvals_filename)

#Make formatting the same between the two dataframes and pivot longer

#ANOVA
anova_freq_pvals <- anova_freq_pvals %>% rename(group_pair = groups) %>% 
  pivot_longer(cols = c(anova_freq_pval_unadj, anova_freq_pval_BH),
               names_to = "stats_test", values_to = "value" )

#ttest
ttest_freq_pvals <- ttest_freq_pvals %>% 
  pivot_longer(cols= c(ttest_groupfreq_pval_unadj, ttest_groupfreq_pval_BH), 
               names_to = "stats_test", values_to = "value")

#Bind the two dataframes together and combine
frequency_pvals <- rbind(ttest_freq_pvals, anova_freq_pvals)


write.csv(frequency_pvals, file = paste(working_directory, "/frequency_pvals.csv", sep = ""), row.names = FALSE) 

###MEDIANMARKER###

#Download ANOVA-frequency file
anova_marker_pvals_filename <- file.path(paste(working_directory, "/output/signaling_markers/sig_markers_anova_p_vals_both.csv", sep = ""))
anova_marker_pvals <- read.csv(anova_marker_pvals_filename)

#Download ttest-frequency file
ttest_marker_pvals_filename <- file.path(paste(working_directory, "/output/signaling_markers/sig_markers_ttest_p_vals_both.csv", sep = ""))
ttest_marker_pvals <- read.csv(ttest_marker_pvals_filename)

#Make formatting the same between the two dataframes and pivot longer

#ANOVA
anova_marker_pvals <- anova_marker_pvals %>% mutate(group_pair = "Control mRNA Protein") %>% 
  pivot_longer(cols = c(anova_marker_pval_unadj, anova_marker_pval_BH),
               names_to = "stats_test", values_to = "value" )

#ttest
ttest_marker_pvals <- ttest_marker_pvals %>% 
  pivot_longer(cols= c(ttest_groupmarker_pval_unadj, ttest_groupmarker_pval_BH), 
               names_to = "stats_test", values_to = "value")

#Bind the two dataframes together and combine
marker_pvals <- rbind(ttest_marker_pvals, anova_marker_pvals)


write.csv(marker_pvals, file = paste(working_directory, "/marker_pvals.csv", sep = ""), row.names = FALSE) 


