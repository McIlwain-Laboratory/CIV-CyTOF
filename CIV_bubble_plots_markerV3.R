#Cameron's Code for CIV Analysis
#Aug 2024
#Run this to create a bubble plot of the top ten hits for t-test between stimulated
#and unstimulated.
#Each bubble will include the p-value of the ANOVA by size
#And a z-score coloring of stim/unstim
#Note Oct2024 - This code is rewritten so that the top 3 hits per ANOVA by 
#group for each tissue type and stimulation are 
#explored for the all of the tissue type and stimulation pairs. This will help
#us understand the underlying biology that may be different between stimulations
#and tissue type for a given vaccine and cellpop/marker pair.

#Note: p-value Overall is the ANOVA p-value

library(tidyverse)
library(ggdendro)
library(ggplot2)
library(grid)
library(gridExtra)
library(svglite)

#Download ANOVA p-values  
setwd("Replace me")
volc_file_path <- file.path(getwd(), "output/signaling_markers/volcano_plot")
volc_data_path <- paste0(volc_file_path, "/mark_overall_info_both.csv", sep = "")
volc_data <- read.csv(volc_data_path)

#Download full combined data set and clean
current_directory <- getwd()

folder_name <- "preprocessed data for volc plot"
out_path <- file.path(current_directory, folder_name)
comb_csv_file_path <- file.path(out_path, "volc_combined_preprocessed.csv")

preprocessed_data <- read.csv(comb_csv_file_path)

all_mark_data <- preprocessed_data %>% filter(reagent != "frequency")

#changing NA's to 0 and Inf to the max value
all_mark_data$feature[is.na(all_mark_data$feature)] <- 0
max_value <- max(all_mark_data$feature[is.finite(all_mark_data$feature)], na.rm = TRUE)
all_mark_data$feature[is.infinite(all_mark_data$feature)] <- max_value

all_mark_data <- all_mark_data %>% rename(cell_pop = population, marker = reagent, stim_cond = stimulation)

#------------------------------------------------------------------------
#Filter significance data into four groupings based on tissue type and stimulation
tissue_types = unique(all_mark_data$tissue_type)
stim_types = c("P", "R")

#Creates a dataframe that will take in all of the names of the top 20% cell_type 
#marker pairs and split them into the tissue and stimulation
#This will be further divided into a dataframe that is only the unique
#cell_type marker pairs for all tissue/stimulation.
#Another pair of data frames will be only the unique cell_type marker pairs for
#A given stimulation regardless of tissue type (one for PMA and the other for R848)
top20percent_pop_mark_total <- data.frame(cell_pop = character(), marker = character(),
                                          tissue_type = character(), stim_cond = character())

for(stim in stim_types){
  for(tissue in tissue_types){

  
    Top_Hits_SigList <- volc_data %>% 
      filter(tissue_type == tissue, stim_cond == stim) %>%
      mutate(abs_mean_feature = abs(mean_feature)) %>%
      mutate(PValue_Rank = rank(ttest_stim_pval_BH, ties.method = "average")) %>%
      mutate(Mgtde = rank(-abs_mean_feature, ties.method = "average")) %>%
      mutate(combRank = PValue_Rank + Mgtde) %>%
      arrange((combRank)) %>%
      slice(1:39) %>% 
      select(cell_pop, marker, tissue_type, stim_cond)
    #Currently doing the top20% by ranking (mean feature + significance value) 195 cell-pop/marker pairs for 39 samples
    
    top20percent_pop_mark_total <- rbind(top20percent_pop_mark_total, Top_Hits_SigList)
    
  }
}

#Unique total pairs 
unique_all_stim_pairs <- unique(top20percent_pop_mark_total[, c("cell_pop", "marker")])

#Unique total pairs split by stimulation
for (stim in stim_types){
  top20_stim_df <- top20percent_pop_mark_total %>% filter (stim_cond == stim)
  
  top20_stim_pop_mark <- unique(top20_stim_df[, c("cell_pop", "marker")])
  
  #dynamically create new dataframe with name
  data_frame_name <- paste0("unique_", stim, "_stim_pairs")
  assign(data_frame_name, top20_stim_pop_mark)
}

#We now have the top 20% stim_pairs aggregated for any stim, tissue pair we will create a dataframe
#that has the p-value ANOVA between samples and a Z-score. Once we have this dataframe
#we can reorganize it with regards to both the heatmap and the dendrogram.

#I can use merge(all_mark_data, unique_set) to get only the data that we need 
#from each unique set.

#For the set that includes both stimulation conditions and tissue types
unique_all_stim_pairs_all_data <- merge(all_mark_data, unique_all_stim_pairs) %>% 
  filter(stim_cond != "U")

#Create dataframe to get all of the ANOVA values from and then do three loops to sum

unique_all_stim_pairs_stats <- data.frame(cell_pop = character(), marker = character(),
                                          tissue_type = character(), stim_cond = character(),
                                          anova_marker_pval_unadj = numeric())

top12dataframe <- data.frame()

for(stim in stim_types){
  for(tissue in tissue_types){
    
    subtop12dataframe <- data.frame()
    
    for (i in 1:nrow(unique_all_stim_pairs)){
      #create filtering
      cell_population <- unique_all_stim_pairs$cell_pop[i]
      stimulation_marker <- unique_all_stim_pairs$marker[i]
      
      #Filter subset df to do the ANOVA with
      subset_df <- unique_all_stim_pairs_all_data %>% filter(stim_cond == stim) %>%
        filter(tissue_type == tissue) %>% 
        filter(cell_pop == cell_population & marker == stimulation_marker)
      
      #Run the ANOVA 
      overall_anova_result <- aov(feature ~ group, data = subset_df)
      overall_p_value <- summary(overall_anova_result)[[1]]$`Pr(>F)`[1]
      
      #Put all of the values into a dataframe and combine with unique_all_stim_pairs_stats df
      anova_df <- data.frame(cell_pop = cell_population, marker = stimulation_marker,
                             stim_cond = stim, tissue_type = tissue, anova_marker_pval_unadj = overall_p_value)
      
      #For overall dataframe
      unique_all_stim_pairs_stats <- rbind(unique_all_stim_pairs_stats, anova_df)
      
      #For finding top 3 per
      subtop12dataframe <- rbind(subtop12dataframe, anova_df)
    }
    #Create 
    subtop12dataframe <- subtop12dataframe %>% arrange(subtop12dataframe$anova_marker_pval_unadj) %>%
      slice(1:3) %>% select(cell_pop, marker)
    top12dataframe <- rbind(top12dataframe, subtop12dataframe)
  }
}

#Create Z-score for all of the pairs
unique_all_stim_pairs_zscores <- unique_all_stim_pairs_all_data %>%
  group_by(cell_pop, marker, group, tissue_type, stim_cond) %>% 
  summarise(mean_feature_byGroup = mean(feature), .groups = 'drop') %>%
  group_by(cell_pop, marker, stim_cond, tissue_type) %>%
  mutate(mean_z = c(scale(mean_feature_byGroup))) %>%
  ungroup()

#Create a master dataframe
unique_all_stim_pairs_master <- merge(unique_all_stim_pairs_stats, 
                                      unique_all_stim_pairs_zscores, 
                                      all.y = TRUE)

#Get only unique pairs from top 12 dataframe
top12dataframe <- unique(top12dataframe[, c("cell_pop", "marker")])

#organize into dataframes dynamically so that we get the top 3 by ANOVA for each

top12_unique_all_stim_pairs_master <- merge(unique_all_stim_pairs_master, top12dataframe)

rownames_for_dendo <- unique(top12_unique_all_stim_pairs_master[,c("cell_pop", "marker")])

#Create dendogram for clustering
Top_Hits_DataPoints_Dendo <- top12_unique_all_stim_pairs_master %>% 
  pivot_wider(id_cols = c("cell_pop", "marker"),
              names_from = c(group, tissue_type, stim_cond), 
              values_from = mean_z,names_sep = "_") %>% 
  select(-cell_pop, -marker) %>%
  as.matrix()

rownames(Top_Hits_DataPoints_Dendo) <- paste(rownames_for_dendo$cell_pop, rownames_for_dendo$marker, sep = "_")

Top_Hits_DataPoints_Dendogram <- as.dendrogram(hclust(d = dist(x = Top_Hits_DataPoints_Dendo)))

#dendro_plot <- ggdendrogram(data = Top_Hits_DataPoints_Dendogram, rotate = TRUE)
#print(dendro_plot) #Just to check that the function is working

##Extract order from the dendrogram
name_order <- dendro_data(Top_Hits_DataPoints_Dendogram, type = "rectangle")
label_order <- name_order$labels
label_order2 <- label_order$label

#Make combined name for both (cell_pop, marker) + (tissue_type, stim_cond, group)
top12_unique_all_stim_pairs_master$full_name <- paste(top12_unique_all_stim_pairs_master$cell_pop, top12_unique_all_stim_pairs_master$marker, sep = "_")
top12_unique_all_stim_pairs_master$metaname <- paste(top12_unique_all_stim_pairs_master$stim_cond, 
                                                     top12_unique_all_stim_pairs_master$tissue_type,
                                                     top12_unique_all_stim_pairs_master$group, sep = "_")


#Set levels of combined names to match dendrogram for heatmap

top12_unique_all_stim_pairs_master$full_name <- as.factor(top12_unique_all_stim_pairs_master$full_name)

top12_unique_all_stim_pairs_master$full_name <- factor(top12_unique_all_stim_pairs_master$full_name,
                                                  levels = label_order2, 
                                                  ordered = TRUE)

#metaname order
metaname_order <- c("P_LN_Control", "P_LN_mRNA", "P_LN_Protein","P_PB_Control", "P_PB_mRNA", "P_PB_Protein",
                    "R_LN_Control", "R_LN_mRNA", "R_LN_Protein","R_PB_Control", "R_PB_mRNA", "R_PB_Protein")

top12_unique_all_stim_pairs_master$metaname <- as.factor(top12_unique_all_stim_pairs_master$metaname)

top12_unique_all_stim_pairs_master$metaname <- factor(top12_unique_all_stim_pairs_master$metaname,
                                                      levels = metaname_order,
                                                      ordered = TRUE)

#Do BH correction for ANOVA

top12_unique_all_stim_pairs_master <- top12_unique_all_stim_pairs_master %>% 
  group_by(tissue_type, stim_cond) %>%
  mutate(anova_marker_pval_BH = p.adjust(anova_marker_pval_unadj, method = "BH", n = 80)) %>% 
  ungroup()

#Print heat maps - can change the size with scale_size
top12_heatmap <- ggplot(top12_unique_all_stim_pairs_master, aes(metaname, full_name, col = mean_z)) +
  geom_tile(col = "black", fill = "white") +
  geom_point(aes(size = abs(-log10(anova_marker_pval_BH))), shape = 15) +
  scale_color_gradient2(name = "Mean Z-score", mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1.5,1.5)) +
  scale_size_area(name = "-log10(FDR)", max_size = 20, limits = c(0,4)) +
  theme_bw() #+
#theme(    axis.text.y = element_blank())

print(top12_heatmap)


full_bubbleplot_path <- file.path(paste(current_directory, "/output/signaling_markers/heatmap_fig4b.svg", sep = ""))

#Save the plot with specified dimensions
ggsave(full_bubbleplot_path, plot = top12_heatmap, width = 11, height = 6.77, bg = "white") 

#Reorganize and save csv file with all of the p-value info
top12_unique_all_stim_pairs_master <- top12_unique_all_stim_pairs_master %>% relocate(anova_marker_pval_BH, .after = anova_marker_pval_unadj) %>%
  relocate(group, .after = tissue_type) %>% 
  left_join(volc_data %>% 
              select(ttest_stim_pval_unadj, ttest_stim_pval_BH, cell_pop, stim_cond, marker, tissue_type), 
            by = c("cell_pop", "stim_cond", "marker", "tissue_type")) %>%
  relocate(full_name, .after = ttest_stim_pval_BH) %>%
  relocate(metaname, .after = full_name)

csv_directory <- paste(current_directory, "/output/signaling_markers", sep = "")
write.csv(top12_unique_all_stim_pairs_master, file = file.path(csv_directory, "fig4b_pValues_updatedOct29.csv"))
