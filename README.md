The scripts in this repo are designed to be user-friendly, even for those with limited coding experience. Each script includes instructions and comments at the top to guide you through any adjustments that need to be made before running.

Importantly, first follow the instructions in Script #1 to create a folder where this analysis will take place and download the raw CellEngine data into it.

Script Overview:

Script #1: Arcsinh_JH_sasha_plug_and_play.R  (for preprocessing data)

Description: Takes the data from Cell Engine and performs data transformations. Normalizes stimulation values by subtracting Stim - Unstim. Unstim values remain untouched. The output is a “preprocessed” csv file for each tissue ready for further analysis.

Script #2: CIV_abund_p_vals_V5.R (for Fig 2 supp)

Description: This script generates cell abundance CSV files with p-values, filtered for UNSTIM samples only. It performs an ANOVA to calculate an overall p-value for comparing the three vaccine groups and conducts pairwise comparisons using a t-test after the ANOVA, with BH corrections applied globally within each tissue type.

Script #3: CIV_Fig2A_2C.R (for Fig 2A)

Description: Uses the preprocessed data from Arcsinh_JH_sasha_plug_and_play.R to create heatmaps of populations by frequency (value = log2(x/control)).

Script #4: CIV_abund_plots_V4.R  (for Fig 2B and 2D)

Description: This script produces cell abundance box plots for both PB and LN. Filtered for unstim and loops through each vaccine treatment.

Script #5: CIV_marker_p_vals_V2.R (for Fig 2 supp)

Description: This script produces two different dataframes of p-values. Similar to script #2, it performs an ANOVA to get an overall pvalue for comparing the 3 vaccine groups for each stimulation marker-cell population pair. Afterwards, it performs t-test between different vaccine groups. Values normalized to unstim (by subtracting Stim - Unstim).

Script #6: CIV_alluvial_plots_V5.R (for Fig 3A)

Description: Generates alluvial plots to visualize the distribution of cell populations across tissue types, vaccine groups, and stimulations. Filters for stimulation "P," calculates mean frequencies, and corrects values for nested cell populations (e.g., B cells, CD4T, CD8T). The corrected data is then used to produce an alluvial plot, where flow represents the mean frequencies, scaled and color-coded by tissue type, vaccine group, and cell population.

Script #7: CIV_comp_bw_tissues_freq_V6 (for Fig 3B)

Description: Outputs CSV file with mean frequencies and adjusted p-values for cell populations in PB and LN tissues, calculated using ANOVA tests and BH correction. Also generates plots (log scaled as well as unscaled) that visualize these frequencies, highlighting statistically significant populations with shaded polygons and labels (corrected p < 0.05). Calculates correlation values between LN and PB for populations given a vaccine grouping. Any grouping has a triangle if a t-test in either LN or PB is significant.

Script #8: CIV_marker_volcano_plots_v6_CM.R (for Fig 4A)

Description: First section of code is the same as Arcsinh_JH_sasha_plug_and_play.R for PB and then for LN. Subtracts (stim - unstim) to normalize and performs inverse hyperbolic sine transformation. Middle section of code performs statistical tests. Performs a t-test for the stim versus unstim and combining them into a full list. Final section of code makes volcano plots using the p-value stat from this t-test. Outputs 4 volcano plots, 2 stim cond x 2 tissues. Additionally, outputs CSV with ttest results and magnitude change between stim versus unstim.

Script #9: CIV_bubble_plots_markerV4.R (for Fig 4B)

Description: The first section of code takes the statistical tests from Script #8 and the processed data from Script #1. It finds the top 20% of cell type-phospho marker pairs by the lowest ranking sums of both significance and magnitude for each tissue type-stimulation pair. The code then finds the list of these cell type-phospho marker pairs that are inclusive to all tissue type-stimulation pairs. The code then runs an ANOVA on all of these cell type-phospho marker pairs by ANOVA and finds the top 3 hits for each tissue type-stimulation pair. The top 3 hits for each tissue type-stimulation pair are all plotted in a heat map for all tissue type-stimulation pairs, making each tile a distinct cellular type-phospho marker for a given tissue type, stimulation condition, and vaccine group. A Z-score is calculated as the number of standard deviations a vaccine group mean value is away from the mean for a cell type-phospho marker regardless of vaccine group and plotted as the color of individual tiles in the heatmap. Size of individual tiles correlates to the Benjamini-Hochberg procedure corrected significance values for the ANOVA against vaccine groups.

Script #10: CIV_correl_plots_freq_V3.R  (for Fig 5A)

Description: Uses raw data from Cell Engine and frequency p-values (specifically the ANOVA p-value) from Table S2. Script merges these two based on sampleID and selects the top 10 populations that have the best p-value hits. Then, it reads in metadata and makes plots for the following: pathology score vs freq of significant pop, cov rna scores vs freq of significant pop, radiography scores from day0 to day7 vs freq of significant pop, and radiography scores from day7 vs freq of significant pop. P values reflect the p value of the correlation, adjusted for mult hypothesis with BH correction. Creates csv output with all R and p values.

Script #11: CIV_correl_plots_marker_V4.R  (for Fig 5B)

Description: Uses top three statistical hits by ANOVA for cell pop-phospho marker pair for every tissue-stim pair (as previously determined in script #11), corrects for mult hypothesis testing for 10 hypotheses by BH correction - 3 hits by 4 tissue-stim pairs minus 2 redundancies = 10 tests. Correlation test is conducted between given cell pop-marker pair’s marker intensities and given metadata scores, correlation coefficient R and p value of correlation is calculated and adjusted for mult hypothesis with BH correction. The correlation is plotted. CSV’s are created with R value and adjusted p value. Note: Variables at the top of script need to be adjusted for which tissue and which stim to run the script for.

Script #12: CIV_dataframe_merger.R (for table)

Description: This takes the csv files from script #2 (CIV_abund_plots_V5.R) and the csv files from script #4 (CIV_marker_p_vals_V2.R) and combines them into two dataframes. One for frequency p-values and one for marker p-values. The combined dataframes are written as CSVs in the main file.

Script #13: CIV_Heatmap_Correlations.R (for Fig 5C,D)

Description: This uses the previously made files from the correlation plots to create heatmaps. Heatmaps show the R-value for the correlation and are shaded by significance (p<0.05 are cross-hatched). It prints two PDFs. One for the frequency data and one for the marker data.

Variable naming conventions:
cell_pop: cell population
group: vaccine group (Control for Mock group, Protein for Protein Group, mRNA for mRNA group)
stim_cond: stimulation condition (P for PMA/I Stimulation, R for R848 Stimulation, U for unstimulated)
marker: signaling marker
pop_mark: combination of given cell population and signaling marker
time: d7 (all animals in this study were sacrificed and necropsied at 7 days post-challenge)
tissue_type: tissue (PB for peripheral blood mononuclear cells and LN for dissociated lymph node)
mean_feature/frequency: value for a given mean_feature (signal intensity) or population %
group_pairs: Groups compared during a t-test. Will include all three vaccination groups if an ANOVA
