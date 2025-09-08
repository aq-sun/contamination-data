# contamination-data

This pipeline produces a summary report and a set of site-specific reports for negative controls. 

## Processes
For the summary report, we perform the following analyses:
- A comparison of the volume of contamination in each site by averaging read counts per sample
- Grouping types of contamination by clade (bacteria, fungi, viruses, archaea, protists)
- Separating DNA and RNA contaminants into two files, and calculating Shannon’s diversity and evenness index for each set using R vegan
- A comparison of taxa richness in each site on a species, genus, and family level
- A comparison of the proportion of classified reads to unclassified reads for every sample on both an absolute scale and relative scale
- Sorting for the top 20 contaminants by percentage of reads, and the top contaminants by read count (that meet the threshold of bacteria > 500 or other clades > 50)
    - Plotted as a series of heatmaps, and arranged by date of sampling (x-axis) and total read count (y-axis).

For the site report, we perform the following analyses:
- Sorting for the top 10 most common contaminants across samples by percentage of reads within each sample
    - Calculating their respective read counts, in addition to the total read count, to visualise the distribution of the most abundant contaminants in each sample
- Filtering taxa for those appearing in more than 3 samples with over 5 reads per sample
    - Displaying known taxa on a heatmap by their degree of pathogenicity, according to a reference sheet of contaminants
    - Displaying known taxa on a heatmap by their recognised niche or known environment, according to a reference sheet of contaminants
- Calculating statistical significance of each taxa by comparing all the taxa in the current site to all other sites by performing Welch’s t-test, accounting for possibility of unequal variation
- Displaying a table and heatmap of all taxa which have been positively significant by comparison to the same taxa in at least one other site

All heatmaps display viruses on a genus-level and all other taxa on a family-level.
All spike-ins are removed from analyses.

An example command would be:
```
files=$(cat /neg_control_reports/all_controls_list.txt )
nextflow run main.nf --reports "$files" --metadata "/neg_control_reports/metadata.csv"  --profile docker --site_key site_key.json
```
### Input
```
--reports "$files"
--metadata "/neg_control_reports/metadata.csv"
--profile docker 
--site_key site_key.json
```
--reports provides a list of filepaths to Kraken2 reports and --metadata provides a file with information regarding each sample, including climb_id, site, control_type, and date of collecting/receiving it. --site_key anonymises each site by assigning them aliases. 

### Output
```
pipeline_info	site_reports	summary_reports
```
- summary_reports
    - dataframes
        - ```count_df.csv percentage_df.csv``` - 2 csv files for a percentage and a count-based dataframe for every sample on file
        - ```total_microbe_counts.csv``` - a csv file for total counts found in every site across different clades (e.g. bacteria, virus, fungi)
        - ```G_level_richness.csv S_level_richness.csv F_level_richness.csv``` - 3 csv files for richness data on a species, genus, and family level
        - ```top_perc_heatmap.txt top_count_heatmap.txt``` - 2 text files for the percentage and count-based heatmaps
    - summary_report as an html file
- site_reports
    - dataframes
        - A list of csv files for each site detailing the t-test output for every taxa that is below the p<0.05 threshold
    - site_reports
        - A list of site_reports as an html file
