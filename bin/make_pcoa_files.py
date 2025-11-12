#!/usr/bin/env python

import pandas as pd
import os
import argparse
import json
from natsort import natsorted

from utils import define_heatmap_datasets, make_heatmap_df, make_date_list, filter_multiples, make_empty_df

def sort_by_date(count_df, date_list):
    #add and sort samples by date
    count_df.loc[len(count_df)] = date_list
    row_number = count_df.index.get_loc(count_df[count_df["Scientific_Name"] == "Date"].index[0])
    count_df.iloc[row_number, 2:-1] = pd.to_datetime(count_df.iloc[row_number, 2:-1], format='%Y-%m-%d')
    #Sort all samples by date
    count_df.iloc[row_number, 2:-1] = natsorted(count_df.iloc[row_number, 2:-1])

    #get timestamps line and turn to list
    timestamp_df = count_df.drop(columns="Counts_Overall")
    row_number = timestamp_df.index.get_loc(timestamp_df[timestamp_df["Scientific_Name"] == "Date"].index[0])
    date_df = timestamp_df.iloc[row_number]
    dates = list(date_df)
    #remove timestamps
    count_df = count_df[~count_df.Scientific_Name.str.contains("Date")]
    return count_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse text files and a CSV file.")
    parser.add_argument('--reports', nargs='+', help="List of text files", required=True)
    parser.add_argument('--metadata', help="CSV file path", required=True)
    parser.add_argument('--site_key', help="JSON file specifying site name to number", required=True)
    parser.add_argument('--reference', help="excel file for dictionary of contaminants", required=True)
    parser.add_argument('--r_dir', help="directory for plots to be made in R code", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    r_path = args.r_dir
    os.makedirs(r_path, exist_ok=True) #make directory path into real directory


    #site name to number
    site_key = {}
    with open(args.site_key, 'r') as f:
        site_key = json.load(f)

    grouped_metadata = metadata.groupby(['site', 'control_type_details'])

    mscape_datasets, mscape_samples, mscape_dates = define_heatmap_datasets(grouped_metadata, site_key)

    exclude_cols = ["Scientific_Name"]

    present_dfs = []
    diff_loop = 0
    for needed_samples in mscape_samples:
        #for dataset in mscape_datasets:
        count_df, og_df = make_heatmap_df(needed_samples, reports, "count")

        sample_date_df = mscape_dates[diff_loop]
        date_list = make_date_list(sample_date_df)

        #name_line = ["domain", "scientific_name"] + sample_names[loop] + ["counts_overall"]
        #ßcount_df.loc[len(count_df)] = name_line
        #add and sort samples by date
        count_df.loc[len(count_df)] = date_list
        row_number = count_df.index.get_loc(count_df[count_df["Scientific_Name"] == "Date"].index[0])
        count_df.iloc[row_number, 2:-1] = pd.to_datetime(count_df.iloc[row_number, 2:-1], format='%Y-%m-%d')
        #Sort all samples by date
        count_df.iloc[row_number, 2:-1] = natsorted(count_df.iloc[row_number, 2:-1])

        #get timestamps line and turn to list
        timestamp_df = count_df.drop(columns="Counts_Overall")
        row_number = timestamp_df.index.get_loc(timestamp_df[timestamp_df["Scientific_Name"] == "Date"].index[0])
        date_df = timestamp_df.iloc[row_number]
        dates = list(date_df)
        #remove timestamps
        count_df = count_df[~count_df.Scientific_Name.str.contains("Date")]

        count_df = count_df.drop(columns=["Domain", "Counts_Overall"])

        include_cols = [col for col in count_df.columns if col not in exclude_cols]

        # Boolean indexing to filter rows showing only genus/family in the rank column
        filtered_df = count_df[count_df[include_cols].apply(lambda x: filter_multiples(x), axis=1)]

        present_dfs.append(filtered_df)
        diff_loop += 1
    
    # Merge the DataFrames on a specific column
    all_df = pd.concat(present_dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
    all_df = all_df.groupby('Scientific_Name', as_index=False).first()
    all_df.fillna(0, inplace=True)

    # Sort by count
    numeric = all_df.drop(columns=['Scientific_Name'])
    numeric_sums = numeric.sum(axis=1)
    all_df["counts_overall"] = numeric_sums
    all_df = all_df.sort_values(by="counts_overall", ascending=False)

    #find unlabeled contaminants and save as csv file
    niche_table = pd.read_excel(args.reference)
    niche_table.fillna("NaN", inplace=True)

    edit_df = pd.DataFrame(columns=niche_table.columns)

    niche_table.set_index("Taxa", inplace=True)
    nan_rows = niche_table.isna().any(axis=1)

    counts = []
    for index, row in all_df.iterrows():
        taxa_name = row["Scientific_Name"]
        
        if taxa_name in niche_table.index:
            niche_list = list(niche_table.loc[taxa_name])
            pathogenicity = niche_list[0]
            niches = [x for x in niche_list if x is not pathogenicity]
            all_niches = [taxa_name] + niche_list
            if pathogenicity == "NaN" or niches == ["NaN", "NaN", "NaN", "NaN"]:
                edit_df = pd.concat([pd.DataFrame([all_niches], columns=edit_df.columns), edit_df], ignore_index=True)
                counts.append(row["counts_overall"])
        else:
            niche_list = [taxa_name, "NaN", "NaN", "NaN", "NaN", "NaN"]
            edit_df = pd.concat([pd.DataFrame([niche_list], columns=edit_df.columns), edit_df], ignore_index=True)
            counts.append(row["counts_overall"])

    edit_df = edit_df.iloc[::-1]
    edit_df["Counts_Overall"] = counts

    col_order = [col for col in edit_df.columns if col != "Counts_Overall"] + ["Counts_Overall"]
    edit_df = edit_df[col_order]

    edit_df = edit_df.replace('NaN', '', regex=True)
    edit_df["All taxa here are present in the samples, and either missing a pathogenicity label, or missing entries in all 3 niche labels (Lab, Human, and Industry). Please update them in bin/contaminant_literature.xlsx accordingly. If there are multiple entries in one category, please make sure to put no spaces between the commas when listing them."] = ''
    edit_df.to_csv(os.path.join(r_path, "unlabeled_contaminants.csv"), index=False, na_rep='')

    all_df = all_df.drop(columns="counts_overall")
    # Depict site-specific plots
    site_loop = 0
    for dataset_name in mscape_datasets:

        #make new dfs using information from our previous merged all_df to ensure empty taxa are also accounted for
        new_dfs = []
        for df in present_dfs:
            new_df = all_df[df.columns]
            new_dfs.append(new_df)

        current_df = new_dfs[site_loop]
        if len(current_df.index) > 0:
            current_df.set_index("Scientific_Name", inplace=True)
        else:
            current_df = make_empty_df(current_df)
            current_df.set_index("Scientific_Name", inplace=True)

        other_dfs = [df for df in new_dfs if df is not current_df]
        # Merge the DataFrames on a specific column
        other_df = pd.concat(other_dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
        if len(other_df.index) > 0:
            #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
            other_df = other_df.groupby('Scientific_Name', as_index=True).first()
        else:
            other_df = make_empty_df(other_df)
            current_df.set_index("Scientific_Name", inplace=True)
        other_df.fillna(0, inplace=True)


        transposed_current = current_df.transpose()
        transposed_current["site"] = dataset_name
        transposed_other = other_df.transpose()
        transposed_other["site"] = "other_sites"

        pcoa_df = pd.concat([transposed_current, transposed_other], axis=0)
        pcoa_df.fillna(value=0, inplace=True)

        # Sum columns
        columns_to_sum = pcoa_df.drop(columns=['Scientific_Name', 'site'])
        column_sums = columns_to_sum.sum(axis=1)
        pcoa_df["Counts_Overall"] = column_sums
        pcoa_df = pcoa_df[~pcoa_df["Counts_Overall"].astype(float).isin([0])]
        pcoa_df = pcoa_df.drop(columns="Counts_Overall")
        
        pcoa_df.to_csv(os.path.join(r_path, f"{dataset_name}.pcoa.txt"), sep=',', index=False)
        site_loop += 1



