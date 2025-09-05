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
    parser.add_argument('--r_dir', help="directory for plots to be made in R code", required=True)
    args = parser.parse_args()

    reports = args.reports
    metadata = pd.read_csv(args.metadata)
    output_path = args.r_dir
    os.makedirs(output_path, exist_ok=True) #make directory path into real directory

    #site name to number
    site_key = {}
    with open(args.site_key, 'r') as f:
        site_key = json.load(f)

    grouped_metadata = metadata.groupby(['site', 'control_type_details'])

    datasets, mscape_datasets, samples, mscape_samples, sample_dates, mscape_dates = define_heatmap_datasets(grouped_metadata, site_key)

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
    
    # Depict site-specific plots
    site_loop = 0
    for dataset_name in mscape_datasets:
        # Merge the DataFrames on a specific column
        all_df = pd.concat(present_dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
        #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
        all_df = all_df.groupby('Scientific_Name', as_index=False).first()
        all_df.fillna(0, inplace=True)

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

        pcoa_df.to_csv(os.path.join(output_path, f"{dataset_name}.pcoa.txt"), sep='\t', index=False)
        site_loop += 1



