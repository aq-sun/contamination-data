#!/usr/bin/env python
import pandas as pd

spikeins = {
        12242:["Virgaviridae","Tobamovirus","Tobacco_mosaic_virus"]
    }

def get_label(ids, site_key, run_id=None):
    ids_list = [id.lower() for id in ids]

    ids_list_anon = []
    for id in ids_list:
        if id in site_key:
            ids_list_anon.append(site_key[id])
        else:
            ids_list_anon.append(id)
    ids_list = ids_list_anon

    ids_list = ''.join(ids_list)
    if "public" in ids and run_id:
        ids_list = ids_list.replace('public', 'public_'+run_id) #add run_id to name
    ids_list = ids_list.replace('water_extraction_control', '(water)')
    ids_list = ids_list.replace('resp_matrix_mc110', '(matrix)')
    return ids_list

def define_datasets(grouped_metadata, site_key):
    datasets = [] #total datasets
    samples = [] #samples grouped by datasets

    for sets in grouped_metadata:
            ids = list(sets[0])
            if "public" in ids[0].lower(): #if dataset is public
                subgroup = sets[1].groupby(['run_id'])
                for subset in subgroup:
                    run_id = ''.join(subset[0])
                    #turn scientific_name from a list to a string
                    ids_list = get_label(ids, site_key, run_id=run_id)

                    datasets.append(ids_list)
                    table = subset[1] #list of all ids in sub_dataset
                    samples.append(list(table['climb_id'])) #climb id
            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id
    return datasets, samples

def define_heatmap_datasets(grouped_metadata, site_key):
    #group by site
    datasets = []
    mscape_datasets = []
    samples = []
    mscape_samples = []
    sample_dates = []
    mscape_dates = []

    for sets in grouped_metadata:
            ids = list(sets[0])
            if "public" in ids[0].lower(): #if dataset is public
                subgroup = sets[1].groupby(['run_id'])
                for subset in subgroup:
                    run_id = ''.join(subset[0])
                    #turn scientific_name from a list to a string
                    ids_list = get_label(ids, site_key, run_id=run_id)
                    datasets.append(ids_list)
                    table = subset[1] #list of all ids in sub_dataset
                    samples.append(list(table['climb_id'])) #climb id
                    dates_columns = ['collection_date', 'received_date']
                    sample_dates.append(table[dates_columns]) #sample dates

            else:
                #turn scientific_name from a list to a string
                ids_list = get_label(ids, site_key)
                datasets.append(ids_list)
                mscape_datasets.append(ids_list)
                table = sets[1] #list of all ids in dataset
                samples.append(list(table['climb_id'])) #climb id
                mscape_samples.append(list(table['climb_id']))
                dates_columns = ['collection_date', 'received_date']
                sample_dates.append(table[dates_columns]) #sample dates
                mscape_dates.append(table[dates_columns])
    
    return datasets, mscape_datasets, samples, mscape_samples, sample_dates, mscape_dates

def convert_to_numeric(column):
        if column.name not in ['Rank', 'Scientific_Name', 'Domain']:
            return pd.to_numeric(column, errors='coerce')
        else:
            return column

def make_date_list(sample_date_df):

    date_list = ['', 'Date']
    for index, row in sample_date_df.iterrows():
        if row['collection_date'] == row['collection_date']:
            date_list.append(row['collection_date'])
        elif row['received_date'] == row['received_date']:
            date_list.append(row['received_date'])
        else:
            date_list.append('NaN')
    date_list.append(100)

    return date_list

def split_dfs(mscape_datasets, total_df):
    mscape_dfs = []
    #make specific dfs for all dataframes in mscape_datasets
    for dataset_name in mscape_datasets:
        short_name = dataset_name.split("(")[0]
        short_name_p = short_name + "\("
        category = dataset_name.split("(")[1].replace(")","")
        name_col = ["Domain", "Scientific_Name"]
        df = total_df.loc[:, total_df.columns.str.contains(short_name_p) & total_df.columns.str.contains(category)]
        df[name_col] = total_df[name_col].copy()
        mscape_dfs.append(df)
    
    loop_count = 0
    sorted_mscapes = []
    sorted_m_counts = []

    for mscape_df in mscape_dfs:
        #Get rid of dataset as prefixes
        mscape_df = mscape_df.rename(columns={c: c.replace(f"{mscape_datasets[loop_count]}_", "") for c in mscape_df.columns if c not in ['Scientific_Name', 'Domain']})

        final_df = mscape_df[~mscape_df.Scientific_Name.str.contains("Date")]
        sorted_mscapes.append(final_df)
        loop_count += 1
    
    return sorted_mscapes

def save_labelled_df(total_df, all_dates, output_path, df_type):
        labelled_df = total_df.copy()
        labelled_df.loc[len(labelled_df)] = all_dates

        set_row = ["", "Dataset"]
        count_row = ["", "Total Read Count"]
        name_row = ["Domain", "Scientific_Name"]
        set_name = list(labelled_df.columns)
        set_row.extend([i.split('_', 1)[0] for i in set_name[2:]])
        name_row.extend([i.split('_', 1)[1] for i in set_name[2:]])
    
        labelled_df.columns = name_row #label each sample column by sample name
        labelled_df.loc[len(labelled_df)] = set_row #label the dataset the sample is from

        #Change datetime timestamp back to string
        row_number = labelled_df.index.get_loc(labelled_df[labelled_df["Scientific_Name"] == "Date"].index[0])
        labelled_df.iloc[row_number, 2:] = [i.strftime('%Y-%m-%d') for i in labelled_df.iloc[row_number, 2:]]
        labelled_df.to_csv(f"{output_path}/dataframes/top_{df_type}_heatmap.txt", sep='\t', index=False)

def make_count_and_perc_dfs(needed_samples, reports, df_type):
    # Initialize an empty list to store dataframes
    dfs = []

    #   Loop over each kraken file in reports directory
    for sample in needed_samples:
        found = False
        for filename in reports:
            if sample in filename:
                found = True

                #open new file and read it line by line
                file = open(filename)

                #create 3 lists for count of sequences then rank and scientific name
                perc_seqs = []
                read_counts = []
                taxid = []
                rank = []
                sci_name = []
                domain = []

                #only start reading each file when it starts listing our domain
                current_domain = None

                # Iterate over each line in the file
                for line in file:
                    if line.startswith("%"):
                        continue

                    # Split the line into columns
                    columns = line.split()

                    #record the current scientific name of this line
                    current_name = columns[5]
                    current_rank = columns[3]

                    # Check for the microbial category in scientific name column and set the flag to start readin
                    if current_name in ["Sar", "Discoba", "Metamonada"]:
                        current_domain = "Protists"
                    elif current_name == "Fungi":
                        current_domain = "Fungi"
                    elif current_name == "Bacteria":
                        current_domain = "Bacteria"
                    elif current_name == "Viruses":
                        current_domain = "Viruses"
                    elif current_name == "Archaea":
                        current_domain = "Archaea"
                    elif current_name == "Eukaryota":
                        current_domain = "Eukaryota"
                    elif (current_domain == "Protists" and current_rank == "D1") \
                        or (current_domain == "Fungi" and current_rank == "K") \
                        or (current_domain in ["Bacteria","Viruses","Archaea","Eukaryota"] and current_rank == "D"):
                        current_domain = None

                    perc_seqs.append(columns[0])
                    read_counts.append(columns[1])
                    taxid.append(columns[4])
                    rank.append(current_rank)
                    sci_name.append('_'.join(columns[5:])) #this turns scientific name (which sometimes have multiple words) into a list within a list
                    domain.append(current_domain)

                # Turn the four lists into a dataframe, using sample ID in place of "% of seqs" or "read counts", depending on whether you want counts or percentages
                if df_type == "perc":
                    df = pd.DataFrame({sample: perc_seqs, "Rank": rank, "Scientific_Name": sci_name, "Domain":domain, "Taxon_ID": taxid})
                else:
                    df = pd.DataFrame({sample: read_counts, "Rank": rank, "Scientific_Name": sci_name, "Domain":domain, "Taxon_ID": taxid})


                #add the new dataframe to the list of dataframes
                dfs.append(df)

        if not found:
            print(f"No kraken report for sample {sample} has been provided!")

    # Merge the DataFrames on a specific column
    merged_df = pd.concat(dfs, axis = 0, join= "outer")  # Change join to 'outer' for outer join
    #merge on scientific name - this means that no scientific name is repeated if it's present in different dataframes
    merged_df = merged_df.groupby('Scientific_Name', as_index=False).first()
    # Merge the DataFrames on a specific column
    merged_df = merged_df.apply(convert_to_numeric).fillna(0)
    return merged_df

def make_heatmap_df(needed_samples, reports, df_type):
    df = make_count_and_perc_dfs(needed_samples, reports, df_type)
    no_taxid_df = df.drop(columns=["Taxon_ID"])

    #remove spikeins
    for spike in spikeins:
        no_spike_df = no_taxid_df.loc[~no_taxid_df["Scientific_Name"].astype(str).isin(spikeins[spike])]

    #Filtering the dataframe    
    # Define keywords and columns to search
    if df_type in ("perc", "count"):
        keywords = ['G'] #the taxonomy level I want to filter for
    elif df_type == "thresh":
        keywords = ['G', 'R', 'U']
    columns_to_search = ['Rank']
    
    #filter for virus families
    virus_df = no_spike_df[no_spike_df['Domain'] == "Viruses"]
    virusfam_df = virus_df[virus_df['Rank'] == "F"]

    # Boolean indexing to filter rows showing only genus/family in the rank column
    genus_df = no_spike_df[no_spike_df[columns_to_search].apply(lambda x: x.isin(keywords).any(), axis=1)]
    genus_df = genus_df[~genus_df.Domain.str.contains("Viruses", na=False)] #remove viruses from dataframe

    #Add viruses back in as family
    genusfam = pd.concat([genus_df, virusfam_df], axis = 0)

    # Select columns to sum (excluding 'Rank' and 'Scientific Name')
    columns_to_sum = genusfam.drop(columns=['Rank', 'Scientific_Name', 'Domain'])
    # Calculate the sum of values in each row
    column_sums = columns_to_sum.sum(axis=1)
    # Count the number of columns
    num_columns = columns_to_sum.shape[1]
    #either get a sum of all seqs, or normalise it by dividing by num_columns
    if df_type == "perc":
        genusfam["Perc_Seqs_Overall"] = column_sums/num_columns
        filtered_df = genusfam[~genusfam["Perc_Seqs_Overall"].astype(float).isin([0])] # Boolean indexing to filter rows in % of seqs overall that contain 0
    elif df_type == "thresh" or df_type == "count":
        genusfam["Counts_Overall"] = column_sums
        filtered_df = genusfam[~genusfam["Counts_Overall"].astype(float).isin([0])]

    
    # Rearrange column 'Scientific Name' to the first position
    wordy_columns = ['Scientific_Name', 'Domain', "Counts_Overall", "Perc_Seqs_Overall"]
    # check if columns are not in the wordy_columns list
    if df_type == "perc":
        column_order = ['Domain'] + ['Scientific_Name'] + [col for col in filtered_df.columns if col not in wordy_columns] + ["Perc_Seqs_Overall"]
    elif df_type == "thresh" or df_type == "count":
        column_order = ['Domain'] + ['Scientific_Name'] + [col for col in filtered_df.columns if col not in wordy_columns] + ["Counts_Overall"]
    filtered_df = filtered_df[column_order]
    
    #drop column “rank” since they are all Gs 
    no_rank_df = filtered_df.drop(columns=['Rank'])
    return no_rank_df, df

def filter_multiples(x):
    filtered_x = x[x>=5] #find columns with over or equal to 5 counts
    return len(filtered_x) >= 3 #filter for rows with >= 5 counts in >=3 columns

def make_empty_df(df):
    empty_row = [0]
    dummy_loop = 1
    while dummy_loop < len(df.columns):
        empty_row.append(0)
        dummy_loop += 1

    dummy_df = pd.DataFrame(columns=df.columns)
    dummy_df = pd.concat([pd.DataFrame([empty_row], columns=dummy_df.columns), dummy_df], ignore_index=True)
    
    for col_type in dummy_df.columns: #label string columns with "none"
        if col_type in ["Domain", "Scientific_Name", "Niche"]:
            dummy_df[col_type] = "None"

    return dummy_df
