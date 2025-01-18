#!/usr/bin/python3

import pandas as pd
from gtfparse import read_gtf
import os

tx2gene = snakemake.input.tx2gene
gtf = snakemake.input.gtf
out_tpm = snakemake.output.tpm
out_counts = snakemake.output.counts


final_df_tpm = pd.DataFrame()
final_df_counts = pd.DataFrame()
cycle = 1
sample_list = []

for q in snakemake.input.quant:

    data_tpm = pd.read_csv(q, sep='\t', usecols=['target_id','tpm'])
    data_counts = pd.read_csv(q, sep='\t', usecols=['target_id','est_counts'])

    # get sample name
    sample = os.path.basename(os.path.dirname(q))

    # reformat dataframe
    data_tpm.set_index('target_id', inplace=True)
    data_tpm.rename(columns={"tpm":sample}, inplace=True)
    data_counts.set_index('target_id', inplace=True)
    data_counts.rename(columns={"est_counts":sample}, inplace=True)

    # Merge dataframes
    if cycle == 1:
        final_df_tpm = data_tpm
        final_df_counts = data_counts
    else:
        final_df_tpm = pd.merge(final_df_tpm, data_tpm, left_index=True, right_index=True)
        final_df_counts = pd.merge(final_df_counts, data_counts, left_index=True, right_index=True)
    
    sample_list += [sample]
    cycle+=1

# transcript ID --> gene ID
ids = pd.read_csv(tx2gene, sep='\t', names=['transcript','gene'])
ids.set_index('transcript', inplace=True, drop=False)

final_df_tpm = pd.merge(final_df_tpm, ids, left_index=True, right_index=True)
final_df_tpm.set_index('gene', inplace=True)
final_df_counts = pd.merge(final_df_counts, ids, left_index=True, right_index=True)
final_df_counts.set_index('gene', inplace=True)

# Add gene name
df_gtf = read_gtf(gtf)
id_name = df_gtf[['gene_id','gene_name']].to_pandas().drop_duplicates(ignore_index=True)
#df_gtf = pd.read_csv(gtf, sep='\t') # for annotations in tsv format
#id_name = df_gtf[['gene_id','gene_name']].drop_duplicates(ignore_index=True)

index = final_df_tpm.index.tolist()
names =[]

for i in range(len(index)):
    id = index[i]
    name = id_name[id_name['gene_id']==id].iloc[0]['gene_name']
    names.append(name)

final_df_tpm['gene_name'] = names
final_df_counts['gene_name'] = names

cols_tpm = final_df_tpm.columns.tolist()
cols_tpm = cols_tpm[-1:] + cols_tpm[:-1]
cols_tpm = cols_tpm[-1:] + cols_tpm[:-1]
final_df_tpm = final_df_tpm[cols_tpm]

cols_counts = final_df_counts.columns.tolist()
cols_counts = cols_counts[-1:] + cols_counts[:-1]
cols_counts = cols_counts[-1:] + cols_counts[:-1]
final_df_counts = final_df_counts[cols_counts]

# Write to file
final_df_tpm.to_csv(out_tpm, sep='\t')
final_df_counts.to_csv(out_counts, sep='\t')
