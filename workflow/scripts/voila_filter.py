#!/usr/bin/python3

import pandas as pd
import os

voila_file = snakemake.input.tsv
voila_df = pd.DataFrame()
if 'deltapsi' in voila_file:
    voila_df = pd.read_csv(snakemake.input.tsv, header=10, sep='\t')
else:
    voila_df = pd.read_csv(snakemake.input.tsv, header=7, sep='\t')

genes = pd.read_csv(snakemake.input.exp_genes, sep='\t', usecols=['gene'])
outfile = snakemake.output.filtered_tsv

# Drop genes that are not in the TPM filtered list
voila_filtered = voila_df[voila_df.gene_id.isin(genes.gene)] 

# Write to file
voila_filtered.to_csv(outfile, sep='\t', index=False)
