#!/usr/bin/python3

import pandas as pd
import os

raw_dir = snakemake.params.indir
out_dir = snakemake.params.outdir

for filename in os.listdir(raw_dir):
    f = os.path.join(raw_dir, filename)
    # checking if it is a file
    if os.path.isfile(f):

        voila_df = pd.read_csv(f, sep='\t', comment='#')

        if len(voila_df) > 0:
            genes = pd.read_csv(snakemake.input.exp_genes, sep='\t', usecols=['gene'])
            outfile = os.path.join(out_dir, filename)

            # Drop genes that are not in the TPM filtered list
            voila_filtered = voila_df[voila_df.gene_id.isin(genes.gene)] 

            # Write to file
            voila_filtered.to_csv(outfile, sep='\t', index=False)
