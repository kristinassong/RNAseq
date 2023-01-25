#!/usr/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import os

raw_dir = snakemake.params.indir
out_dir = snakemake.params.outdir

count_dic = {}

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
        
            if filename not in ['heatmap.tsv','junctions.tsv','summary.tsv']:
                event = filename.split('.')[0]
                count_dic[event] = len(pd.unique(voila_filtered['event_id']))

count_df = pd.DataFrame({'Event': count_dic.keys() , 'Count': list(count_dic.values())})
count_df = count_df.sort_values(by=['Count'], ascending=True)

# Create a bar chart showing the number of splicing events per event type
plt.figure(figsize=(6,4))

plt.barh(count_df['Event'],count_df['Count'], color='dimgray')
#bars = pd.DataFrame({'num_of_events':count_df['Count']}, index=count_df['Event'])

#bars['num_of_events'].plot(kind="barh", width=0.7)
plt.title('Number of splicing events by event type')
plt.xlabel('Number of events')

plt.tight_layout()
plt.savefig(snakemake.output.bar_chart)