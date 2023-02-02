#!/usr/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import os


# Access snakemake inputs, outputs and params
raw_dir = snakemake.params.indir
out_dir = snakemake.params.outdir
exp_genes = snakemake.input.exp_genes
DE_up = snakemake.input.DE_up
DE_down = snakemake.input.DE_down
outfile = snakemake.output.bar_chart
sno = snakemake.params.sno[0]


def filter_genes(raw_dir,out_dir,exp_genes):
    # Group by event_id and gene_id
    event_id_count_dic = {}
    gene_id_count_dic = {}
    gene_lists = {}

    for filename in os.listdir(raw_dir):
        f = os.path.join(raw_dir, filename)
        # checking if it is a file
        if os.path.isfile(f) and ".tsv" in filename:

            voila_df = pd.read_csv(f, sep='\t', comment='#')

            if len(voila_df) > 0:
                genes = pd.read_csv(exp_genes, sep='\t', usecols=['gene'])
                outfile = os.path.join(out_dir, filename)

                # Drop genes that are not in the TPM filtered list
                voila_filtered = voila_df[voila_df.gene_id.isin(genes.gene)] 

                # Write to file
                voila_filtered.to_csv(outfile, sep='\t', index=False)
            
                if filename not in ['heatmap.tsv','junctions.tsv','summary.tsv']:
                    event = filename.split('.')[0]
                    event_id_count_dic[event] = len(pd.unique(voila_filtered['event_id']))
                    unique_genes = pd.unique(voila_filtered['gene_id'])
                    gene_id_count_dic[event] = len(unique_genes)
                    gene_lists[event] = list(unique_genes)

    event_id_count_df = pd.DataFrame({'Event': event_id_count_dic.keys() , 'Count': list(event_id_count_dic.values())})
    event_id_count_df = event_id_count_df.sort_values(by=['Count'], ascending=True)
    gene_id_count_df = pd.DataFrame({'Event': gene_id_count_dic.keys() , 'Count': list(gene_id_count_dic.values())})
    gene_id_count_df = gene_id_count_df.sort_values(by=['Count'], ascending=True)

    return event_id_count_df, gene_id_count_df, gene_lists

def splicing_and_DE(spliced_genes, DE_genes):

    common_count = {}

    # List of genes that are DE
    DE_df = pd.read_csv(DE_genes, sep='\t', usecols=['gene'])
    
    # For each splicing event type, find genes that are also DE
    for e in list(spliced_genes.keys()):
        common_count[e] = len(set(spliced_genes[e]) & set(DE_df['gene'].values.tolist()))

    df_common = pd.DataFrame({'Event': common_count.keys() , 'Count': list(common_count.values())})
    df_common = df_common.sort_values(by=['Count'], ascending=True)

    return df_common

# Create various horizontal bar charts
def create_figures(event_df, gene_df, gene_lst, DE_up_file, DE_down_file, outfile, sno):

    fig, axs = plt.subplots(2, 2, figsize=(10,9))

    # Number of splicing events by event_id
    axs[0,0].barh(event_df['Event'],event_df['Count'], color='dimgray')
    axs[0,0].set_title('Splicing events',size=14)
    axs[0,0].set(xlabel="Number of events")

    # Number of genes per splicing event
    axs[0,1].barh(gene_df['Event'],gene_df['Count'], color='dimgray')
    axs[0,1].set_title('Spliced genes',size=14)
    axs[0,1].set(xlabel="Number of genes")

    # Number of genes that are spliced and DE (downregulated) 
    df_spliced_DE_down = splicing_and_DE(gene_lst, DE_down_file)
    axs[1,0].barh(df_spliced_DE_down['Event'],df_spliced_DE_down['Count'], color='dimgray')
    axs[1,0].set_title('Spliced and downregulated genes',size=14)
    axs[1,0].set(xlabel="Number of genes")

    # Number of genes that are spliced and DE (upregulated)
    df_spliced_DE_up = splicing_and_DE(gene_lst, DE_up_file)
    axs[1,1].barh(df_spliced_DE_up['Event'],df_spliced_DE_up['Count'], color='dimgray')
    axs[1,1].set_title('Spliced and upregulated genes',size=14)
    axs[1,1].set(xlabel="Number of genes")

    fig.suptitle(sno+' KD',fontsize=16)

    plt.tight_layout()
    plt.savefig(outfile)

    return


# Call functions
event_id_count_df, gene_id_count_df, gene_lists = filter_genes(raw_dir,out_dir,exp_genes)
create_figures(event_id_count_df, gene_id_count_df, gene_lists, DE_up, DE_down, outfile, sno)