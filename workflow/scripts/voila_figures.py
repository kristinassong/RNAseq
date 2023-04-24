#!/usr/bin/python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pybedtools import BedTool
import os


# Access snakemake inputs, outputs and params
raw_dir = snakemake.params.indir
out_dir = snakemake.params.outdir
exp_genes = snakemake.input.exp_genes
DE_up = snakemake.input.DE_up
DE_down = snakemake.input.DE_down
outfile = snakemake.output.bar_chart
sno = snakemake.params.sno[0]
sno_interactions_dir = snakemake.params.sno_interactions

SPLICING_EVENTS = ['alt3and5prime','alt3prime','alt5prime','alternate_first_exon','alternate_last_exon','alternative_intron',
                    'cassette','multi_exon_spanning','mutually_exclusive','p_alt3prime','p_alt5prime',
                    'p_alternate_first_exon','p_alternate_last_exon','tandem_cassette']


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
            
                if filename not in ['heatmap.tsv','junctions.tsv','summary.tsv','other.tsv']:
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

def parse_sno(sno_dir,sno_id):

    # Prep sno interactions file for bedtools intersect
    sno_file = os.path.join(sno_dir, sno_id + '.bed')
    sno_df = pd.read_csv(sno_file, sep='\t', 
        names=['seqname', 'target_window_start', 'target_window_end', 'snoid', 'count', 'strand',
            'sno_window_start', 'sno_window_end', 'mean_score', 'min_score', 'max_score', 'target'])

    # Drop unnecessary columns for bedtools intersect
    sno_df.drop(columns=['snoid', 'sno_window_start', 'sno_window_end', 'mean_score', 'min_score', 'max_score'], inplace=True)
    sno_df = sno_df[['seqname', 'target_window_start', 'target_window_end', 'target', 'count', 'strand']]

    # Add 'chr' to seqname if not already present
    sno_df = sno_df.astype({'seqname':'string'})
    if 'chr' not in sno_df.iloc[0]['seqname']:
        sno_df['seqname'] = 'chr' + sno_df['seqname']

    return sno_df

def parse_voila(df):

    # Create a dictionary of lists
    # KEY: event_id  VALUE: list of start and end positions of exons
    pos_dic = {}

    for i in range(len(df)):
        
        event_id = df.iloc[i]['event_id']
        for t in ['reference_exon_coord','spliced_with_coord','junction_coord']:
            if isinstance(df.iloc[i][t], str):
                coord = df.iloc[i][t].split("-") # list
                for el in coord:
                    if el in ['na','']:
                        coord.remove(el)
                    else:
                        el = int(el) # convert strings to ints
                if '1' in coord:
                    coord.remove('1')
                if event_id in pos_dic:
                    pos_dic[event_id] = pos_dic[event_id] + coord
                else:
                    pos_dic[event_id] = coord

    # Reformat data to produce a bed file
    bed_df = pd.DataFrame(columns=['seqid', 'start', 'end', 'gene_name', 'strand'])
    
    info_df = df.drop_duplicates(subset=['event_id']) # df to get info from

    for e in pos_dic.keys():

        # Grep info
        event_info = info_df[info_df['event_id'] == e]
        gene_name = event_info.iloc[0]['gene_name']
        seq_id = event_info.iloc[0]['seqid']
        strand = event_info.iloc[0]['strand']
        new_event = pd.DataFrame({"gene_name":gene_name,'seqid':seq_id,'strand':strand,'start':[min(pos_dic[e])],'end':[max(pos_dic[e])]})
        bed_df = pd.concat([bed_df, new_event])

    bed_df['count']=3 # arbitrary value to meet bed format
    bed_df = bed_df[['seqid', 'start', 'end', 'gene_name', 'count', 'strand']] # reorder columns

    # Add 'chr' to seqid if not already present
    bed_df = bed_df.astype({'seqid':'string'})
    if 'chr' not in bed_df.iloc[0]['seqid']:
        bed_df['seqid'] = 'chr' + bed_df['seqid']

    return bed_df

def splicing_and_sno_binding(sno_interactions_dir, sno, splicing_dir):

    common_count_event_id = {}
    common_count_gene = {}

    sno_df = parse_sno(sno_interactions_dir, sno) # sno interaction file

    # For each splicing event type, find splicing events/genes that are also bound by the snoRNA
    print('Splicing events bound by snoRNA')
    for filename in os.listdir(splicing_dir):
        if '.tsv' in filename and filename.split(".")[0] in SPLICING_EVENTS:
            
            f = os.path.join(splicing_dir, filename)
            df = pd.read_csv(f, sep='\t', usecols=['gene_id','gene_name','seqid','strand','event_id',
                                                    'reference_exon_coord','spliced_with_coord','junction_coord'])

            if len(df) > 0:

                splicing_df = parse_voila(df)

                splicing_bed = BedTool.from_dataframe(splicing_df)
                sno_bed = BedTool.from_dataframe(sno_df)

                df_intersect = splicing_bed.intersect(sno_bed, s=True).to_dataframe()
                if len(df_intersect) > 0:
                    # drop duplicate rows
                    df_intersect.drop_duplicates(inplace=True)
                    pd.set_option('display.max_rows', None) # no limit for displaying df rows
                    print(filename)
                    print(df_intersect)
                    common_count_event_id[filename.split(".")[0]] = len(df_intersect)
                    df_intersect.columns = ['seqid','start','end','gene_name','count','strand']
                    df_intersect_by_gene = df_intersect.drop_duplicates(subset=['gene_name'])
                    common_count_gene[filename.split(".")[0]] = len(df_intersect_by_gene)
                else:
                    common_count_event_id[filename.split(".")[0]] = 0
                    common_count_gene[filename.split(".")[0]] = 0

            else:
                common_count_event_id[filename.split(".")[0]] = 0
                common_count_gene[filename.split(".")[0]] = 0

    df_common_event_id = pd.DataFrame({'Event': common_count_event_id.keys() , 'Count': list(common_count_event_id.values())})
    df_common_event_id = df_common_event_id.sort_values(by=['Count'], ascending=True)
    df_common_gene = pd.DataFrame({'Event': common_count_gene.keys() , 'Count': list(common_count_gene.values())})
    df_common_gene = df_common_gene.sort_values(by=['Count'], ascending=True)

    return df_common_event_id, df_common_gene

# Create various horizontal bar charts
def create_figures(event_df, gene_df, gene_lst, DE_up_file, DE_down_file, outfile, sno, sno_interactions_dir, splicing_dir):

    fig, axs = plt.subplots(3, 2, figsize=(10,9))

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

    # Number of splicing events/genes that are bound by sno
    df_common_event_id, df_common_gene = splicing_and_sno_binding(sno_interactions_dir, sno, splicing_dir)
    axs[2,0].barh(df_common_event_id['Event'],df_common_event_id['Count'], color='dimgray')
    axs[2,0].set_title('Splicing events bound by snoRNA',size=14)
    axs[2,0].set(xlabel="Number of overlaps")
    axs[2,1].barh(df_common_gene['Event'],df_common_gene['Count'], color='dimgray')
    axs[2,1].set_title('Spliced genes bound by snoRNA',size=14)
    axs[2,1].set(xlabel="Number of genes")

    fig.suptitle(sno+' KD',fontsize=16)

    plt.tight_layout()
    plt.savefig(outfile)

    return


# Call functions
event_id_count_df, gene_id_count_df, gene_lists = filter_genes(raw_dir,out_dir,exp_genes)
create_figures(event_id_count_df, gene_id_count_df, gene_lists, DE_up, DE_down, outfile, sno, sno_interactions_dir, out_dir)