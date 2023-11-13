#!/usr/bin/python3

import pandas as pd
from pybedtools import BedTool
import argparse
import os

SPLICING_EVENTS = ['SE.MATS.JC.txt','A5SS.MATS.JC.txt','A3SS.MATS.JC.txt','MXE.MATS.JC.txt','RI.MATS.JC.txt']


def filter_by_threshold(df):
    # Filter rMATS output by PValue, FDR and IncLevelDifference
    filtered_df = df[df['PValue']<0.05]
    filtered_df = filtered_df[filtered_df['FDR']<0.05]
    filtered_df = filtered_df[filtered_df['IncLevelDifference'].abs()>0.20]
    return filtered_df

def filter_by_tpm(splicing_df,tpm_df):
    # Keep genes that are expressed in at least one condition (NC5 & SNORD22 KD in kallisto TPM matrix)
    expressed_df = pd.DataFrame(columns=splicing_df.columns)
    # List of unique spliced genes
    uniq = set(splicing_df['GeneID'].values.tolist())
    for gene in uniq:
        gene_tpm_df = tpm_df[tpm_df['gene']==gene]
        if len(gene_tpm_df)>0: # TPM >=1 in at least one sample
            gene_splicing_df = splicing_df[splicing_df['GeneID']==gene]
            expressed_df = pd.concat([expressed_df,gene_splicing_df],ignore_index=True)
    return expressed_df

def parse_splicing(df,type):
    # Transform splicing df to bed
    df['score']=3 # arbitary score for bed format
    if type == 'SE' or type == 'MXE':
        df = df[['ID','chr','upstreamES','downstreamEE','geneSymbol','score','strand','GeneID','PValue','FDR','IncLevelDifference']]
        df.rename(columns={"upstreamES": "start", "downstreamEE": "end"},inplace=True)
    elif type == 'RI':
        df = df[['ID','chr','riExonStart_0base','riExonEnd','geneSymbol','score','strand','GeneID','PValue','FDR','IncLevelDifference']]
        df.rename(columns={"riExonStart_0base": "start", "riExonEnd": "end"},inplace=True)
    else: # A3SS or A5SS
        pos_df = df[df['strand']=='+']
        neg_df = df[df['strand']=='-']
        if type == 'A3SS':
            pos_df = pos_df[['ID','chr','flankingES','shortEE','geneSymbol','score','strand','GeneID','PValue','FDR','IncLevelDifference']]
            pos_df.rename(columns={"flankingES": "start", "shortEE": "end"},inplace=True)
            neg_df = neg_df[['ID','chr','longExonStart_0base','flankingEE','geneSymbol','score','strand','GeneID','PValue','FDR','IncLevelDifference']]
            neg_df.rename(columns={"longExonStart_0base": "start", "flankingEE": "end"},inplace=True)
        else: # A5SS
            pos_df = pos_df[['ID','chr','longExonStart_0base','flankingEE','geneSymbol','score','strand','GeneID','PValue','FDR','IncLevelDifference']]
            pos_df.rename(columns={"longExonStart_0base": "start", "flankingEE": "end"},inplace=True)
            neg_df = neg_df[['ID','chr','flankingES','shortEE','geneSymbol','score','strand','GeneID','PValue','FDR','IncLevelDifference']]
            neg_df.rename(columns={"flankingES": "start", "shortEE": "end"},inplace=True)
        df = pd.concat([pos_df,neg_df],ignore_index=True)
    return df

def binding_ovlp(df1,df2):
    # bedtools intersect
    bed1 = BedTool.from_dataframe(df1)
    bed2 = BedTool.from_dataframe(df2)
    intersection_df = bed1.intersect(bed2,s=True,u=True).to_dataframe(disable_auto_names=True, names=['chr','start','end','ID','score','strand'])
    return intersection_df

def get_full_entry(bed_df,org_df):
    # Get full information on filtered splicing events
    final_df = pd.DataFrame(columns=org_df.columns)
    for i in range(len(bed_df)):
        id = bed_df.iloc[i]['ID']
        row = org_df[org_df['ID']==id]
        final_df = pd.concat([final_df,row],ignore_index=True)
    return final_df

def spliced_and_bound(spliced_df, event_type, binding_factors_dict, working_dir):
    # bedtools intersect AS splicing events with binding interactions
    spliced_df_bed6 = spliced_df[['chr','start','end','ID','score','strand']] # reformat rMATS output to bed6 format
    
    for key in binding_factors_dict:
        intersection_df = binding_ovlp(spliced_df_bed6, binding_factors_dict[key])

        # get original entry from rMATS output
        event_file = os.path.join(working_dir, event_type + '.tsv')
        outpath = os.path.join(working_dir, event_type+'_'+key+'.tsv')
        get_full_entry(intersection_df, pd.read_csv(event_file, sep='\t')).to_csv(outpath,sep='\t',index=None)
    return

binding_bool = snakemake.params.binding

if binding_bool:
    # Access snakemake inputs, outputs and params
    working_dir = snakemake.params.working_dir
    outfile = snakemake.output.result
    infile = snakemake.input.binding_factors

    # Create an upset plot showing AS genes bound by snoRNAs and RBPs

    # Get dfs of all binding factors
    binding_factors_paths = open(infile, 'r').readline().split(",")
    binding_factors_dict = {}
    for b in binding_factors_paths:
        binding_factor = os.path.basename(b).split(".")[0]
        binding_factors_dict[binding_factor] = pd.read_csv(b, sep='\t')

    # Find AS events bound by binding factors
    for event in SPLICING_EVENTS:
        event_type = event.split(".")[0]
        event_file = os.path.join(working_dir, event_type + '.tsv')
        event_df = parse_splicing(pd.read_csv(event_file, sep='\t'), event_type) # AS events df reformatted
        spliced_and_bound(event_df, event_type, binding_factors_dict, working_dir)
        
else:
    # Access snakemake inputs, outputs and params
    raw_dir = snakemake.params.indir # splicing event dir
    out_dir = snakemake.params.outdir
    tpm = snakemake.params.tpm

    # Kallisto TPM matrix
    tpm_df = pd.read_csv(tpm,sep='\t')

    # For each splicing event type, FILTER by multiple thresholds
    for event in SPLICING_EVENTS:
        print(event)
        file = os.path.join(raw_dir, event)
        event_type = event.split(".")[0]
        df = pd.read_csv(file, sep='\t')
        filtered = filter_by_threshold(df)
        filtered = filter_by_tpm(filtered,tpm_df)

        # Splicing events filtered by threshold and TPM
        print('Original filtered splicing events')
        filtered.to_csv(os.path.join(out_dir,event_type+'.tsv'),sep='\t',index=None)