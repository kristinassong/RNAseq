#!/usr/bin/python3

import pandas as pd
from pybedtools import BedTool
import os


# Access snakemake inputs, outputs and params
raw_dir = snakemake.params.indir # splicing event dir
out_dir = snakemake.params.outdir
tpm = snakemake.params.tpm
sno_file = snakemake.input.sno
eftud2_file = snakemake.input.eftud2
prpf8_file = snakemake.input.prpf8

SPLICING_EVENTS = ['SE.MATS.JC.txt','A5SS.MATS.JC.txt','A3SS.MATS.JC.txt','MXE.MATS.JC.txt']


def filter_by_threshold(df):
    # Filter rMATS output by FDR and IncLevelDifference
    filtered_df = df[df['FDR']<0.05]
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
        df = df[['chr','upstreamES','downstreamEE','geneSymbol','score','strand','GeneID','FDR','IncLevelDifference']]
        df.rename(columns={"upstreamES": "start", "downstreamEE": "end"},inplace=True)
    else: # A3SS or A5SS
        pos_df = df[df['strand']=='+']
        neg_df = df[df['strand']=='-']
        if type == 'A3SS':
            pos_df = pos_df[['chr','flankingES','shortEE','geneSymbol','score','strand','GeneID','FDR','IncLevelDifference']]
            pos_df.rename(columns={"flankingES": "start", "shortEE": "end"},inplace=True)
            neg_df = neg_df[['chr','longExonStart_0base','flankingEE','geneSymbol','score','strand','GeneID','FDR','IncLevelDifference']]
            neg_df.rename(columns={"longExonStart_0base": "start", "flankingEE": "end"},inplace=True)
        else: # A5SS
            pos_df = pos_df[['chr','longExonStart_0base','flankingEE','geneSymbol','score','strand','GeneID','FDR','IncLevelDifference']]
            pos_df.rename(columns={"longExonStart_0base": "start", "flankingEE": "end"},inplace=True)
            neg_df = neg_df[['chr','flankingES','shortEE','geneSymbol','score','strand','GeneID','FDR','IncLevelDifference']]
            neg_df.rename(columns={"flankingES": "start", "shortEE": "end"},inplace=True)
        df = pd.concat([pos_df,neg_df],ignore_index=True)
    return df

def spliced_and_bound(splicing_df,event_type,binding_df):
    # bedtools intersect splicing events with binding interactions
    if 'start' not in splicing_df.columns:
        splicing_df = parse_splicing(splicing_df,event_type)
    splicing_bed = BedTool.from_dataframe(splicing_df)
    binding_bed = BedTool.from_dataframe(binding_df)
    intersection_df = splicing_bed.intersect(binding_bed, s=True,u=True).to_dataframe()
    intersection_df.columns = splicing_df.columns
    intersection_df = intersection_df.sort_values(by=['IncLevelDifference','FDR'],key=abs,ascending=[False,True]).drop_duplicates(subset=['start','end'],keep='last')
    return intersection_df


# Kallisto TPM matrix
tpm_df = pd.read_csv(tpm,sep='\t')
# snoRNA binding interactions
sno_df = pd.read_csv(sno_file,sep='\t')
# RBP binding interactions
eftud2_df = pd.read_csv(eftud2_file,sep='\t')
prpf8_df = pd.read_csv(prpf8_file,sep='\t')


# For each splicing event type, FILTER and find events bound by snoRNA and RBPs
for event in SPLICING_EVENTS:
    print(event)
    file = os.path.join(raw_dir, event)
    event_type = event.split(".")[0]
    df = pd.read_csv(file, sep='\t')
    filtered = filter_by_threshold(df)
    filtered = filter_by_tpm(filtered,tpm_df)

    # Bound by snoRNA
    sno_bound = spliced_and_bound(filtered,event_type,sno_df)
    sno_bound.to_csv(os.path.join(out_dir,event_type+'_SNORD22.tsv'),sep='\t',index=None)
    # Bound by snoRNA and EFTUD2
    sno_eftud2_bound = spliced_and_bound(sno_bound,event_type,eftud2_df)
    sno_eftud2_bound.to_csv(os.path.join(out_dir,event_type+'_SNORD22_EFTUD2.tsv'),sep='\t',index=None)
    # Bound by snoRNA and PRPF8
    sno_prpf8_bound = spliced_and_bound(sno_bound,event_type,prpf8_df)
    sno_prpf8_bound.to_csv(os.path.join(out_dir,event_type+'_SNORD22_PRPF8.tsv'),sep='\t',index=None)
    # Bound by snoRNA, EFTUD2 and PRPF8
    sno_eftud2_prpf8_bound = spliced_and_bound(sno_eftud2_bound,event_type,prpf8_df)
    sno_eftud2_prpf8_bound.to_csv(os.path.join(out_dir,event_type+'_SNORD22_EFTUD2_PRPF8.tsv'),sep='\t',index=None)
    print('Done')