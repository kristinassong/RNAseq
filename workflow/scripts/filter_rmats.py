#!/usr/bin/python3

import pandas as pd
import shutil
import argparse
import os
import pyranges as pr
import numpy as np


SPLICING_EVENTS = ['SE.MATS.JC.txt','A5SS.MATS.JC.txt','A3SS.MATS.JC.txt','MXE.MATS.JC.txt','RI.MATS.JC.txt']

raw_dir = snakemake.params.dir+"/raw" # rmats output dir
out_dir = snakemake.params.dir+"/filtered" # filtered rmats output dir
tpm = snakemake.input.tpm # kallisto tpm matrix
tpm_df = pd.read_csv(tpm,sep='\t')
fdr = snakemake.params.fdr
deltapsi = snakemake.params.deltapsi
rc = snakemake.params.rc
bp_low = snakemake.params.basepsi_low
bp_high = snakemake.params.basepsi_high
gtf = snakemake.params.gtf


# Get all protein coding genes in gtf
df_gtf = pr.read_gtf(gtf).df
id_biotype = df_gtf[['gene_id','gene_biotype']].drop_duplicates(ignore_index=True)
pc_genes_list = id_biotype[id_biotype['gene_biotype']=='protein_coding'].gene_id.tolist()


def filter_by_readcounts(df,rc):
    """
    Keep events with average RNA-seq read count>=threshold in both sample groups
    """
    for col in ['IJC_SAMPLE_1','SJC_SAMPLE_1','IJC_SAMPLE_2','SJC_SAMPLE_2']:
        df[col+'_sum'] = df[col].str.split(',').apply(lambda x: np.sum([int(i) for i in x]))
    df['avg_read_count_SAMPLE_1'] = (df['IJC_SAMPLE_1_sum'] + df['SJC_SAMPLE_1_sum'])/len(df.loc[0,'IJC_SAMPLE_1'].split(','))
    df['avg_read_count_SAMPLE_2'] = (df['IJC_SAMPLE_2_sum'] + df['SJC_SAMPLE_2_sum'])/len(df.loc[0,'IJC_SAMPLE_2'].split(','))
    df = df[(df['avg_read_count_SAMPLE_1']>=rc) & (df['avg_read_count_SAMPLE_2']>=rc)]
    df.drop(columns=['IJC_SAMPLE_1_sum','SJC_SAMPLE_1_sum','IJC_SAMPLE_2_sum','SJC_SAMPLE_2_sum','avg_read_count_SAMPLE_1','avg_read_count_SAMPLE_2'],inplace=True)
    return df


def filter_by_basePSI(df,bp_low,bp_high):
    """
    Filter out events with average PSI <0.05 or >0.95 in both sample groups
    """
    for col in ['IncLevel1','IncLevel2']:
        df[col+'_avg'] = df[col].str.split(',').apply(lambda x: np.mean([float(i) for i in x]) if 'NA' not in x else -100)
    df = df[(df['IncLevel1_avg']>=bp_low) & (df['IncLevel1_avg']<=bp_high) & (df['IncLevel2_avg']>=bp_low) & (df['IncLevel2_avg']<=bp_high)]
    df.drop(columns=['IncLevel1_avg','IncLevel2_avg'],inplace=True)
    return df


def filter_by_tpm(rmats_df,tpm_df):
    """
    Keep genes that are expressed in at least one sample based on the kallisto TPM matrix
    """

    # Drop all rows in tpm_df where max TPM<1
    samples = tpm_df.columns.values.tolist()
    for i in ['gene','transcript','gene_name']:
        samples.remove(i)
    exp_tpm_df = tpm_df[~((tpm_df[samples]<1).all(axis=1))]

    filtered_rmats_df = pd.DataFrame(columns=rmats_df.columns)
    all_genes = set(rmats_df['GeneID'].values.tolist())
    
    for gene in all_genes:
        gene_tpm_df = exp_tpm_df[exp_tpm_df['gene']==gene]
        if len(gene_tpm_df)>0 and gene in pc_genes_list:
            gene_rmats_df = rmats_df[rmats_df['GeneID']==gene]
            filtered_rmats_df = pd.concat([filtered_rmats_df,gene_rmats_df],ignore_index=True)

    return filtered_rmats_df


"""
For each splicing event type, FILTER by multiple thresholds
"""
for event in SPLICING_EVENTS:
    file = os.path.join(raw_dir, event)
    event_type = event.split(".")[0]
    df = pd.read_csv(file, sep='\t')

    # All significant AS events
    filtered_df = filter_by_readcounts(df,rc)
    filtered_df = filter_by_basePSI(filtered_df,bp_low,bp_high)
    all_sig_df = filter_by_tpm(filtered_df,tpm_df)

    # Differential AS events
    diff_df = all_sig_df[all_sig_df['FDR']<=fdr]
    diff_df = diff_df[diff_df['IncLevelDifference'].abs()>=deltapsi]

    # Splicing events filtered by threshold and TPM
    all_sig_df.to_csv(os.path.join(out_dir,'all_sig_'+event_type+'.tsv'),sep='\t',index=None)
    diff_df.to_csv(os.path.join(out_dir,'diff_'+event_type+'.tsv'),sep='\t',index=None)


# Delete tmp directory
tmp_dir = snakemake.params.dir+"/tmp"
if os.path.isdir(tmp_dir):
    shutil.rmtree(tmp_dir)