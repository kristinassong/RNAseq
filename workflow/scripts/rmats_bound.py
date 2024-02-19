#!/usr/bin/python3

import pandas as pd
import os
from pybedtools import BedTool


SPLICING_EVENTS = ['SE.tsv','A5SS.tsv','A3SS.tsv','MXE.tsv','RI.tsv']

rmats_dir = "/home/kristinasong/Desktop/git/RNAseq/results/rmats/NC5-22/filtered"
sno_interaction = "/home/kristinasong/Desktop/interactions/snoGloBe/ENSG00000277194.bed"
sno_df = pd.read_csv(sno_interaction, sep='\t')
sno_bed = BedTool.from_dataframe(sno_df)


def parse_splicing_df(df,type):
    """
    Transform splicing df into bed format
    """
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

    df = df[['chr','start','end','ID','score','strand']]

    return df


def pad_splicing_events(df):
    """
    Pad splicing events by +- 100 nts.
    """
    df['start'] -= 100
    df['end'] += 100
    df['chr'] = df['chr'].map(lambda x: x.lstrip('chr')) # remove 'chr' string
    return df


def get_full_entry(bed_df,org_df):
    """
    Get full information on selected splicing events.
    """
    final_df = pd.DataFrame(columns=org_df.columns)

    for i in range(len(bed_df)):
        id = bed_df.iloc[i]['ID']
        row = org_df[org_df['ID']==id]
        final_df = pd.concat([final_df,row],ignore_index=True)

    return final_df


"""
For each splicing event type, get splicing events bound by a snoRNA
"""
for event in SPLICING_EVENTS:
    file = os.path.join(rmats_dir, event)
    event_type = event.split(".")[0]
    print(event_type+"...")
    df_org = pd.read_csv(file, sep='\t')
    df_parsed = parse_splicing_df(df_org,event_type)
    padded_df = pad_splicing_events(df_parsed)
    padded_bed = BedTool.from_dataframe(padded_df)

    intersection_df = padded_bed.intersect(sno_bed,s=True,u=True).to_dataframe(disable_auto_names=True, names=['chr','start','end','ID','score','strand'])
    sno_name = os.path.basename(sno_interaction).split(".")[0]
    outfile = os.path.join(rmats_dir, event_type+"_"+sno_name+".tsv")
    get_full_entry(intersection_df, df_org).to_csv(outfile,sep='\t',index=None)