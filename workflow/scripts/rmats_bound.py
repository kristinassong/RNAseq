#!/usr/bin/python3

import pandas as pd
import os
from pybedtools import BedTool


SPLICING_EVENTS = ['SE.tsv','A5SS.tsv','A3SS.tsv','MXE.tsv','RI.tsv']

rmats_dir = "/home/kris98/scratch/RNAseq/results/rmats/NC5-22/filtered"
sno_interaction = "/home/kris98/scratch/snoFlak/results/interactions/snoGloBe/ENSG00000277194.bed"
sno_df = pd.read_csv(sno_interaction, sep='\t')
sno_bed = BedTool.from_dataframe(sno_df)


def pad_splicing_events(df):
    """
    Pad splicing events by +- 100 nts.
    """
    df['score'] = 3
    df = df[['chr','upstreamES','downstreamEE','GeneID','score','strand']]
    df['upstreamES'] -= 100
    df['downstreamEE'] += 100
    df['chr'] = df['chr'].map(lambda x: x.lstrip('chr')) # remove 'chr' string
    return df


def get_full_entry(bed_df,org_df):
    """
    Get full information on selected splicing events.
    """
    final_df = pd.DataFrame(columns=org_df.columns)

    for i in range(len(bed_df)):
        id = bed_df.iloc[i]['ID']
        start = bed_df.iloc[i]['start']
        end = bed_df.iloc[i]['end']
        # go back to values before padding
        start += 100
        end -= 100

        row = org_df[(org_df['ID']==id) & (org_df['start']==start) & (org_df['end']==end)]
        final_df = pd.concat([final_df,row],ignore_index=True)

    return final_df


"""
For each splicing event type, get splicing events bound by a snoRNA
"""
for event in SPLICING_EVENTS:
    file = os.path.join(rmats_dir, event)
    event_type = event.split(".")[0]
    print(event_type)
    df = pd.read_csv(file, sep='\t')
    padded_df = pad_splicing_events(df)
    padded_bed = BedTool.from_dataframe(padded_df)

    intersection_df = padded_bed.intersect(sno_bed,s=True,u=True).to_dataframe(disable_auto_names=True, names=['chr','start','end','ID','score','strand'])
    sno_name = os.path.basename(sno_interaction).split(".")[0]
    outfile = os.path.join(rmats_dir, event_type+"_"+sno_name+".tsv")
    get_full_entry(intersection_df, df).to_csv(outfile,sep='\t',index=None)