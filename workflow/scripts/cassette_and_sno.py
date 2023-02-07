#!/usr/bin/python3

import pandas as pd
from pybedtools import BedTool
import os

sno = snakemake.params.snoRNA[0]
sno_interactions_dir = snakemake.params.sno_interactions
splicing_dir = snakemake.params.splicing
cassette_simp = snakemake.output.cassette
intersect_file = snakemake.output.intersect

def df_to_tsv(cassette,pos,f):

    # Reformat parsed data to produce a simplified version of cassette.tsv
    new_df = pd.DataFrame(columns=['gene_id','gene_name','seqid','strand','event_id',
        'exon1_start','exon1_end','exon2_start','exon2_end','exon3_start','exon3_end'])

    cassette_simp = cassette.drop_duplicates(subset=['event_id'])
    
    for e in pos.keys():

        # Remove duplicates and numerically sort start and end positions of exons
        pos[e] = sorted(list(set(pos[e])))

        # Grep info from cassette_simp
        event_info = cassette_simp[cassette_simp['event_id'] == e]
        gene_id = event_info.iloc[0]['gene_id']
        gene_name = event_info.iloc[0]['gene_name']
        seq_id = event_info.iloc[0]['seqid']
        strand = event_info.iloc[0]['strand']

        new_event = pd.DataFrame({"gene_id":gene_id,"gene_name":gene_name,'seqid':seq_id,'strand':strand,
            'event_id':[e],'exon1_start':[pos[e][0]],'exon1_end':[pos[e][1]],'exon2_start':[pos[e][2]],
            'exon2_end':[pos[e][3]],'exon3_start':[pos[e][4]],'exon3_end':[pos[e][5]]})
        new_df = pd.concat([new_df, new_event])
    
    new_df.to_csv(f, sep='\t', index=False)

    return new_df

def parse_voila(spl_dir,f):

    # Parse cassette.tsv to extract start and end positions of exons
    cassette_file = os.path.join(spl_dir, 'cassette.tsv')
    cassette_df = pd.read_csv(cassette_file, sep='\t',
        usecols=['gene_id','gene_name','seqid','strand','event_id','spliced_with_coord'])

    # Create a dictionary of lists
    # KEY: event_id  VALUE: list of start and end positions of exons
    pos_dic = {}

    for i in range(len(cassette_df)):
        
        event_id = cassette_df.iloc[i]['event_id']
        coord = cassette_df.iloc[i]['spliced_with_coord'].split("-") # list
        coord = [int(j) for j in coord] # convert strings to ints

        if event_id in pos_dic:
            pos_dic[event_id] = pos_dic[event_id] + coord
        else:
            pos_dic[event_id] = coord
    
    # Save parsed data to output files
    new_df = df_to_tsv(cassette_df,pos_dic,f)

    return new_df

def parse_sno(sno_dir,sno_id):

    # Reformat snoRNA interaction bed file to run bedtools intersect
    sno_file = os.path.join(sno_dir, sno_id+'.bed')
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

def prep_bed(df):
    
    # Reformat df for bedtools intersect

    # Add 'chr' to seqid if not already present
    df = df.astype({'seqid':'string'})
    if 'chr' not in df.iloc[0]['seqid']:
        df['seqid'] = 'chr' + df['seqid']

    df.drop(columns=['gene_id','event_id','exon1_end','exon2_start','exon2_end','exon3_start'],inplace=True)
    df['count']=3 # arbitrary value to meet bed format
    df = df[['seqid', 'exon1_start', 'exon3_end', 'gene_name', 'count', 'strand']]
    df.drop_duplicates(inplace=True)

    return df

def bedtools_intersect(cassette_df, sno_df, outfile):

    # Find cassette exons that are bound by the sno
    cassette_bed = BedTool.from_dataframe(cassette_df)
    sno_bed = BedTool.from_dataframe(sno_df)

    df_intersect = cassette_bed.intersect(sno_bed, s=True).to_dataframe()

    if len(df_intersect) > 0:
        df_intersect.drop_duplicates(inplace=True)
        df_intersect.columns = ['seqid','start','end','gene_name','count','strand']
        df_intersect.to_csv(outfile, sep='\t', index=False)
    else: # no overlaps found
        f = open(outfile, "a")
        f.write("No overlaps found!")
        f.close()

    return

# Call functions
sno_df = parse_sno(sno_interactions_dir, sno)
spl_df = parse_voila(splicing_dir, cassette_simp)
bedtools_intersect(prep_bed(spl_df), sno_df, intersect_file)