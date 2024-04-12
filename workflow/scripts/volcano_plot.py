#!/usr/bin/python3

### Adapted from the script written by Danny Bergeron

"""
Create a volcano plot to show DE genes identified by DESeq2.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from gtfparse import read_gtf
import numpy as np


# Configure comparison, pval_threshold and colors
comp1, comp2 = str(snakemake.wildcards.comp).split('-')

# For ASO KD
"""
if 'NC' not in comp1:
        comp1 = 'SNORD'+comp1+' KD'
if 'NC' not in comp2:
        comp2 = 'SNORD'+comp2+' KD'
"""

pval_threshold = snakemake.params.pval_threshold
colors = {'log2FC < -1 & padj < '+ str(pval_threshold): 'cornflowerblue','log2FC > 1 & padj < '+ str(pval_threshold): 'firebrick','n.s.': 'grey'}

# Load DE df
df = pd.read_csv(snakemake.input.DE_output)

# Rename columns
df = df.set_axis(['gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'], axis=1)

# Drop genes/transcripts with NaN in log2FoldChange, pvalue and/or padj
df = df.dropna(subset=['log2FoldChange', 'pvalue', 'padj'])

# Filter genes by biotype
# Only keep protein coding genes
gtf = snakemake.params.gtf
df_gtf = read_gtf(gtf).to_pandas()
df_gtf = df_gtf[df_gtf['gene_biotype']=='protein_coding']
df = df[df.gene.isin(df_gtf.gene_id)]
print(len(df))

# Create -log10padj column
df['-log10padj'] = -np.log10(df['padj'])

# Create hue column for significant points (|log2FC|>1 & padj<0.05)
df.loc[(df['padj'] < pval_threshold) & (df['log2FoldChange'] > 1),
        'sig.'] = 'log2FC > 1 & padj < '+str(pval_threshold)
df.loc[(df['padj'] < pval_threshold) & (df['log2FoldChange'] < -1),
        'sig.'] = 'log2FC < -1 & padj < '+str(pval_threshold)
df['sig.'] = df['sig.'].replace('nan','n.s.')
df['sig.'] = df['sig.'].fillna('n.s.')

# Extract genes that are significant
outfile_up = snakemake.output.up_genes
outfile_down = snakemake.output.down_genes
df_genes = df[~df['sig.'].str.contains('n.s.')]
df_genes = df_genes[['gene','log2FoldChange','padj']]

up = df_genes[(df_genes['log2FoldChange'] > 1) & (df_genes['padj'] < pval_threshold)]
print(len(up))
up.sort_values('log2FoldChange',inplace=True,ascending=False)
down = df_genes[df_genes['log2FoldChange'] < -1 & (df_genes['padj'] < pval_threshold)]
down.sort_values('log2FoldChange',inplace=True,ascending=False)
print(len(down))
up.to_csv(outfile_up, sep='\t', index=False)
down.to_csv(outfile_down, sep='\t', index=False)

# Create volcano function
def volcano(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            pval_threshold, **kwargs):
    """
    Creates a volcano plot (using a x, y and hue column).
    """
    #sns.set_theme()
    plt.figure(figsize=(5.5,4))
    plt.rcParams['svg.fonttype'] = 'none'

    plt.suptitle(title, fontsize=16)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    # Add threshold lines (padj)
    plt.axhline(y=-np.log10(pval_threshold), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(2), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(0.5), color='black', ls='--', lw=0.5)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    leg = plt.legend(fontsize="medium",bbox_to_anchor=(1.6, 0.4), loc='lower right', ncol=1,
        facecolor='white',framealpha=1)
    leg.get_frame().set_linewidth(0)
    plt.savefig(path, bbox_inches='tight', dpi=600)
    
    return

# Create volcano
volcano(df, 'log2FoldChange', '-log10padj', 'sig.',
        'log2(Fold Change)', '-log10(FDR-adjusted p-value)',
        f'{comp1} vs {comp2} using DESeq2',
        colors, snakemake.output.volcano, pval_threshold)