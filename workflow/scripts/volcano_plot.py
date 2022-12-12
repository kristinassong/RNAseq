#!/usr/bin/python3

### Adapted from the script written by Danny Bergeron

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Get list of genes after TPM filtering
genes_file = snakemake.input.filtered_genes
genes = pd.read_csv(genes_file, sep='\t', usecols=['gene'])

# Configure comparison, pval_threshold and colors
comp1, comp2 = str(snakemake.wildcards.comp).split('-')

if 'NC' not in comp1:
        comp1 = 'SNORD'+comp1
if 'NC' not in comp2:
        comp2 = 'SNORD'+comp2

#DE_tool, quantifier = str(snakemake.wildcards.DE_tool), str(snakemake.wildcards.quant)
pval_threshold = snakemake.params.pval_threshold
colors = {'|log2FC| > 1 & padj < '+str(pval_threshold): 'blue',
            'n.s.': 'grey'}

# Load DE df
df = pd.read_csv(snakemake.input.DE_output)

# Rename columns
df.set_axis(['gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'], axis=1, inplace=True)

# Drop genes/transcripts with NaN in log2FoldChange, pvalue and/or padj
df = df.dropna(subset=['log2FoldChange', 'pvalue', 'padj'])

# Drop genes that are not in the TPM filtered list
df = df[df.gene.isin(genes.gene)] 

# Create -log10padj column
df['-log10padj'] = -np.log10(df['padj'])

# Create hue column for significative points (|log2FC| > 1 & padj<0.05)
df.loc[(df['padj'] < pval_threshold) & (np.abs(df['log2FoldChange']) > 1),
        'sig.'] = '|log2FC| > 1 & padj < '+str(pval_threshold)
df['sig.'] = df['sig.'].fillna('n.s.')

# Create volcano function
def volcano(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            pval_threshold, **kwargs):
    """
    Creates a volcano plot (using a x, y and hue column).
    """
    sns.set_theme()
    plt.rcParams['svg.fonttype'] = 'none'

    plt.suptitle(title, fontsize=25)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col,
                    palette=color_dict, edgecolor='face',
                    s=25, **kwargs)

    # Add threshold lines (padj)
    plt.axhline(y=-np.log10(pval_threshold), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(2), color='black', ls='--', lw=0.5)
    plt.axvline(x=np.log2(0.5), color='black', ls='--', lw=0.5)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    leg = plt.legend(fontsize="medium",bbox_to_anchor=(1.6, 0.4), loc='lower right', ncol=1,
        facecolor='white',framealpha=1)
    leg.get_frame().set_linewidth(0)
    plt.savefig(path, bbox_inches='tight', dpi=600)
    
    return

# Create volcano
volcano(df, 'log2FoldChange', '-log10padj', 'sig.',
        'log2(Fold Change)', '-log10(FDR-adjusted p-value)',
        f'Comparison between {comp1} and {comp2} using DESeq2 after Kallisto',
        colors, snakemake.output.volcano, pval_threshold)