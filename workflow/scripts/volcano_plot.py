#!/usr/bin/python3

### Script written by Danny Bergeron

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Configure comparison, pval_threshold and colors
comp1, comp2 = str(snakemake.wildcards.comp).split('-')
#DE_tool, quantifier = str(snakemake.wildcards.DE_tool), str(snakemake.wildcards.quant)
pval_threshold = snakemake.params.pval_threshold
colors = {'|log2FC| > 1 & padj<'+str(pval_threshold): 'lightgreen',
            'Not significative': 'grey'}

# Load DE df
df = pd.read_csv(snakemake.input.DE_output)

# Rename columns to standardize their names across DE tools (Convert EdgeR and
# Sleuth cols into the same DESeq2 cols)
std_names = {'logFC': 'log2FoldChange', 'b': 'log2FoldChange',
            'PValue': 'pvalue', 'FDR': 'padj', 'pval': 'padj'}
df = df.rename(columns=std_names)

# Drop genes/transcripts with NaN in log2FoldChange, pvalue and/or padj
df = df.dropna(subset=['log2FoldChange', 'pvalue', 'padj'])

# Create -log10padj column
df['-log10padj'] = -np.log10(df['padj'])

# Create hue column for significative points (|log2FC| > 1 & padj<0.05)
df.loc[(df['padj'] < pval_threshold) & (np.abs(df['log2FoldChange']) > 1),
        'Significative data'] = '|log2FC| > 1 & padj<'+str(pval_threshold)
df['Significative data'] = df['Significative data'].fillna('Not significative')

# Create volcano function
def volcano(df, x_col, y_col, hue_col, xlabel, ylabel, title, color_dict, path,
            pval_threshold, **kwargs):
    """
    Creates a violin plot (using a x, y and hue column).
    """
    sns.set_theme()
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["legend.loc"] = 'upper right'

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
    plt.legend(fontsize="small")
    plt.savefig(path, bbox_inches='tight', dpi=600)

# Create volcano
volcano(df, 'log2FoldChange', '-log10padj', 'Significative data',
        'log2(Fold change)', '-log10(FDR-adjusted p-value)',
        f'Comparison between {comp1} and {comp2} using DESeq2 after Kallisto',
        colors, snakemake.output.volcano, pval_threshold)