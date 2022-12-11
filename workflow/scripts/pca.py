#!/usr/bin/python3

### Adapted from Étienne Fafard-Couture's PCA script

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re

df = pd.read_csv(snakemake.input.tpm, sep='\t', index_col='gene')
df = df.T

# Standardize the values (remove mean and divide by stdev)
X = StandardScaler().fit_transform(df)

# Initialize pca
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X)
principal_df = pd.DataFrame(data = principal_components, columns = ['PC1', 'PC2'])
principal_df['sample'] = df.index

# Modify labels for legend
def legend_text(tuples):
    labels = []
    for t in tuples:
        if 'NC' in t[0]:
            labels.append(t[0])
        else:
            snoKD = 'SNORD'+t[0]+'_ASO'+str(t[1])
            labels.append(snoKD)
    return labels

# Add condition and sample information to the PCA dataframe
design = pd.read_csv(snakemake.params.design, sep='\s+')
tup = design[['condition','ASO']].apply(tuple, axis=1)
principal_df['label'] = legend_text(tup)

var1, var2 = round(pca.explained_variance_ratio_[0], 4) * 100, round(pca.explained_variance_ratio_[1], 4) * 100

# Create color palette for the samples
def color_palette(labels):
    palette = []
    for l in labels:
        if 'NC' in l and 'dimgray' not in palette:
            palette.append('dimgray')
        elif 'ASO1' in l and 'blue' not in palette:
            palette.append('blue')
        elif 'ASO2' in l and 'royalblue' not in palette:
            palette.append('royalblue')
        else:
            continue
    return palette

# Create pca_plot function
def pca_plot(df, x_col, y_col, hue_col, xlabel, ylabel, title, path, **kwargs):
    
    # Creates a PCA (scatter) plot (using a x, y and hue column).
    
    sns.set_theme()
    plt.rcParams['svg.fonttype'] = 'none'

    plt.suptitle(title, fontsize=20)
    ax = sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col, palette=color_palette(df[hue_col]), edgecolor='face',
                    alpha=0.7, s=50, **kwargs)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    ax.legend(fontsize='small', loc='upper right')

    plt.savefig(path, bbox_inches='tight', dpi=600)

# Create PCA scatter plot
pca_plot(principal_df, 'PC1', 'PC2', 'label', f'PC1 ({var1}%)', f'PC2 ({var2}%)',
        'PCA plot based on scaled TPM', snakemake.output.plot)