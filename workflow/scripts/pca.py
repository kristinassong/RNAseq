#!/usr/bin/python3

### Adapted from Étienne Fafard-Couture's PCA script

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


df = pd.read_csv(snakemake.input.tpm, sep='\t', index_col='gene')
df = df.drop(columns=['transcript', 'gene_name'])
df = df.T

# Standardize the values (remove mean and divide by stdev)
X = StandardScaler().fit_transform(df)

# Initialize pca
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X)
principal_df = pd.DataFrame(data = principal_components, columns = ['PC1', 'PC2'])
principal_df['sample'] = df.index

# Modify labels for legend
def legend_text(design,principal_df):
    labels = []
    for s in principal_df['sample'].values.tolist():
        labels.append(design[design['sample']==s].iloc[0]['condition'])
    return labels

# Add condition and sample information to the PCA dataframe
design = pd.read_csv(snakemake.params.design, sep='\s+')
principal_df['label'] = legend_text(design,principal_df)

var1, var2 = round(pca.explained_variance_ratio_[0], 4) * 100, round(pca.explained_variance_ratio_[1], 4) * 100

# Create pca_plot function
def pca_plot(df, x_col, y_col, hue_col, xlabel, ylabel, title, path, **kwargs):
    
    # Creates a PCA (scatter) plot (using a x, y and hue column).
    
    plt.figure(figsize=(5.5,4))
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams["legend.loc"] = 'upper right'

    plt.suptitle(title, fontsize=16)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col, edgecolor='face',
                    alpha=0.7, s=50, **kwargs)

    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)
    plt.legend(fontsize='medium')

    plt.savefig(path, bbox_inches='tight', dpi=600)

# Create PCA scatter plot
pca_plot(principal_df, 'PC1', 'PC2', 'label', f'PC1 ({var1:.2f}%)', f'PC2 ({var2:.2f}%)',
        'PCA plot based on scaled TPM', snakemake.output.plot)

principal_df.to_csv(snakemake.output.tsv, sep='\t', index=False)