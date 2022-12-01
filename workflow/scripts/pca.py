#!/usr/bin/python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re

tpm = snakemake.input.tpm

df = pd.read_csv(tpm, sep='\t', index_col='gene')
df = df.T

# Standardize the values (remove mean and divide by stdev)
X = StandardScaler().fit_transform(df)

# Initialize pca
pca = PCA(n_components=2)
principal_components = pca.fit_transform(X)
principal_df = pd.DataFrame(data = principal_components, columns = ['PC1', 'PC2'])
var1, var2 = round(pca.explained_variance_ratio_[0], 4) * 100, round(pca.explained_variance_ratio_[1], 4) * 100

# Create pca_plot function
def pca_plot(df, x_col, y_col, hue_col, xlabel, ylabel, title, path, **kwargs):
    
    #Creates a PCA (scatter) plot (using a x, y and hue column).
    
    sns.set_theme()
    plt.rcParams['svg.fonttype'] = 'none'

    plt.suptitle(title, fontsize=25)
    sns.scatterplot(data=df, x=x_col, y=y_col, hue=hue_col, edgecolor='face',
                    alpha=0.7, s=50, **kwargs)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    plt.savefig(path, bbox_inches='tight', dpi=600)

# Create PCA scatter plot
pca_plot(principal_df, 'PC1', 'PC2', 'sample', f'PC1 ({var1}%)', f'PC2 ({var2}%)',
        'PCA plot based on scaled TPM', snakemake.output.plot)