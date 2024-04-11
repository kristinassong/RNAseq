#!/usr/bin/env python

"""
Perform GO analysis on DE genes and visualize through a bar chart.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gprofiler import GProfiler


# Get list of genes
genes = list(pd.read_csv(snakemake.input.genes, sep='\t')['gene'])

# Run gprofiler for enrichment analysis of functional (GO and other) terms
### NEED INTERNET CONNECTION!
gp = GProfiler(return_dataframe=True)
go_results = gp.profile(organism='hsapiens', query=genes)

# Log-transform p-values
go_results['negative_log10_of_p_value'] = np.negative(np.log10(go_results['p_value']))

# Extract GO terms of interest
go_mf = go_results[go_results['source']=='GO:MF'].sort_values(by=['negative_log10_of_p_value'])
go_cc = go_results[go_results['source']=='GO:CC'].sort_values(by=['negative_log10_of_p_value'])
go_bp = go_results[go_results['source']=='GO:BP'].sort_values(by=['negative_log10_of_p_value'])

go_terms = {"GO:BP":go_bp,"GO:MF":go_mf,"GO:CC":go_cc}
colours = {"GO:BP":"indianred","GO:MF":"limegreen","GO:CC":"steelblue"}

go_concat = pd.DataFrame()
colours_present = {}

# Only keep GO term types that contain at least one GO term of that type
for term in go_terms.keys():
    if len(go_terms[term]) > 0:
        go_concat = pd.concat([go_concat,go_terms[term]], ignore_index=True)
        colours_present.update({term:colours[term]})

# Represent GO results with a bar chart
plt.figure(figsize=(9,4))

bars = pd.DataFrame({'source':go_concat.source.values.tolist(),'p_val':go_concat.negative_log10_of_p_value.values.tolist()},
    index=go_concat.name.values.tolist())

bars['p_val'].plot(kind="barh",color=bars['source'].replace(colours_present),width=0.7)
plt.title('GO Analysis of genes')
plt.xlabel('$-log_{10}$(p-value)')

labels = list(colours_present.keys())
handles = [plt.Rectangle((0,0),1,1, color=colours_present[label]) for label in labels]
plt.legend(handles,labels,loc='lower right')

plt.tight_layout()
plt.savefig(snakemake.output.bar_chart)