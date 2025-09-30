# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:02:27 2024

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)

# Preprocess the data (optional: filter, normalize, log transform, etc.)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find differentially expressed genes
sc.tl.rank_genes_groups(adata, 'source', groups=['glioblastoma'], reference='PBMC', method='wilcoxon')

# Extract results
deg_results = pd.DataFrame({
    'gene': adata.uns['rank_genes_groups']['names']['glioblastoma'],
    'log2_fold_change': adata.uns['rank_genes_groups']['logfoldchanges']['glioblastoma'],
    'pvals': adata.uns['rank_genes_groups']['pvals']['glioblastoma'],
    'pvals_adj': adata.uns['rank_genes_groups']['pvals_adj']['glioblastoma']
})

# Calculate -log10 p-values
deg_results['-log10_pvals_adj'] = -np.log10(deg_results['pvals_adj'])

# Plot the volcano plot
plt.figure(figsize=(10, 8))
sns.scatterplot(x='log2_fold_change', y='-log10_pvals_adj', data=deg_results, edgecolor=None)

# Highlight significant genes (adjust thresholds as necessary)
significant_genes = deg_results[(deg_results['pvals_adj'] < 0.05) & (abs(deg_results['log2_fold_change']) > 1)]
plt.scatter(significant_genes['log2_fold_change'], significant_genes['-log10_pvals_adj'], color='red')

# Add labels
for i, row in significant_genes.iterrows():
    plt.text(row['log2_fold_change'], row['-log10_pvals_adj'], row['gene'], fontsize=8)

plt.title('Volcano Plot of Differential Gene Expression: Glioblastoma vs PBMC')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-Value')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.axvline(x=1, color='black', linestyle='--')
plt.axvline(x=-1, color='black', linestyle='--')
plt.grid(True)
plt.show()
