# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 04:17:01 2024

@author: AROCHE6
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text  # Correct import

# Load the Excel file
file_path = r"C:\Users\aroche6\abc\LOG DFS PBMC VS TUMOR.xlsx"
deg_results = pd.read_excel(file_path)

# Ensure the data has the necessary columns
required_columns = ['gene', 'log2_fold_change', 'pvals', 'pvals_adj']
if not all(col in deg_results.columns for col in required_columns):
    raise ValueError(f"Input file must contain the following columns: {required_columns}")

# Calculate -log10 p-values
deg_results['-log10_pvals_adj'] = -np.log10(deg_results['pvals_adj'])

# Plot the volcano plot for all genes
plt.figure(figsize=(15, 10))
sns.scatterplot(x='log2_fold_change', y='-log10_pvals_adj', data=deg_results, edgecolor=None)

# Highlight significant genes (adjust thresholds as necessary)
significant_genes = deg_results[(deg_results['pvals_adj'] < 0.05) & (abs(deg_results['log2_fold_change']) > 1)]
plt.scatter(significant_genes['log2_fold_change'], significant_genes['-log10_pvals_adj'], color='red')

# Add labels to significant genes, ensuring they don't overlap
texts = []
for i, row in significant_genes.iterrows():
    texts.append(plt.text(row['log2_fold_change'], row['-log10_pvals_adj'], row['gene'], fontsize=8))

# Adjust text to prevent overlap
adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

plt.title('Volcano Plot of Differential Gene Expression: Glioblastoma vs PBMC')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-Value')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.axvline(x=1, color='black', linestyle='--')
plt.axvline(x=-1, color='black', linestyle='--')
plt.grid(True)
plt.show()
output_file_path = r"C:\Users\aroche6\abc\FullVolcanoPlot.xlsx"
significant_genes.to_excel(output_file_path, index=False)