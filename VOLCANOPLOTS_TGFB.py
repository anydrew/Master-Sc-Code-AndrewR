# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:52:49 2024

@author: AROCHE6
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the Excel file
file_path = r"C:\Users\aroche6\abc\LOG DFS PBMC VS TUMOR.xlsx"
deg_results = pd.read_excel(file_path)

# Ensure the data has the necessary columns
required_columns = ['gene', 'log2_fold_change', 'pvals', 'pvals_adj']
if not all(col in deg_results.columns for col in required_columns):
    raise ValueError(f"Input file must contain the following columns: {required_columns}")

# Calculate -log10 p-values
deg_results['-log10_pvals_adj'] = -np.log10(deg_results['pvals_adj'])

# Define the TGFB genes of interest
conor_genes = ["HIF1-a", "SOD1", "SOD2", "GSH", "GSSH", "TRX-1", "NOX1", "NOX2", "CD38"]

# Filter the results to include only the TGFB genes
tgfb_results = deg_results[deg_results['gene'].isin(conor_genes)]

# Plot the volcano plot for TGFB genes
plt.figure(figsize=(10, 8))
sns.scatterplot(x='log2_fold_change', y='-log10_pvals_adj', data=tgfb_results, edgecolor=None)

# Highlight significant genes (adjust thresholds as necessary)
significant_genes = tgfb_results[(tgfb_results['pvals_adj'] < 0.05) & (abs(tgfb_results['log2_fold_change']) > 1)]
plt.scatter(significant_genes['log2_fold_change'], significant_genes['-log10_pvals_adj'], color='red')

# # Add labels to significant genes
# for i, row in significant_genes.iterrows():
#     plt.text(row['log2_fold_change'], row['-log10_pvals_adj'], row['gene'], fontsize=8)

plt.title('Volcano Plot of Differential Gene Expression for Conor Genes: HC-PBMC-NK vs GBM-NK')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-Value')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.axvline(x=1, color='black', linestyle='--')
plt.axvline(x=-1, color='black', linestyle='--')
plt.grid(True)
plt.show()
