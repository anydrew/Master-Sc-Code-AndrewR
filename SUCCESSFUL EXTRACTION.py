# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 18:32:51 2024

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd

# File path to the .h5ad file
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"

# Load the data from the .h5ad file
adata = sc.read_h5ad(file_path)

# Ensure that 'source' is in the observations (obs) of the AnnData object
if 'source' not in adata.obs.columns:
    raise ValueError("The 'source' column is not found in the AnnData object.")

# Define the subsets
subsets = ['PBMC', 'glioblastoma']

# Dictionary to hold mean expression dataframes
mean_expression_dfs = {}

# Process each subset
for subset in subsets:
    # Filter the AnnData object for the current subset
    subset_data = adata[adata.obs['source'] == subset]

    # Calculate mean expression for each gene
    mean_expression = subset_data.X.mean(axis=0)

    # Create a DataFrame for the mean expression
    mean_expression_df = pd.DataFrame(mean_expression.A1, index=subset_data.var_names, columns=['Mean Expression'])

    # Store the DataFrame in the dictionary
    mean_expression_dfs[subset] = mean_expression_df

# Create an Excel writer object
excel_path = r"C:\Users\aroche6\abc\TEST250924.xlsx"
with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    # Write each subset's mean expression DataFrame to a separate sheet
    for subset, df in mean_expression_dfs.items():
        df.to_excel(writer, sheet_name=subset)

print(f"Mean expression values have been saved to {excel_path}")