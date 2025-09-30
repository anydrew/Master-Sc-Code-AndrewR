# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:13:10 2024

@author: AROCHE6
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:08:27 2024

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd

# File path to the .h5ad file
file_path =  r"C:\Users\aroche6\abc\all_nk_cells.h5ad"

# Load the data from the .h5ad file
adata = sc.read_h5ad(file_path)

#myfilter4 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')


# Ensure that 'source' is in the observations (obs) of the AnnData object
if 'subset_source' not in adata.obs.columns:
    raise ValueError("The 'subset_source' column is not found in the AnnData object.")

# Define the subsets
subsets = ['KIR+, PBMC', 'CD57+, PBMC', 'NKG2A+, PBMC', 'CD56bright, PBMC', 'Adaptive, PBMC', 'CD56bright, glioblastoma_tumor', 'CD56dim, glioblastoma_tumor']

# Dictionary to hold mean expression dataframes
mean_expression_dfs = {}

# Process each subset
for subset in subsets:
    # Filter the AnnData object for the current subset
    subset_data = adata[adata.obs['subset_source'] == subset]

    # Calculate mean expression for each gene
    mean_expression = subset_data.X.mean(axis=0)

    # Create a DataFrame for the mean expression
    mean_expression_df = pd.DataFrame(mean_expression.A1, index=subset_data.var_names, columns=['Mean Expression'])
    
    # Store the DataFrame in the dictionary
    mean_expression_dfs[subset] = mean_expression_df

# Create an Excel writer object
excel_path = r"C:\Users\aroche6\abc\qwert.xlsx"
with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    # Write each subset's mean expression DataFrame to a separate sheet
    for subset, df in mean_expression_dfs.items():
        df.to_excel(writer, sheet_name=subset)

print(f"Mean expression values have been saved to {excel_path}")
