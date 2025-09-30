# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 00:28:57 2025

@author: AROCHE6
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 21:29:41 2025

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

# Dictionary to hold fraction of cells dataframes
fraction_of_cells_dfs = {}

# Process each subset
for subset in subsets:
    # Filter the AnnData object for the current subset
    subset_data = adata[adata.obs['source'] == subset]
    
    # Calculate fraction of cells expressing each gene
    fraction_of_cells = (subset_data.X > 0).mean(axis=0)
    
    # Create a DataFrame for the fraction of cells
    fraction_of_cells_df = pd.DataFrame(
        fraction_of_cells.A1, index=subset_data.var_names, columns=['Fraction of Cells']
    )
    
    # Store the DataFrame in the dictionary
    fraction_of_cells_dfs[subset] = fraction_of_cells_df

# Create an Excel writer object
excel_path = r"C:\Users\aroche6\abc\FractionOfCells_AllGenes.xlsx"
with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    # Write each subset's fraction of cells DataFrame to a separate sheet
    for subset, df in fraction_of_cells_dfs.items():
        df.to_excel(writer, sheet_name=subset)

print(f"Fraction of cells values for all genes have been saved to {excel_path}")
