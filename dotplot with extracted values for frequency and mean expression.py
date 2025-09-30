# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 23:10:46 2025

@author: AROCHE6
"""


import scanpy as sc
import pandas as pd

# File paths
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
output_file = r"C:\Users\aroche6\abc\REAL_SIGOXSTRESS_with_Fractions.xlsx"

# Load the data
adata = sc.read_h5ad(file_path)

# Filter data for PBMC and glioblastoma groups
myfilter12 = (adata.obs['source'] == 'PBMC') | (adata.obs['source'] == 'glioblastoma')

# Define genes of interest
sigglycgluc = ["APOE", "AREG", "HSPA1B", "CDK1", "BCL2", "AIF1", "HSPA1A", "EPAS1", "RUNX2", "PLK1", "SMPD3", "NCF1", "FOS", "PLCG2", "APP", "PLCB1", "ABCC1", "NR4A2", "MAP1LC3A", "AIFM2", "FER", "ATP2A2", "FOXO1", "AGAP3", "MAPK13", "SLC25A29", "PRDX1", "MSRA", "ERCC6", "GSTP1", "GSR", "DHFR", "ZFP36", "PDXK", "EZH2", "KAT2B", "MAPK8", "FANCD2", "CD38", "PPIA", "ERCC1", "TXN", "MAP2K4", "NFE2L2", "PLCG1", "CASP3", "MTHFD2", "ABL1", "JAK2", "SETX", "SLC25A12", "SIRT1", "TXNIP", "IPCEF1", "SESN1", "MTHFD1", "STK11", "MAPK9", "HM13", "FUT8", "TNFAIP3", "MAPK3", "GCH1", "JUN", "NUDT1", "NDUFS8", "CPEB2", "SLC25A32", "SLC2A3", "PRDX3", "MAPK1", "AKT1", "KDM6B", "PIK3R5", "GCLC", "PRDX2", "TGFBR3", "DUOX1", "SLC25A6"
]
dotplot_data = sc.pl.dotplot(
    adata[myfilter12], 
    sigglycgluc, 
    groupby='source', 
    title="CHEMOTAXIS_genes Geneset Expression GBM-NK vs HC-PBMC-NK",
    return_fig=True  # Retrieve the data used in the dot plot
)

# Extract mean expression (dot color) data
mean_expression_df = dotplot_data.dot_color_df

# Calculate fraction of cells (dot size) data
fraction_of_cells_df = dotplot_data.dot_size_df

# Combine data into a single DataFrame for export
combined_data = pd.concat(
    [mean_expression_df, fraction_of_cells_df], 
    keys=["Mean Expression", "Fraction of Cells"], 
    axis=1
)

# Save data to Excel
combined_data.to_excel(output_file)

print(f"Dot plot data with fractions of cells extracted and saved to {output_file}")
