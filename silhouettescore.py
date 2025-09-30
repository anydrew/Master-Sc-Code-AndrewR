# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 20:55:00 2025

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
# Filter to PBMC and glioblastoma
# # Define the groups
# dim_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
# bright_source = ['CD56bright_PBMC']

# # Filter cells belonging to either dim_sources or bright_source
# myfilter6 = adata.obs['subset_source'].isin(dim_sources + bright_source)
# adata_filtered = adata[myfilter6].copy()

# # Label groups: 'Dim' for dim_sources and 'Bright' for bright_source
# adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
#     lambda x: 'Dim' if x in dim_sources else 'Bright'
# )

pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Label PBMC vs Glioblastoma
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Make sure UMAP coordinates exist in filtered data
adata_filtered.obsm['X_umap'] = adata[myfilter7].obsm['X_umap']

# Prepare data for silhouette calculation
X_umap = adata_filtered.obsm['X_umap']
labels = adata_filtered.obs['source_group'].values

# Convert labels to numeric for silhouette_score function
label_numeric = np.array([0 if l == 'Dim' else 1 for l in labels])

# Calculate Silhouette Score using UMAP coordinates and group labels
sil_score = silhouette_score(X_umap, label_numeric)

print(f"Silhouette Score for Dim vs Bright groups on UMAP embeddings: {sil_score:.4f}")
