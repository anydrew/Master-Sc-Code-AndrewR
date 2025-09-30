# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 17:02:18 2025

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load data
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)

# Step 1: Filter to PBMC and glioblastoma
myfilter = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_filtered = adata[myfilter].copy()

# Step 2: Preprocessing
# sc.pp.normalize_total(adata_filtered, target_sum=1e4)
# sc.pp.log1p(adata_filtered)
sc.pp.highly_variable_genes(adata_filtered, n_top_genes=2000, subset=True)
sc.pp.scale(adata_filtered, max_value=10)
sc.pp.pca(adata_filtered, n_comps=20)

# Step 3: Neighbors + clustering
sc.pp.neighbors(adata_filtered, n_neighbors=10, use_rep='X_pca')
sc.tl.umap(adata_filtered)

# Step 4: Run Leiden clustering with resolution adjusted to get 2 clusters
sc.tl.leiden(adata_filtered, resolution=0.1)  # Lower resolution → fewer clusters

# Step 5: Plot UMAP colored by cluster
sc.pl.umap(
    adata_filtered,
    color='leiden',
    title='UMAP: Unbiased Clustering of PBMC and Glioblastoma Cells (2 clusters)',
    palette='Set1'
)
