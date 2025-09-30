# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 17:15:23 2025

@author: AROCHE6
"""

import scanpy as sc
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
# Step 1: Filter to PBMC and Glioblastoma
myfilter = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_filtered = adata[myfilter].copy()

# Step 2: (Skip normalization and log1p — assuming it's already done)

# Step 3: Highly variable gene selection (optional but recommended if not done yet)
# # You can comment this out if you know HVGs are already selected.
# sc.pp.highly_variable_genes(adata_filtered, n_top_genes=2000, subset=True)

# # Step 4: Scale & PCA
# sc.pp.scale(adata_filtered, max_value=10)
# sc.pp.pca(adata_filtered, n_comps=20)

# Step 5: Neighbors and UMAP for visualization
sc.pp.neighbors(adata_filtered, n_neighbors=10, use_rep='X_pca')
sc.tl.umap(adata_filtered)

# Step 6: Leiden clustering with low resolution to force 2 clusters
sc.tl.leiden(adata_filtered, resolution=0.00000000000000000000000000000000000001)  # adjust resolution if >2 clusters

# Step 7: Plot UMAP with clusters
sc.pl.umap(adata_filtered, color='leiden', title='Unbiased Clusters (Leiden)')
