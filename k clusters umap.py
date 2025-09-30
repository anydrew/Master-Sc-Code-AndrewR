# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 17:44:59 2025

@author: AROCHE6
"""
import scanpy as sc
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# Step 1: Load the AnnData object
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)

# Step 2: Filter to PBMC and glioblastoma
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Step 3: (Optional) Subset to highly variable genes for faster computation
# sc.pp.highly_variable_genes(adata_filtered, n_top_genes=1000)
# adata_filtered = adata_filtered[:, adata_filtered.var['highly_variable']].copy()

# Step 4: Run KMeans directly on raw expression matrix
# Run KMeans on UMAP coordinates (2D)
kmeans = KMeans(n_clusters=8, random_state=0).fit(adata_filtered.obsm['X_umap'])
adata_filtered.obs['kmeans_umap'] = kmeans.labels_.astype(str)

# Plot
sc.pl.umap(
    adata_filtered,
    color='kmeans_umap',
    title='KMeans Clustering (UMAP Coordinates, 2 Clusters)',
    legend_loc='on data'
)


# # Step 5: Use existing UMAP from full dataset (if available)
# if 'X_umap' in adata.obsm:
#     adata_filtered.obsm['X_umap'] = adata.obsm['X_umap'][myfilter]
# else:
#     # If UMAP not present, compute it minimally
#     sc.pp.neighbors(adata_filtered)
#     sc.tl.umap(adata_filtered)

# # Step 6: Plot KMeans clusters on UMAP
# sc.pl.umap(
#     adata_filtered,
#     color='kmeans_raw',
#     title='KMeans Clustering (Raw Expression, 2 Clusters)',
#     legend_loc='on data'
# )

# # # Run KMeans on UMAP coordinates (2D)
# kmeans = KMeans(n_clusters=2, random_state=0).fit(adata_filtered.obsm['X_umap'])
# adata_filtered.obs['kmeans_umap'] = kmeans.labels_.astype(str)

# # Plot
# sc.pl.umap(
#     adata_filtered,
#     color='kmeans_umap',
#     title='KMeans Clustering (UMAP Coordinates, 2 Clusters)',
#     legend_loc='on data'
# )
