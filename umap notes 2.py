# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 15:29:29 2025

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd
import umap
import umap.plot
import matplotlib.pyplot as plt
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
# Step 1: Filter to dim and bright groups
dim_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
bright_source = ['CD56bright_PBMC']

myfilter = adata.obs['subset_source'].isin(dim_sources + bright_source)
adata_filtered = adata[myfilter].copy()

# Step 2: Label groups
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'dim' if x in dim_sources else 'bright'
)

# Step 3: PCA before UMAP
sc.pp.pca(adata_filtered, n_comps=20)

# Step 4: Train UMAP with umap-learn
reducer = umap.UMAP()
embedding = reducer.fit_transform(adata_filtered.obsm['X_pca'])

# Step 5: Plot connectivity structure with edge bundling
umap.plot.connectivity(reducer, edge_bundling='hammer')

# Optional Step 6: Plot UMAP embedding colored by group
# (Useful for visualizing node positions separately)
group_colors = adata_filtered.obs['source_group'].map({'dim': 'blue', 'bright': 'red'})

plt.figure()
plt.scatter(embedding[:, 0], embedding[:, 1], c=group_colors, s=10, alpha=0.7)
plt.title('UMAP: dim vs bright (colored)')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.show()
