# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 14:37:56 2025

@author: AROCHE6
"""
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd, numpy as np, matplotlib.pyplot as plt, scipy.stats as ss, math
import scipy.spatial.distance as distance
from sklearn.preprocessing import normalize
from umap import UMAP
from sklearn.manifold import TSNE, SpectralEmbedding
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
import scanpy as sc
import pandas as pd

# Step 1: Define sources
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']

# Step 2: Filter data
myfilter = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter].copy()

# Step 3: Label groups for coloring
adata_filtered.obs['color_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Step 4: Use existing UMAP if available
adata_filtered.obsm['X_umap'] = adata[myfilter].obsm['X_umap']

# Step 5: Compute PCA and neighbors for connectivity
sc.pp.pca(adata_filtered, n_comps=20)
sc.pp.neighbors(adata_filtered, use_rep='X_pca')

# Step 6: Plot UMAP with connectivity edges
color_dict = {'PBMC': 'red', 'Glioblastoma': 'blue'}
sc.pl.umap(
    adata_filtered,
    color='color_group',
    palette=color_dict,
    edges=True,
    edges_width=0.5,
    title='UMAP: PBMC Dim vs Glioblastoma Dim with Connectivity'
)
