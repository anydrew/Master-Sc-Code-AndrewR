# -*- coding: utf-8 -*-
"""
UMAP of TiCD56bright vs TiCD56dim Glioblastoma NK cells with connectivity
(no umap.plot)
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Load AnnData object
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)

# Step 1: Filter for glioblastoma bright and dim NK subsets
myfilter = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
    (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
)
adata_filtered = adata[myfilter].copy()

# Step 2: Label bright vs dim
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'bright' if 'bright' in x else 'dim'
)

# Step 3: PCA + neighbors + UMAP
sc.pp.pca(adata_filtered, n_comps=20)
sc.pp.neighbors(adata_filtered, use_rep='X_pca')
sc.tl.umap(adata_filtered)

# Step 4: Plot UMAP with edges
sc.pl.umap(
    adata_filtered,
    color='source_group',
    palette={'dim': 'blue', 'bright': 'red'},
    edges=True,
    edges_width=0.5,
    title='UMAP: TiCD56bright vs TiCD56dim (glioblastoma) with connectivity'
)

