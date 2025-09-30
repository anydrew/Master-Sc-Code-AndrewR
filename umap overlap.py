# -*- coding: utf-8 -*-
"""
Created on Wed Jul  2 18:45:56 2025

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd
import numpy as np
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
# Subset relevant cells
# pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
# gbm_source = ['TiCD56dim_glioblastoma']
myfilter5 = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_filtered = adata[myfilter5].copy()

# Step 2: Label source
adata_filtered.obs['source_group'] = adata_filtered.obs['source'].apply(
    lambda x: 'PBMC' if x == 'PBMC' else 'glioblastoma'
)

# Reuse UMAP from original (if available)
adata_filtered.obsm['X_umap'] = adata[myfilter5].obsm['X_umap']

# Round UMAP to find overlaps
umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Group by rounded coordinates and source_group
group_counts = adata_filtered.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)

# Define overlapping coordinate keys (both PBMC and Glioblastoma present)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['glioblastoma'] > 0)].index

# Assign color labels for visualization
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# Define color map
color_dict = {'PBMC': 'blue', 'glioblastoma': 'red', 'Overlapping': 'green'}

# Plot UMAP with overlaps highlighted
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')

# --- Calculate Jaccard index ---

# Get sets of unique rounded UMAP coordinates for each group
pbmc_coords = set(group_counts[group_counts['PBMC'] > 0].index)
glioblastoma_coords = set(group_counts[group_counts['glioblastoma'] > 0].index)

# Calculate intersection and union
intersection = pbmc_coords.intersection(glioblastoma_coords)
union = pbmc_coords.union(glioblastoma_coords)

# Calculate Jaccard index
jaccard_index = len(intersection) / len(union) if len(union) > 0 else np.nan

print(f"Jaccard Index for PBMC vs Glioblastoma UMAP overlaps: {jaccard_index:.4f}")
