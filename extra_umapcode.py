# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 19:34:14 2025

@author: AROCHE6
"""

file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
filtered_adata = adata[myfilter6]
sc.pl.umap(filtered_adata, color='subset_source')
filtered_adata.uns['subset_source_colors'] = ['red', 'blue']  # Make sure order aligns with subset categories

# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')
runfile('C:/Users/aroche6/abc/umap_file.py', wdir='C:/Users/aroche6/abc')
myfilter7 = ((adata.obs['subset_source'] == 'Adaptive_PBMC') | (adata.obs['subset_source'] == 'CD57+_PBMC') | (adata.obs['subset_source'] == 'KIR+_PBMC') | (adata.obs['subset_source'] == 'NKG2A+_PBMC') | (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma'))
filtered_adata = adata[myfilter7]
filtered_adata.uns['subset_source_colors'] = ['red', 'red', 'red', 'blue']  # Make sure order aligns with subset categories
sc.pl.umap(filtered_adata, color='subset_source')
sc.pl.umap(filtered_adata)
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
filtered_adata = adata[myfilter6]
filtered_adata.uns['subset_source_colors'] = ['red', 'blue']  # Make sure order aligns with subset categories

# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')
runfile('C:/Users/aroche6/abc/umap_file.py', wdir='C:/Users/aroche6/abc')
filtered_adata.uns['bright_dim_subset_colors'] = ['red', 'blue']
sc.pl.umap(filtered_adata, color='bright_dim_subset')
filtered_adata.uns['bright_dim_subset_colors', 'TiCD56dim_glioblastoma_colors'] = ['red', 'blue']

sc.pl.umap(filtered_adata)
filtered_adata.uns['bright_dim_subset_colors', 'TiCD56dim_glioblastoma_colors'] = ['red', 'blue']

sc.pl.umap(filtered_adata, color='bright_dim_subset, TiCD56dim_glioblastoma_colors')
sc.pl.umap(filtered_adata, color= ['red', 'blue'])
myfilter7 = ((adata.obs['subset_source'] == 'Adaptive_PBMC') | (adata.obs['subset_source'] == 'CD57+_PBMC') | (adata.obs['subset_source'] == 'KIR+_PBMC') | (adata.obs['subset_source'] == 'NKG2A+_PBMC') | (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma'))
filtered_adata = adata[myfilter7]
import scanpy as sc
import matplotlib.pyplot as plt

# Define the PBMC and glioblastoma sources
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']

# Create a filtered AnnData object with relevant cells
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Create a new column to define the group (PBMC or GBM)
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)
sc.pp.pca(adata_filtered)
sc.pp.neighbors(adata_filtered)
sc.tl.umap(adata_filtered)

# Plot UMAP colored by group
sc.pl.umap(adata_filtered, color='source_group', palette=['blue', 'red'], title='UMAP: PBMC vs Glioblastoma dim cells')
import scanpy as sc
import matplotlib.pyplot as plt

# Define the PBMC and glioblastoma sources
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']

# Create a filtered AnnData object with relevant cells
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Create a new column to define the group (PBMC or GBM)
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)

# Plot UMAP colored by group
sc.pl.umap(adata_filtered, color='source_group', palette=['blue', 'red'], title='UMAP: PBMC vs Glioblastoma dim cells')
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
filtered_adata = adata[myfilter6]
sc.pl.umap(filtered_adata, color='subset_source')
# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)
sc.pp.pca(adata_filtered)
sc.pp.neighbors(adata_filtered)
sc.tl.umap(adata_filtered)  # Make sure order aligns with subset categories
filtered_adata.uns['subset_source_colors'] = ['red', 'blue']
# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')
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
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
filtered_adata = adata[myfilter6]
sc.pl.umap(filtered_adata, color='subset_source')
# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)
sc.pp.pca(adata_filtered)
sc.pp.neighbors(adata_filtered)
sc.tl.umap(adata_filtered)  # Make sure order aligns with subset categories
filtered_adata.uns['subset_source_colors'] = ['red', 'blue']
# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
filtered_adata = adata[myfilter6]
#sc.pl.umap(filtered_adata, color='subset_source')
# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)
sc.pp.pca(filtered_adata)
sc.pp.neighbors(filtered_adata)
sc.tl.umap(filtered_adata)  # Make sure order aligns with subset categories
filtered_adata.uns['subset_source_colors'] = ['red', 'blue']
# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')
myfilter6 = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
    (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
)
filtered_adata = adata[myfilter6]
#sc.pl.umap(filtered_adata, color='subset_source')
# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)
sc.pp.pca(filtered_adata)
sc.pp.neighbors(filtered_adata)
sc.tl.umap(filtered_adata)  # Make sure order aligns with subset categories
filtered_adata.uns['subset_source_colors'] = ['red', 'blue']
# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')

## ---(Thu Jun 26 13:04:09 2025)---
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd, numpy as np, matplotlib.pyplot as plt, scipy.stats as ss, math
import scipy.spatial.distance as distance
from sklearn.preprocessing import normalize
from umap import UMAP
from sklearn.manifold import TSNE, SpectralEmbedding

file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
myfilter6 = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
)
myfilter6 = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
adata = sc.read_h5ad(file_path)
myfilter6 = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
filtered_adata = adata[myfilter6]
filtered_adata.uns['subset_source_colors'] = ['blue']
sc.pl.umap(filtered_adata, color='subset_source')
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']

# Create a filtered AnnData object with relevant cells
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Create a new column to define the group (PBMC or GBM)
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)

# Plot UMAP colored by group
sc.pl.umap(adata_filtered, color='source_group', palette=['blue', 'red'], title='UMAP: PBMC vs Glioblastoma dim cells')
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']

# Create a filtered AnnData object with relevant cells
myfilter8 = adata.obs['subset_source'].isin(pbmc_sources)
adata_filtered = adata[myfilter8].copy()

# Create a new column to define the group (PBMC or GBM)
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources
)

# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)

# Plot UMAP colored by group
sc.pl.umap(adata_filtered, color='source_group', palette=['red'], title='UMAP: PBMC vs Glioblastoma dim cells')
myfilter8 = adata.obs['subset_source'].isin(pbmc_sources)
adata_filtered = adata[myfilter8].copy()

# Create a new column to define the group (PBMC or GBM)
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)

# Plot UMAP colored by group
sc.pl.umap(adata_filtered, color='source_group', palette=['red'],title='UMAP: PBMC vs Glioblastoma dim cells')
myfilter9 = adata.obs['subset_source'].isin(gbm_source)
adata_filtered = adata[myfilter9].copy()

# Create a new column to define the group (PBMC or GBM)
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)

# Plot UMAP colored by group
sc.pl.umap(adata_filtered, color='source_group', palette=['blue'],title='UMAP: PBMC vs Glioblastoma dim cells')
import scanpy as sc
import pandas as pd
import numpy as np

# Subset relevant cells
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Label PBMC vs Glioblastoma
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Reuse UMAP from original (if available)
adata_filtered.obsm['X_umap'] = adata[myfilter7].obsm['X_umap']

# Round UMAP to find overlaps
umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Group by rounded coordinates and source_group
group_counts = adata_filtered.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)

# Define overlapping coordinate keys (both PBMC and Glioblastoma present)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['Glioblastoma'] > 0)].index

# Assign color labels
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# Define color map
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'green'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'yellow'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'black'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'gold'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'red', 'Glioblastoma': 'blue', 'Overlapping': 'gold'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
import scanpy as sc
import pandas as pd
import numpy as np
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
# Subset relevant cells

adata_filtered = adata[myfilter6].copy()


# Reuse UMAP from original (if available)
adata_filtered.obsm['X_umap'] = adata[myfilter7].obsm['X_umap']

# Round UMAP to find overlaps
umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Group by rounded coordinates and source_group
group_counts = adata_filtered.obs.groupby(['umap_rounded', 'subset_source']).size().unstack(fill_value=0)

# Define overlapping coordinate keys (both PBMC and Glioblastoma present)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['Glioblastoma'] > 0)].index

# Assign color labels
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# Define color map
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'green'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
import scanpy as sc
import pandas as pd
import numpy as np

# Subset relevant cells
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source = ['TiCD56dim_glioblastoma']
myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_filtered = adata[myfilter7].copy()

# Label PBMC vs Glioblastoma
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
)

# Reuse UMAP from original (if available)
adata_filtered.obsm['X_umap'] = adata[myfilter7].obsm['X_umap']

# Round UMAP to find overlaps
umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Group by rounded coordinates and source_group
group_counts = adata_filtered.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)

# Define overlapping coordinate keys (both PBMC and Glioblastoma present)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['Glioblastoma'] > 0)].index

# Assign color labels
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# Define color map
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'green'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'yellow'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'black'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'gold'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'PBMC': 'red', 'Glioblastoma': 'blue', 'Overlapping': 'gold'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
sc.pl.umap(adata_filtered, color='Overlapping', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
# Subset only overlapping cells
overlapping_adata = adata_filtered[adata_filtered.obs['color_group'] == 'Overlapping'].copy()

# Plot UMAP with only overlapping cells in green
sc.pl.umap(
    overlapping_adata,
    color='color_group',
    palette={'Overlapping': 'green'},
    title='UMAP: Overlapping PBMC & Glioblastoma dim cells (Green only)'
)
# Subset only overlapping cells
overlapping_adata = adata_filtered[adata_filtered.obs['color_group'] == 'Overlapping'].copy()

# Plot UMAP with only overlapping cells in green
sc.pl.umap(
    overlapping_adata,
    color='color_group',
    palette={'Overlapping': 'yellow'},
    title='UMAP: Overlapping PBMC & Glioblastoma dim cells (Green only)'
)
sc.pl.umap(
    overlapping_adata,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='UMAP: Overlapping PBMC & Glioblastoma dim cells (Green only)'
)
import scanpy as sc
import pandas as pd

# Step 1: Filter data
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
adata_bright = adata[myfilter6].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['subset_source'].apply(
    lambda x: 'PBMC' if x == 'CD56bright_PBMC' else 'Glioblastoma'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['Glioblastoma'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'PBMC': 'red', 'Glioblastoma': 'blue', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
import scanpy as sc
import pandas as pd

# Step 1: Filter data
myfilter6 = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
    (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
)
adata_bright = adata[myfilter6].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['subset_source'].apply(
    lambda x: 'Bright' if x == 'TiCD56bright_glioblastoma' else 'Dim'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['Bright'] > 0) & (group_counts['Dim'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'Bright': 'red', 'Dim': 'blue', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
import scanpy as sc
import pandas as pd

# Step 1: Filter data
myfilter5 = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_bright = adata[myfilter5].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['source'].apply(
    lambda x: 'PBMC' if x == 'HC-PBMC-NK' else 'GBM-NK'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['HC-PBMC-NK'] > 0) & (group_counts['GBM-NK'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'HC-PBMC-NK': 'blue', 'GBM-NK': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)

## ---(Thu Jun 26 17:19:43 2025)---
import scanpy as sc
import pandas as pd
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)# Step 1: Filter data
myfilter5 = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_bright = adata[myfilter5].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['source'].apply(
    lambda x: 'PBMC' if x == 'HC-PBMC-NK' else 'GBM-NK'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['HC-PBMC-NK'] > 0) & (group_counts['GBM-NK'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'HC-PBMC-NK': 'blue', 'GBM-NK': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
import scanpy as sc
import pandas as pd
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)# Step 1: Filter data
myfilter5 = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_bright = adata[myfilter5].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['source'].apply(
    lambda x: 'PBMC' if x == 'HC-PBMC-NK' else 'GBM-NK'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter5].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['HC-PBMC-NK'] > 0) & (group_counts['GBM-NK'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'HC-PBMC-NK': 'blue', 'GBM-NK': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
import scanpy as sc
import pandas as pd
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)# Step 1: Filter data
myfilter5 = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_bright = adata[myfilter5].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['source'].apply(
    lambda x: 'PBMC' if x == 'PBMC' else 'glioblastoma'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter5].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['glioblastoma'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
import scanpy as sc
import pandas as pd
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)# Step 1: Filter data
import scanpy as sc
import pandas as pd
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)# Step 1: Filter data
myfilter5 = (
    (adata.obs['source'] == 'PBMC') |
    (adata.obs['source'] == 'glioblastoma')
)
adata_bright = adata[myfilter5].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['source'].apply(
    lambda x: 'PBMC' if x == 'PBMC' else 'glioblastoma'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter5].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['glioblastoma'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
import scanpy as sc
import pandas as pd

# Step 1: Filter cells by source
myfilter_source = (adata.obs['source'] == 'PBMC') | (adata.obs['source'] == 'glioblastoma')
adata_source = adata[myfilter_source].copy()

# Step 2: Copy UMAP coordinates from full object (assuming precomputed)
adata_source.obsm['X_umap'] = adata[myfilter_source].obsm['X_umap']

# Step 3: Round UMAP coordinates to identify visual overlaps
umap_coords = pd.DataFrame(adata_source.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_source.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 4: Identify overlaps based on coordinate grouping
group_counts = adata_source.obs.groupby(['umap_rounded', 'source']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['glioblastoma'] > 0)].index

# Step 5: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source']

adata_source.obs['color_group'] = adata_source.obs.apply(assign_color, axis=1)

# Step 6: Plot 1 – All cells
color_dict = {'PBMC': 'blue', 'glioblastoma': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_source,
    color='color_group',
    palette=color_dict,
    title='UMAP: All PBMC vs Glioblastoma Cells (Overlap = Gold)'
)

# Step 7: Plot 2 – Only overlapping cells
overlapping_only = adata_source[adata_source.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='UMAP: Overlapping Cells Only (Gold)'
)

color_dict = {'PBMC': 'blue', 'glioblastoma': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    'PBMC',
    color='color_group',
    palette=color_dict,
    title='UMAP: All PBMC vs Glioblastoma Cells (Overlap = Gold)'
)
# Step 1: Filter cells by source
myfilter_source = (adata.obs['source'] == 'PBMC')
adata_source = adata[myfilter_source].copy()

# Step 2: Copy UMAP coordinates from full object (assuming precomputed)
#adata_source.obsm['X_umap'] = adata[myfilter_source].obsm['X_umap']

# Step 3: Round UMAP coordinates to identify visual overlaps
#umap_coords = pd.DataFrame(adata_source.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
#umap_coords_rounded = umap_coords.round(2)
#adata_source.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 4: Identify overlaps based on coordinate grouping
#group_counts = adata_source.obs.groupby(['umap_rounded', 'source']).size().unstack(fill_value=0)
#overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['glioblastoma'] > 0)].index

# Step 5: Assign color group
#def assign_color(row):
 #   if row['umap_rounded'] in overlapping_keys:
  #      return 'Overlapping'
   # return row['source']

adata_source.obs['color_group'] = adata_source.obs.apply(assign_color, axis=1)

# Step 6: Plot 1 – All cells
color_dict = {'PBMC': 'blue', 'glioblastoma': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_source,
    color='color_group',
    palette=color_dict,
    title='UMAP: All PBMC vs Glioblastoma Cells (Overlap = Gold)'
)
# Step 1: Filter cells by source
myfilter_source = (adata.obs['source'] == 'PBMC')
adata_source = adata[myfilter_source].copy()


adata_source.obs['color_group'] = adata_source.obs.apply(assign_color, axis=1)

# Step 6: Plot 1 – All cells
color_dict = {'PBMC': 'blue'}
sc.pl.umap(
    adata_source,
    color='color_group',
    palette=color_dict,
    title='UMAP: All PBMC vs Glioblastoma Cells (Overlap = Gold)'
)
myfilter70 = ((adata.obs['source'] == 'PBMC'))
filtered_adata = adata[myfilter70]
filtered_adata.uns['subset_source_colors'] = ['blue'] 
sc.pl.umap(filtered_adata, color='subset_source')
sc.pl.umap(filtered_adata, color='source')
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

myfilter70 = ((adata.obs['source'] == 'PBMC'))
filtered_adata = adata[myfilter70]
filtered_adata.uns['source_colors'] = ['blue'] 
sc.pl.umap(filtered_adata, color='source')
myfilter70 = ((adata.obs['source'] == 'glioblastoma'))
filtered_adata = adata[myfilter70]
filtered_adata.uns['source_colors'] = ['red'] 
sc.pl.umap(filtered_adata, color='source')
dim_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
bright_source = ['CD56bright_PBMC']
myfilter7 = adata.obs['subset_source'].isin(dim_sources + bright_source)
adata_filtered = adata[myfilter7].copy()

# Label PBMC vs Glioblastoma
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'dim' if x in dim_sources else 'bright'
)

# Reuse UMAP from original (if available)
adata_filtered.obsm['X_umap'] = adata[myfilter7].obsm['X_umap']

# Round UMAP to find overlaps
umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Group by rounded coordinates and source_group
group_counts = adata_filtered.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)

# Define overlapping coordinate keys (both PBMC and Glioblastoma present)
overlapping_keys = group_counts[(group_counts['dim'] > 0) & (group_counts['bright'] > 0)].index

# Assign color labels
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# Define color map
color_dict = {'dim': 'blue', 'bright': 'red', 'Overlapping': 'green'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'dim': 'blue', 'bright': 'red', 'Overlapping': 'gold'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
overlapping_only = adata_source[adata_source.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='UMAP: Overlapping Cells Only (Gold)'
)
# Subset to overlapping cells only
overlapping_only = adata_filtered[adata_filtered.obs['color_group'] == 'Overlapping'].copy()

# Plot UMAP with overlapping cells (green)
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'green'},
    title='UMAP: Overlapping dim and bright PBMC Cells (Green Only)'
)
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='UMAP: Overlapping dim and bright PBMC Cells (Green Only)'
)
color_dict = {'dim': 'red', 'bright': 'blue', 'Overlapping': 'gold'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
myfilter6 = (
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
    (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
)
adata_bright = adata[myfilter6].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['subset_source'].apply(
    lambda x: 'Bright' if x == 'TiCD56bright_glioblastoma' else 'Dim'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['Bright'] > 0) & (group_counts['Dim'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'Bright': 'red', 'Dim': 'blue', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)

# Step 8: Plot 2 – Just overlapping cells
overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
sc.pl.umap(
    overlapping_only,
    color='color_group',
    palette={'Overlapping': 'gold'},
    title='CD56bright Overlapping Cells Only (Gold)'
)
color_dict = {'Bright': 'blue', 'Dim': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)
# Subset to just Dim cells
dim_only = adata_bright[adata_bright.obs['color_group'] == 'Dim'].copy()

# Plot UMAP with only Dim cells in red
sc.pl.umap(
    dim_only,
    color='color_group',
    palette={'Dim': 'red'},
    title='UMAP: Dim Cells Only (Red)'
)
sc.pl.umap(
    dim_only,
    color='color_group',
    palette={'Dim': 'blue'},
    title='UMAP: Dim Cells Only (Red)'
)
dim_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
bright_source = ['CD56bright_PBMC']
myfilter7 = adata.obs['subset_source'].isin(dim_sources + bright_source)
adata_filtered = adata[myfilter7].copy()

# Label PBMC vs Glioblastoma
adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
    lambda x: 'dim' if x in dim_sources else 'bright'
)

# Reuse UMAP from original (if available)
adata_filtered.obsm['X_umap'] = adata[myfilter7].obsm['X_umap']

# Round UMAP to find overlaps
umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Group by rounded coordinates and source_group
group_counts = adata_filtered.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)

# Define overlapping coordinate keys (both PBMC and Glioblastoma present)
overlapping_keys = group_counts[(group_counts['dim'] > 0) & (group_counts['bright'] > 0)].index

# Assign color labels
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# Define color map
color_dict = {'dim': 'blue', 'bright': 'red', 'Overlapping': 'green'}

# Plot
sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')
color_dict = {'dim': 'blue', 'bright': 'red', 'Overlapping': 'gold'}
dim_only = adata_bright[adata_bright.obs['color_group'] == 'dim'].copy()

# Plot UMAP with only Dim cells in red
sc.pl.umap(
    dim_only,
    color='color_group',
    palette={'Dim': 'blue'},
    title='UMAP: Dim Cells Only (Red)'
)
# Subset to only 'dim' cells
dim_only = adata_filtered[adata_filtered.obs['color_group'] == 'dim'].copy()

# Plot UMAP for dim cells only
sc.pl.umap(
    dim_only,
    color='color_group',
    palette={'dim': 'blue'},
    title='UMAP: dim Cells Only (Blue)'
)
bright_only = adata_filtered[adata_filtered.obs['color_group'] == 'bright'].copy()

# Plot UMAP for dim cells only
sc.pl.umap(
    bright_only,
    color='color_group',
    palette={'bright': 'blue'},
    title='UMAP: dim Cells Only (Blue)'
)
myfilter6 = (
    (adata.obs['subset_source'] == 'CD56bright_PBMC') |
    (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
)
adata_bright = adata[myfilter6].copy()

# Step 2: Label source
adata_bright.obs['source_group'] = adata_bright.obs['subset_source'].apply(
    lambda x: 'PBMC' if x == 'CD56bright_PBMC' else 'Glioblastoma'
)

# Step 3: Reuse UMAP from full data if available
adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# Step 5: Find overlapping coordinate keys
group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['Glioblastoma'] > 0)].index

# Step 6: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# Step 7: Plot 1 – All cells
color_dict = {'PBMC': 'blue', 'Glioblastoma': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_bright,
    color='color_group',
    palette=color_dict,
    title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
)