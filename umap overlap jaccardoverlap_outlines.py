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
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import pandas as pd

# import matplotlib.pyplot as plt
# import pandas as pd

# import matplotlib.pyplot as plt
# import pandas as pd
# import matplotlib.pyplot as plt
# import pandas as pd

# pbmc_sources = ['KIR+_PBMC','Adaptive_PBMC','CD57+_PBMC','NKG2A+_PBMC']
# gbm_source = ['TiCD56dim_glioblastoma']
# myfilter_pbmc_gbm = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
# adata_pbmc_gbm = adata[myfilter_pbmc_gbm].copy()

# X = adata_pbmc_gbm.obsm['X_umap']
# umap_df = pd.DataFrame(X, columns=['UMAP1','UMAP2']).round(2)
# adata_pbmc_gbm.obs['umap_rounded'] = list(umap_df.itertuples(index=False, name=None))

# adata_pbmc_gbm.obs['source_group'] = adata_pbmc_gbm.obs['subset_source'].apply(
#     lambda x: 'PBMC' if x in pbmc_sources else 'GBM'
# )

# group_counts = adata_pbmc_gbm.obs.groupby(['umap_rounded','source_group']).size().unstack(fill_value=0)
# overlap_keys = group_counts[(group_counts['PBMC']>0) & (group_counts['GBM']>0)].index

# adata_pbmc_gbm.obs['color_group'] = adata_pbmc_gbm.obs['umap_rounded'].apply(
#     lambda r: 'Overlapping' if r in overlap_keys else 'Background'
# )

# mask_ov = (adata_pbmc_gbm.obs['color_group']=='Overlapping')
# X_ov = X[mask_ov.values]
# X_bg = X[~mask_ov.values]

# plt.figure(figsize=(7,6))
# plt.scatter(X_bg[:,0], X_bg[:,1], c='lightgrey', s=20, alpha=0.4, edgecolors='none')
# plt.scatter(X_ov[:,0], X_ov[:,1], c='gold', s=70, alpha=0.9, edgecolors='black', linewidths=0.6)
# plt.title("PBMC NK subsets vs GBM CD56dim — Overlap highlighted")
# plt.gca().set_aspect('auto')
# plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
# plt.show()

# myfilter7 = (
#     (adata.obs['subset_source']=='CD56bright_PBMC') |
#     (adata.obs['subset_source']=='TiCD56bright_glioblastoma')
# )
# adata_bright_pbmc = adata[myfilter7].copy()

# X = adata_bright_pbmc.obsm['X_umap']
# umap_df = pd.DataFrame(X, columns=['UMAP1','UMAP2']).round(2)
# adata_bright_pbmc.obs['umap_rounded'] = list(umap_df.itertuples(index=False, name=None))

# group_counts = adata_bright_pbmc.obs.groupby(['umap_rounded','subset_source']).size().unstack(fill_value=0)
# overlap_keys = group_counts[
#     (group_counts['CD56bright_PBMC']>0) & (group_counts['TiCD56bright_glioblastoma']>0)
# ].index

# adata_bright_pbmc.obs['color_group'] = adata_bright_pbmc.obs['umap_rounded'].apply(
#     lambda r: 'Overlapping' if r in overlap_keys else 'Background'
# )

# mask_ov = (adata_bright_pbmc.obs['color_group']=='Overlapping')
# X_ov = X[mask_ov.values]
# X_bg = X[~mask_ov.values]

# plt.figure(figsize=(7,6))
# plt.scatter(X_bg[:,0], X_bg[:,1], c='lightgrey', s=20, alpha=0.4, edgecolors='none')
# plt.scatter(X_ov[:,0], X_ov[:,1], c='gold', s=70, alpha=0.9, edgecolors='black', linewidths=0.6)
# plt.title("CD56bright PBMC vs GBM CD56bright — Overlap highlighted")
# plt.gca().set_aspect('auto')
# plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
# plt.show()

# # Filter
# myfilter5 = (
#     (adata.obs['source'] == 'PBMC') |
#     (adata.obs['source'] == 'glioblastoma')
# )
# adata_global = adata[myfilter5].copy()

# # UMAP coords + rounded bins
# X = adata_global.obsm['X_umap']
# umap_df = pd.DataFrame(X, columns=['UMAP1','UMAP2']).round(2)
# adata_global.obs['umap_rounded'] = list(umap_df.itertuples(index=False, name=None))

# # Overlap keys (bins containing both groups)
# group_counts = adata_global.obs.groupby(['umap_rounded','source']).size().unstack(fill_value=0)
# overlap_keys = group_counts[(group_counts['PBMC']>0) & (group_counts['glioblastoma']>0)].index

# # Background vs Overlapping
# adata_global.obs['color_group'] = adata_global.obs['umap_rounded'].apply(
#     lambda r: 'Overlapping' if r in overlap_keys else 'Background'
# )

# # Split coords
# mask_ov = (adata_global.obs['color_group'] == 'Overlapping')
# X_ov = X[mask_ov.values]
# X_bg = X[~mask_ov.values]

# # Plot
# plt.figure(figsize=(7,6))
# plt.scatter(X_bg[:,0], X_bg[:,1], c='lightgrey', s=20, alpha=0.4, edgecolors='none')
# plt.scatter(X_ov[:,0], X_ov[:,1], c='gold', s=70, alpha=0.9, edgecolors='black', linewidths=0.6)
# plt.title("PBMC vs Glioblastoma — Overlap highlighted")
# plt.gca().set_aspect('auto')
# plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
# plt.show()


# dim_sources = ['KIR+_PBMC','Adaptive_PBMC','CD57+_PBMC','NKG2A+_PBMC']
# bright_source = ['CD56bright_PBMC']
# myfilter_dimbright = adata.obs['subset_source'].isin(dim_sources + bright_source)
# adata_dimbright = adata[myfilter_dimbright].copy()

# X = adata_dimbright.obsm['X_umap']
# umap_df = pd.DataFrame(X, columns=['UMAP1','UMAP2']).round(2)
# adata_dimbright.obs['umap_rounded'] = list(umap_df.itertuples(index=False, name=None))

# adata_dimbright.obs['source_group'] = adata_dimbright.obs['subset_source'].apply(
#     lambda x: 'Dim' if x in dim_sources else 'Bright'
# )

# group_counts = adata_dimbright.obs.groupby(['umap_rounded','source_group']).size().unstack(fill_value=0)
# overlap_keys = group_counts[(group_counts['Dim']>0) & (group_counts['Bright']>0)].index

# adata_dimbright.obs['color_group'] = adata_dimbright.obs['umap_rounded'].apply(
#     lambda r: 'Overlapping' if r in overlap_keys else 'Background'
# )

# mask_ov = (adata_dimbright.obs['color_group']=='Overlapping')
# X_ov = X[mask_ov.values]
# X_bg = X[~mask_ov.values]

# plt.figure(figsize=(7,6))
# plt.scatter(X_bg[:,0], X_bg[:,1], c='lightgrey', s=20, alpha=0.4, edgecolors='none')
# plt.scatter(X_ov[:,0], X_ov[:,1], c='gold', s=70, alpha=0.9, edgecolors='black', linewidths=0.6)
# plt.title("PBMC CD56dim-like vs CD56bright — Overlap highlighted")
# plt.gca().set_aspect('auto')
# plt.xlabel("UMAP1"); plt.ylabel("UMAP2")
# plt.show()


# dim_sources = ['KIR+_PBMC','Adaptive_PBMC','CD57+_PBMC','NKG2A+_PBMC']
# bright_source = ['CD56bright_PBMC']

# myfilter_dimbright = adata.obs['subset_source'].isin(dim_sources+bright_source)
# adata_dimbright = adata[myfilter_dimbright].copy()

# # myfilter6 = ( 
# #     (adata.obs['subset_source']=='TiCD56bright_glioblastoma') |
# #     (adata.obs['subset_source']=='TiCD56dim_glioblastoma')
# # )
# # adata_gbm_brightdim = adata[myfilter6].copy()

# # UMAP coords
# X = adata_dimbright.obsm['X_umap']
# umap_coords = pd.DataFrame(X, columns=['UMAP1','UMAP2']).round(2)
# adata_dimbright.obs['umap_rounded'] = list(umap_coords.itertuples(index=False, name=None))

# # Overlap detection
# group_counts = adata_dimbright.obs.groupby(
#     ['umap_rounded','subset_source']
# ).size().unstack(fill_value=0)
# overlap_keys = group_counts[
#     (group_counts['dim_sources']>0)&
#     (group_counts['bright_source']>0)
# ].index

# adata_dimbright.obs['color_group'] = adata_dimbright.obs['umap_rounded'].apply(
#     lambda r: 'Overlapping' if r in overlap_keys else 'Background'
# )

# # Split coords
# overlap_mask = adata_dimbright.obs['color_group']=='Overlapping'
# X_overlap = X[overlap_mask.values]
# X_background = X[~overlap_mask.values]

# # Plot background (light grey, less transparent)
# plt.scatter(X_background[:,0], X_background[:,1],
#             c='lightgrey', s=20, alpha=0.4,
#             edgecolors='none')

# # Plot overlaps (gold with black outline)
# plt.scatter(X_overlap[:,0], X_overlap[:,1],
#             c='gold', s=70, alpha=0.9,
#             edgecolors='black', linewidths=0.6)

# plt.title("Glioblastoma CD56bright vs CD56dim Overlap")
# plt.gca().set_aspect('auto')
# plt.show()


# # Step 1: Filter cells
# myfilter6 = (
#     (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
#     (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
# )
# adata_gbm_brightdim = adata[myfilter6].copy()
# adata_gbm_brightdim.obs['source_group'] = adata_gbm_brightdim.obs['subset_source']

# # Step 2: Reuse UMAP from the same cells (keeps size/shape consistent)
# adata_gbm_brightdim.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# # Step 3: Detect overlaps (rounded UMAP bins)
# umap_df = pd.DataFrame(adata_gbm_brightdim.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
# umap_rounded = umap_df.round(2)
# adata_gbm_brightdim.obs['umap_rounded'] = list(umap_rounded.itertuples(index=False, name=None))

# group_counts = (
#     adata_gbm_brightdim.obs
#     .groupby(['umap_rounded', 'source_group'])
#     .size()
#     .unstack(fill_value=0)
# )

# overlap_keys = group_counts[
#     (group_counts['TiCD56bright_glioblastoma'] > 0) &
#     (group_counts['TiCD56dim_glioblastoma'] > 0)
# ].index

# # Step 4: Background vs Overlapping labels
# adata_gbm_brightdim.obs['color_group'] = adata_gbm_brightdim.obs.apply(
#     lambda r: 'Overlapping' if r['umap_rounded'] in overlap_keys else 'Background',
#     axis=1
# )

# # Step 5: Plot with outlines ONLY for overlaps
# coords = adata_gbm_brightdim.obsm['X_umap']
# groups = adata_gbm_brightdim.obs['color_group'].values

# fig, ax = plt.subplots(figsize=(7, 6))

# # Background (no outline)
# mask_bg = (groups == 'Background')
# ax.scatter(
# #     coords[mask_bg, 0], coords[mask_bg, 1],
# #     c='lightgrey', s=20, alpha=0.35, edgecolors='none', label='Background'
# # )

# # # Overlaps (gold with dark outline)
# # mask_ov = (groups == 'Overlapping')
# # ax.scatter(
# #     coords[mask_ov, 0], coords[mask_ov, 1],
# #     c='gold', s=60, alpha=0.9,
# #     edgecolors='black', linewidths=0.6,
# #     label='Overlapping'
# # )

# # ax.set_xlabel('UMAP1')
# # ax.set_ylabel('UMAP2')
# # ax.set_aspect('equal', adjustable='box')
# # ax.set_title('Glioblastoma CD56bright vs CD56dim — Overlap highlighted')
# # ax.legend(frameon=False)
# # plt.show()

# # import scanpy as sc
# # import pandas as pd
# # import matplotlib.pyplot as plt
# # import pandas as pd

# # import scanpy as sc
# # import pandas as pd
# # import matplotlib.pyplot as plt

# # import scanpy as sc
# # import pandas as pd

# # Step 1: Define groups
# dim_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
# bright_source = ['CD56bright_PBMC']

# # Step 2: Filter cells from both groups
# myfilter = adata.obs['subset_source'].isin(dim_sources + bright_source)
# adata_dimbright = adata[myfilter].copy()

# # Step 3: Label source group
# adata_dimbright.obs['source_group'] = adata_dimbright.obs['subset_source'].apply(
#     lambda x: 'Dim' if x in dim_sources else 'Bright'
# )

# # Step 4: Reuse original UMAP coordinates (keeps same size/shape as full dataset)
# adata_dimbright.obsm['X_umap'] = adata[myfilter].obsm['X_umap']

# # Step 5: Detect overlaps (rounded UMAP bins)
# umap_coords = pd.DataFrame(adata_dimbright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
# umap_coords_rounded = umap_coords.round(2)
# adata_dimbright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# group_counts = adata_dimbright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
# overlapping_keys = group_counts[(group_counts['Dim'] > 0) & (group_counts['Bright'] > 0)].index

# # Step 6: Assign color group
# def assign_color(row):
#     if row['umap_rounded'] in overlapping_keys:
#         return 'Overlapping'
#     return row['source_group']

# adata_dimbright.obs['color_group'] = adata_dimbright.obs.apply(assign_color, axis=1)

# # Step 7: Plot all cells with overlaps highlighted
# color_dict = {'Dim': 'blue', 'Bright': 'red', 'Overlapping': 'white'}
# sc.pl.umap(
#     adata_dimbright,
#     color='color_group',
#     palette=color_dict,
#     title='PBMC Dim vs Bright (Overlap = Gold)',
#     size=10
# )

# # Step 8: Plot only overlapping cells but keep UMAP layout
# sc.pl.umap(
#     adata_dimbright,
#     color='color_group',
#     palette={'Dim': 'white', 'Bright': 'white', 'Overlapping': 'gold'},
#     title='PBMC Dim vs Bright Overlapping Cells (Gold)',
#     size=80   # make yellow dots stand out
# )


# Step 1: Filter cells
# myfilter6 = (
#     (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
#     (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
# )
# adata_gbm_brightdim = adata[myfilter6].copy()
# adata_gbm_brightdim.obs['source_group'] = adata_gbm_brightdim.obs['subset_source']

# # Step 2: Reuse UMAP from the same cells (keeps size/shape consistent)
# adata_gbm_brightdim.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# # Step 3: Detect overlaps (rounded UMAP bins)
# umap_df = pd.DataFrame(adata_gbm_brightdim.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
# umap_rounded = umap_df.round(2)
# adata_gbm_brightdim.obs['umap_rounded'] = list(umap_rounded.itertuples(index=False, name=None))

# group_counts = (
#     adata_gbm_brightdim.obs
#     .groupby(['umap_rounded', 'source_group'])
#     .size()
#     .unstack(fill_value=0)
# )

# overlap_keys = group_counts[
#     (group_counts['TiCD56bright_glioblastoma'] > 0) &
#     (group_counts['TiCD56dim_glioblastoma'] > 0)
# ].index

# # Step 4: Background vs Overlapping labels
# adata_gbm_brightdim.obs['color_group'] = adata_gbm_brightdim.obs.apply(
#     lambda r: 'Overlapping' if r['umap_rounded'] in overlap_keys else 'Background',
#     axis=1
# )

# # Step 5: Plot with outlines ONLY for overlaps
# coords = adata_gbm_brightdim.obsm['X_umap']
# groups = adata_gbm_brightdim.obs['color_group'].values

# fig, ax = plt.subplots(figsize=(7, 6))

# # Background (no outline)
# mask_bg = (groups == 'Background')
# ax.scatter(
#     coords[mask_bg, 0], coords[mask_bg, 1],
#     c='lightgrey', s=20, alpha=0.35, edgecolors='none', label='Background'
# )

# # Overlaps (gold with dark outline)
# mask_ov = (groups == 'Overlapping')
# ax.scatter(
#     coords[mask_ov, 0], coords[mask_ov, 1],
#     c='gold', s=60, alpha=0.9,
#     edgecolors='black', linewidths=0.6,
#     label='Overlapping'
# )

# ax.set_xlabel('UMAP1')
# ax.set_ylabel('UMAP2')
# ax.set_aspect('equal', adjustable='box')
# ax.set_title('Glioblastoma CD56bright vs CD56dim — Overlap highlighted')
# ax.legend(frameon=False)
# plt.show()



# myfilter6 = (
#     (adata.obs['subset_source']=='TiCD56bright_glioblastoma') |
#     (adata.obs['subset_source']=='TiCD56dim_glioblastoma')
# )
# adata_gbm_brightdim = adata[myfilter6].copy()
# adata_gbm_brightdim.obs['source_group'] = adata_gbm_brightdim.obs['subset_source']

# # Reuse UMAP
# adata_gbm_brightdim.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# # # Overlaps
# # umap_coords = pd.DataFrame(adata_gbm_brightdim.obsm['X_umap'],columns=['UMAP1','UMAP2']).round(2)
# # adata_gbm_brightdim.obs['umap_rounded'] = list(umap_coords.itertuples(index=False, name=None))
# # group_counts = adata_gbm_brightdim.obs.groupby(['umap_rounded','source_group']).size().unstack(fill_value=0)
# # overlap_keys = group_counts[(group_counts['TiCD56bright_glioblastoma']>0)&
# #                             (group_counts['TiCD56dim_glioblastoma']>0)].index

# # adata_gbm_brightdim.obs['color_group'] = adata_gbm_brightdim.obs.apply(
# #     lambda r: 'Overlapping' if r['umap_rounded'] in overlap_keys else 'Background', axis=1
# # )

# # sc.pl.umap(adata_gbm_brightdim, color='color_group',
# #             palette={'Background':'lightgrey','Overlapping':'gold'},
# #             size=40, alpha=0.4, title='Glioblastoma CD56bright vs CD56dim Overlap')

import scanpy as sc
import pandas as pd

# Define groups: PBMC dim-like vs GBM dim
pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
gbm_source  = ['TiCD56dim_glioblastoma']

# Step 1: Filter cells (dim vs dim)
myfilter6 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
adata_dimdim = adata[myfilter6].copy()

# Step 2: Label source group
adata_dimdim.obs['source_group'] = adata_dimdim.obs['subset_source'].apply(
    lambda x: 'PBMC_dim' if x in pbmc_sources else 'GBM_dim'
)

# Step 3: Reuse original UMAP coordinates (keeps shape/size consistent)
adata_dimdim.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# Step 4: Detect overlaps (rounded UMAP bins)
umap_coords = pd.DataFrame(adata_dimdim.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
umap_coords_rounded = umap_coords.round(2)
adata_dimdim.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

group_counts = adata_dimdim.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
overlapping_keys = group_counts[(group_counts['PBMC_dim'] > 0) & (group_counts['GBM_dim'] > 0)].index

# Step 5: Assign color group
def assign_color(row):
    if row['umap_rounded'] in overlapping_keys:
        return 'Overlapping'
    return row['source_group']

adata_dimdim.obs['color_group'] = adata_dimdim.obs.apply(assign_color, axis=1)

# Step 6: Plot all cells with overlaps highlighted
color_dict = {'PBMC_dim': 'blue', 'GBM_dim': 'red', 'Overlapping': 'gold'}
sc.pl.umap(
    adata_dimdim,
    color='color_group',
    palette=color_dict,
    title='Dim PBMC vs Dim GBM (Overlap = Gold)',
    size=20
)



# # Step 1: Filter cells
# myfilter6 = (
#     (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') |
#     (adata.obs['subset_source'] == 'CD56bright_PBMC')
# )
# adata_gbm = adata[myfilter6].copy()

# # Step 2: Label source group
# adata_gbm.obs['source_group'] = adata_gbm.obs['subset_source'].apply(
#     lambda x: 'Bright' if x == 'TiCD56bright_glioblastoma' else 'Dim'
# )

# # Step 3: Reuse original UMAP coordinates (keeps shape/size consistent)
# adata_gbm.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# # Step 4: Detect overlaps (rounded UMAP bins)
# umap_coords = pd.DataFrame(adata_gbm.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
# umap_coords_rounded = umap_coords.round(2)
# adata_gbm.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# group_counts = adata_gbm.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
# overlapping_keys = group_counts[(group_counts['Bright'] > 0) & (group_counts['Dim'] > 0)].index

# # Step 5: Assign color group
# def assign_color(row):
#     if row['umap_rounded'] in overlapping_keys:
#         return 'Overlapping'
#     return row['source_group']

# adata_gbm.obs['color_group'] = adata_gbm.obs.apply(assign_color, axis=1)

# # Step 6: Plot all cells with overlaps highlighted
# color_dict = {'Bright': 'red', 'Dim': 'blue', 'Overlapping': 'white'}
# sc.pl.umap(
#     adata_gbm,
#     color='color_group',
#     palette=color_dict,
#     title='Glioblastoma: TiCD56bright vs TiCD56dim (Overlap = Gold)',
#     size=20
# )

# # # Step 7: Plot only overlapping cells but keep UMAP layout
# # sc.pl.umap(
# #     adata_gbm,
# #     color='color_group',
# #     palette={'Bright': 'red', 'Dim': 'lightgrey', 'Overlapping': 'gold'},
# #     edgecolors='black', linewidths=6, 
# #     title='Glioblastoma TiCD56bright vs Dim Overlaps (Gold)',
# #     size=80   # highlight overlap with bigger yellow dots
# # )
# # # Step 7: Plot only overlapping cells but keep UMAP layout
# # # sc.pl.umap(
# # #     adata_gbm,
# # #     color='color_group',
# # #     palette={'Bright': 'white', 'Dim': 'white', 'Overlapping': 'gold'},
# # #     title='Glioblastoma TiCD56bright vs Dim Overlaps (Gold)',
# # #     size=20   # highlight overlap with bigger yellow dots
# # # )


# # # myfilter5 = (
# # #     (adata.obs['source'] == 'PBMC') |
# # #     (adata.obs['source'] == 'glioblastoma')
# # # )
# # # adata_filtered = adata[myfilter5].copy()

# # # # Step 2: Label source
# # # adata_filtered.obs['source_group'] = adata_filtered.obs['source'].apply(
# # #     lambda x: 'PBMC' if x == 'PBMC' else 'glioblastoma'
# # # )

# # # # Reuse UMAP from original (if available)
# # # adata_filtered.obsm['X_umap'] = adata[myfilter5].obsm['X_umap']

# # # # Round UMAP to find overlaps
# # # umap_coords = pd.DataFrame(adata_filtered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
# # # umap_coords_rounded = umap_coords.round(2)
# # # adata_filtered.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# # # # Group by rounded coordinates and source_group
# # # group_counts = adata_filtered.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)

# # # Define overlapping coordinate keys (both PBMC and Glioblastoma present)
# # overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['glioblastoma'] > 0)].index

# # # Assign color labels for visualization
# # def assign_color(row):
# #     if row['umap_rounded'] in overlapping_keys:
# #         return 'Overlapping'
# #     return row['source_group']

# # adata_filtered.obs['color_group'] = adata_filtered.obs.apply(assign_color, axis=1)

# # # Define color map
# color_dict = {'PBMC': 'blue', 'glioblastoma': 'red', 'Overlapping': 'green'}

# # Plot UMAP with overlaps highlighted
# sc.pl.umap(adata_filtered, color='color_group', palette=color_dict, title='UMAP: PBMC vs Glioblastoma (Overlap in Green)')

# # --- Calculate Jaccard index ---

# # Get sets of unique rounded UMAP coordinates for each group
# pbmc_coords = set(group_counts[group_counts['PBMC'] > 0].index)
# glioblastoma_coords = set(group_counts[group_counts['glioblastoma'] > 0].index)

# # Calculate intersection and union
# intersection = pbmc_coords.intersection(glioblastoma_coords)
# union = pbmc_coords.union(glioblastoma_coords)

# # Calculate Jaccard index
# jaccard_index = len(intersection) / len(union) if len(union) > 0 else np.nan

# print(f"Jaccard Index for PBMC vs Glioblastoma UMAP overlaps: {jaccard_index:.4f}")
