# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 20:35:39 2024

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
# file_path = r"C:\Users\aroche6\abc\all_nk_cells.h5ad"

# file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
# adata = sc.read_h5ad(file_path)
# myfilter6 = (
#     (adata.obs['subset_source'] == 'CD56bright_PBMC') |
#     (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
# )
# filtered_adata = adata[myfilter6]
# sc.pl.umap(filtered_adata, color='subset_source')
# filtered_adata.uns['subset_source_colors'] = ['red', 'blue']  # Make sure order aligns with subset categories

# # Step 4: Plot UMAP
# sc.pl.umap(filtered_adata, color='subset_source')

# import scanpy as sc
# import pandas as pd

# # Step 1: Filter data
# myfilter6 = (
#     (adata.obs['subset_source'] == 'CD56bright_PBMC') |
#     (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma')
# )
# adata_bright = adata[myfilter6].copy()

# # Step 2: Label source
# adata_bright.obs['source_group'] = adata_bright.obs['subset_source'].apply(
#     lambda x: 'PBMC' if x == 'CD56bright_PBMC' else 'Glioblastoma'
# )

# # Step 3: Reuse UMAP from full data if available
# adata_bright.obsm['X_umap'] = adata[myfilter6].obsm['X_umap']

# # Step 4: Detect overlapping cells (based on rounded UMAP coordinates)
# umap_coords = pd.DataFrame(adata_bright.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
# umap_coords_rounded = umap_coords.round(2)
# adata_bright.obs['umap_rounded'] = list(umap_coords_rounded.itertuples(index=False, name=None))

# # Step 5: Find overlapping coordinate keys
# group_counts = adata_bright.obs.groupby(['umap_rounded', 'source_group']).size().unstack(fill_value=0)
# overlapping_keys = group_counts[(group_counts['PBMC'] > 0) & (group_counts['Glioblastoma'] > 0)].index

# # Step 6: Assign color group
# def assign_color(row):
#     if row['umap_rounded'] in overlapping_keys:
#         return 'Overlapping'
#     return row['source_group']

# adata_bright.obs['color_group'] = adata_bright.obs.apply(assign_color, axis=1)

# # Step 7: Plot 1 – All cells
# color_dict = {'PBMC': 'red', 'Glioblastoma': 'blue', 'Overlapping': 'gold'}
# sc.pl.umap(
#     adata_bright,
#     color='color_group',
#     palette=color_dict,
#     title='CD56bright PBMC vs Glioblastoma (Overlap = Gold)'
# )

# # Step 8: Plot 2 – Just overlapping cells
# overlapping_only = adata_bright[adata_bright.obs['color_group'] == 'Overlapping'].copy()
# sc.pl.umap(
#     overlapping_only,
#     color='color_group',
#     palette={'Overlapping': 'gold'},
#     title='CD56bright Overlapping Cells Only (Gold)'
# )

adata = sc.read_h5ad(file_path)
# sc.pl.umap(adata, color=["subset_source"])
# Step 1: Define your filter
# myfilter7 = (
#     (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') |
#     (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor') | 
#     (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
# )

myfilter7 = ((adata.obs['subset_source'] == 'Adaptive_PBMC') | (adata.obs['subset_source'] == 'CD57+_PBMC') | (adata.obs['subset_source'] == 'KIR+_PBMC') | (adata.obs['subset_source'] == 'NKG2A+_PBMC') | (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma'))
filtered_adata = adata[myfilter7]
filtered_adata.uns['subset_source_colors'] = ['red', 'red', 'red', 'blue']  # Make sure order aligns with subset categories

# Step 2: Ensure 'subset' is categorical
# filtered_adata.obs['subset'] = filtered_adata.obs['subset'].astype('category')

# Step 3: Set colors manually
# filtered_adata.uns['subset_colors'] = ['red', 'blue']  # Make sure order aligns with subset categories

# Step 4: Plot UMAP
sc.pl.umap(filtered_adata, color='subset_source')

# import scanpy as sc
# import matplotlib.pyplot as plt
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
# # Define the PBMC and glioblastoma sources
# pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
# gbm_source = ['TiCD56dim_glioblastoma']

# # Create a filtered AnnData object with relevant cells
# myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
# adata_filtered = adata[myfilter7].copy()

# # Create a new column to define the group (PBMC or GBM)
# adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
#     lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
# )

# # Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)
# sc.pp.pca(adata_filtered)
# sc.pp.neighbors(adata_filtered)
# sc.tl.umap(adata_filtered)

# # Plot UMAP colored by group
# sc.pl.umap(adata_filtered, color='source_group', palette=['blue', 'red'], title='UMAP: PBMC vs Glioblastoma dim cells')
# C:\Users\aroche6\AppData\Roaming\Python\Python311\site-packages\scanpy\plotting\_tools\scatterplots.py:394: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
#   cax = scatter(

# import scanpy as sc
# import matplotlib.pyplot as plt

# # Define the PBMC and glioblastoma sources
# pbmc_sources = ['KIR+_PBMC', 'Adaptive_PBMC', 'CD57+_PBMC', 'NKG2A+_PBMC']
# gbm_source = ['TiCD56dim_glioblastoma']

# # Create a filtered AnnData object with relevant cells
# myfilter7 = adata.obs['subset_source'].isin(pbmc_sources + gbm_source)
# adata_filtered = adata[myfilter7].copy()

# # Create a new column to define the group (PBMC or GBM)
# adata_filtered.obs['source_group'] = adata_filtered.obs['subset_source'].apply(
#     lambda x: 'PBMC' if x in pbmc_sources else 'Glioblastoma'
# )

# # Run UMAP (assuming you already did PCA and neighbors; otherwise, do them)

# # Plot UMAP colored by group
# sc.pl.umap(adata_filtered, color='source_group', palette=['blue', 'red'], title='UMAP: PBMC vs Glioblastoma dim cells')




# myfilter12 = (adata.obs['source'] == 'PBMC') | (adata.obs['source'] == 'glioblastoma')
# # sc.pl.umap(adata, color=["cell_type", 'source'])
# myfilter4 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC')
# #myfilter5 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
# myfilter6 = (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
# # filtered_adata = adata[myfilter4]
# # filtered_adata = adata[myfilter5]

# filtered_adata = adata[myfilter12]

# # # Generate UMAP for the filtered data
# sc.tl.umap(filtered_adata)  # Run UMAP on the filtered data
# # sc.pl.umap(filtered_adata, color=['source', "cell_type", 'subset', 'dataset'])

# # myfilter6 = (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
# # filtered_adata = adata[myfilter6]
# # sc.pl.umap(filtered_adata, color=['subset'])


# # Step 1: Define your filter
# # myfilter6 = (
# #     (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') |
# #     (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
# # # )
# # filtered_adata = adata[myfilter6]

# # Step 2: Ensure 'subset' is categorical
# filtered_adata.obs['source'] = filtered_adata.obs['source'].astype('category')

# # Step 3: Set colors manually
# filtered_adata.uns['source'] = ['red', 'blue']  # Make sure order aligns with subset categories

# # Step 4: Plot UMAP
# sc.pl.umap(filtered_adata, color='source')

# #sc.tl.umap(filtered_adata)

# # myfilter10 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC')
# # filtered_adata = adata[myfilter4]
# # sc.pl.umap(filtered_adata, color=['nhood_group'])





# # # Cluster the data using the Leiden algorithm
# # sc.tl.leiden(filtered_adata)

# # # Generate UMAP plot with Leiden clusters
# # sc.pl.umap(filtered_adata, color=["leiden", "cell_type", 'source'])

# # # Perform differential expression analysis using the Wilcoxon test
# # sc.tl.rank_genes_groups(filtered_adata, 'leiden', method='wilcoxon')

# # # Extract the results and filter for adjusted p-values
# # ranked_genes = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['names'])
# # pvals_adj = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['pvals_adj'])

# # # Find the top 20 genes with the lowest adjusted p-values for each cluster
# # top_genes = {}
# # for cluster in ranked_genes.columns:
# #     top_genes[cluster] = ranked_genes[cluster][pvals_adj[cluster].argsort()][:20]

# # # Plot the top 20 markers based on adjusted p-values
# # sc.pl.rank_genes_groups(filtered_adata, n_genes=20, sharey=False)

# # ranked_genes = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['names'])
# # pvals_adj = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['pvals_adj'])

# # # Find the top 20 genes with the lowest adjusted p-values for each cluster
# # top_genes = {}
# # for cluster in ranked_genes.columns:
# #     top_genes[cluster] = ranked_genes[cluster][pvals_adj[cluster].argsort()][:20]

# # # Plot the top 20 markers based on adjusted p-values
# # sc.pl.rank_genes_groups(filtered_adata, n_genes=20, sharey=False)

# # sc.tl.rank_genes_groups(filtered_adata, 'source', method='wilcoxon')
# # # Extract the results and filter for adjusted p-values
# # ranked_genes = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['names'])
# # pvals_adj = pd.DataFrame(filtered_adata.uns['rank_genes_groups']['pvals_adj'])

# # # Find the top 20 genes with the lowest adjusted p-values for each cluster
# # top_genes = {}
# # for cluster in ranked_genes.columns:
# #     top_genes[cluster] = ranked_genes[cluster][pvals_adj[cluster].argsort()][:20]

# # # Plot the top 20 markers based on adjusted p-values
# # sc.pl.rank_genes_groups(filtered_adata, n_genes=20, sharey=False)