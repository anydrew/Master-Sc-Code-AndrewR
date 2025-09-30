# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 11:50:46 2024

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

file_path = r"C:\Users\aroche6\abc\all_nk_cells.h5ad"
sc.read_h5ad(file_path)
adata = sc.read_h5ad(file_path)
# Create a plot showing genes detected as a function of UMI counts.

fig, ax = plt.subplots(figsize=(10, 7))

x = np.asarray(adata.X.sum(axis=1))[:, 0]
y = np.asarray(np.sum(adata.X > 0, axis=1))[:, 0]

ax.scatter(x, y, color="green", alpha=0.25)
ax.set_xlabel("UMI Counts")
ax.set_ylabel("Genes Detected")
ax.set_xscale('log')
ax.set_yscale('log', nonpositive='clip') # Corrected argument name

ax.set_xlim((0.5, 4500))
ax.set_ylim((0.5, 2000))

plt.show()

#@title Threshold cells according to knee plot { run: "auto", vertical-output: true }
cutoff = 200#@param {type:"integer"}
knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
cell_set = np.arange(len(knee))
num_cells = cell_set[knee > cutoff][::-1][0]

fig, ax = plt.subplots(figsize=(10, 7))


ax.loglog(knee, cell_set, linewidth=5, color="g")
ax.axvline(x=cutoff, linewidth=3, color="k")


ax.axhline(y=num_cells, linewidth=3, color="k")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
plt.show()

print(f"{num_cells:,.0f} cells passed the {cutoff} UMI threshold")

adata

# Filter the cells according to the threshold determined from the knee plot
sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_cells(adata, min_counts=knee[expected_num_cells])

# sc.pp.filter_cells(adata, min_counts=knee[expected_num_cells])


cutoff = 200#@param {type:"integer"}
knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
cell_set = np.arange(len(knee))
num_cells = cell_set[knee > cutoff][::-1][0]
    
fig, ax = plt.subplots(figsize=(10, 7))
    
    
ax.loglog(knee, cell_set, linewidth=5, color="g")
ax.axvline(x=cutoff, linewidth=3, color="k")


ax.axhline(y=num_cells, linewidth=3, color="k")

ax.set_xlabel("UMI Counts")
ax.set_ylabel("Set of Barcodes")

plt.grid(True, which="both")
plt.show()

# mito_ensembl_ids = sc.queries.mitochondrial_genes("mmusculus", attrname="ensembl_gene_id")


# These packages are pre-installed on Google Colab, but are included here to facilitate running this notebook locally
# !pip install --quiet matplotlib
# !pip install --quiet scikit-learn
# !pip install --quiet numpy
# !pip install --quiet scipy

# %%time
#     # `kb` is a wrapper for the kallisto and bustools program, and the kb-python package contains the kallisto and bustools executables.
    # !pip install --quiet kb-python==0.24.1
    
    

# %%time
#     # Install scanpy and other packages needed for single-cell RNA-seq analysis
#!pip install --quiet scanpy python-igraph louvain MulticoreTSNE pybiomart
    
    

# mito_genes = mito_ensembl_ids["ensembl_gene_id"].values
#     # for each cell compute fraction of counts in mito genes vs. all genes
#     # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
# adata.obs['percent_mito'] = np.sum(
# adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
#     # add the total counts per cell as observations-annotation to adata
# adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# mito_genes = mito_ensembl_ids["ensembl_gene_id"].values



# sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
# sc.pp.log1p(adata)

# sc.pl.highest_expr_genes(adata, n_top=20)

# Removes cells with less than 1070 umi counts
adata = adata[np.asarray(adata.X.sum(axis=1)).reshape(-1) > 1070]
# Removes genes with 0 umi counts
adata = adata[:, np.asarray(adata.X.sum(axis=0)).reshape(-1) > 0]

adata


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# adata.obs['n_genes'] = number

mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

adata

sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True)

#examine mitochondrial content 
sc.pl.scatter(adata, x='n_counts', y='percent_mito')
sc.pl.scatter(adata, x='n_counts', y='n_genes')

# Create a mask to filter out cells with more than 6500 genes, less than 200 genes or less than 0.2 mitochondrial umi counts
mask = np.logical_or((adata.obs.n_genes < 6500).values, (adata.obs.n_genes > 200).values, (adata.obs.percent_mito < 0.2).values)

#filter
adata = adata[mask, :]

adata

glyc_marker_genes = ['ACSS1', 'ACSS2', 'ACTG1', 'ADH5', 'ADPGK', 'AKR1A1', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH3B1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALDOC', 'BNIP3', 'BPGM', 'CELSR2', 'CELSR3', 'CETP', 'CNP', 'DLAT', 'DLD', 'DMPK', 'DMWD', 'EFNA4', 'EFNA5', 'EIF2D', 'ENO1', 'ENO2', 'ENO3', 'FBP1', 'G6PC3', 'GALM', 'GAPDH', 'GPI', 'HK1', 'HK2', 'HK3', 'HKDC1', 'LDHA', 'LDHB', 'MEGF6', 'MEGF8', 'MINPP1', 'PCK2', 'PDHA1', 'PDHB', 'PFKL', 'PFKM', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGM1', 'PGM2', 'PKM', 'TPI1']

# ax = sc.pl.stacked_violin(adata, marker_genes, groupby='subset_source', rotation=90)

# ax = sc.pl.stacked_violin(adata, marker_genes, groupby='subset_source', rotation=90)

# ax = sc.pl.stacked_violin(adata, marker_genes, groupby='subset_source', rotation=90)
myfilter2 = (adata.obs['source'] == 'brain_normal') | (adata.obs['source'] == 'glioblastoma_tumor')
ax = sc.pl.dotplot(adata[myfilter2], glyc_marker_genes, groupby='subset')
# brain_normal
# glioblastoma_tumor
glyc_marker_genes = ['ACSS1', 'ACSS2', 'ACTG1', 'ADH5', 'ADPGK', 'AKR1A1', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH3B1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALDOC', 'BNIP3', 'BPGM', 'CELSR2', 'CELSR3', 'CETP', 'CNP', 'DLAT', 'DLD', 'DMPK', 'DMWD', 'EFNA4', 'EFNA5', 'EIF2D', 'ENO1', 'ENO2', 'ENO3', 'FBP1', 'G6PC3', 'GALM', 'GAPDH', 'GPI', 'HK1', 'HK2', 'HK3', 'HKDC1', 'LDHA', 'LDHB', 'MEGF6', 'MEGF8', 'MINPP1', 'PCK2', 'PDHA1', 'PDHB', 'PFKL', 'PFKM', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGM1', 'PGM2', 'PKM', 'TPI1']

myfilter3 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56dim, brain_normal') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
#ax = sc.pl.dotplot(adata[myfilter3], glyc_marker_genes, groupby='subset')
ax = sc.pl.dotplot(adata[myfilter3], glyc_marker_genes, groupby='subse_source,', title="Glycolysis and Gluconeogenesis Expression in Glioblastoma NK Cells vs Controls")
NK_ION_GENES = ['HVCN1', 'GPR68', 'NPEPPS', 'SLC16A3', 'SLC9A3R1', 'CLIC4', 'CLIC3', 'SLC9A1', 'SLC12A9', 'SLC16A7', 'SLC4A1', 'SLC4A7']
# 'KIR+, PBMC'
# 'CD57+, PBMC'
# 'NKG2A+, PBMC'
# 'CD56bright, PBMC'
# 'CD56dim, PBMC'
# 'Adaptive, PBMC'
# 'CD56dim, brain_normal'
# 'CD56bright, glioblastoma_tumor'
# 'CD56dim, glioblastoma_tumor'

# import scipy.stats as stats

# # Assuming 'subset' is the column with different groups
# group_labels = adata.obs['subset'].unique()

# for gene in glyc_marker_genes:
#     plt.figure(figsize=(8, 6))
#     for label in group_labels:
#         subset_data = adata[myfilter2 & (adata.obs['subset'] == label)]
#         expression_values = subset_data.obs_vector(gene)

#         plt.scatter([label] * len(expression_values), expression_values, label=label)

#     # Perform t-test between groups
#     p_value = stats.f_oneway(*[adata[myfilter2 & (adata.obs['subset'] == label)].obs_vector(gene) for label in group_labels]).pvalue

#     # Add significance marker
#     if p_value < 0.05:
#         plt.text(len(group_labels) // 2, max(expression_values) + 1, f'p-value: {p_value:.4f}', ha='center')

#     plt.title(f'Dotplot for {gene}')
#     plt.xlabel('Group')
#     plt.ylabel('Gene Expression')
#     plt.legend()
#     plt.show()

# ax = sc.pl.dotplot(adata[myfilter3], glyc_marker_genes, groupby='subset', title="Glycolysis and Gluconeogenesis Expression in Glioblastoma NK Cells vs Controls")

# #### statistics test for each gene
# import scipy.stats as stats
# import scanpy as sc
# import seaborn as sns
# import matplotlib.pyplot as plt
# import pandas as pd, numpy as np, matplotlib.pyplot as plt, scipy.stats as ss, math
# import scipy.spatial.distance as distance
# from sklearn.preprocessing import normalize
# from umap import UMAP
# from sklearn.manifold import TSNE, SpectralEmbedding

# file_path = r"C:\Users\aroche6\abc\all_nk_cells.h5ad"
# sc.read_h5ad(file_path)
# adata = sc.read_h5ad(file_path)

# glyc_marker_genes = ['ACSS1', 'ACSS2', 'ACTG1', 'ADH5', 'ADPGK', 'AKR1A1', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH3B1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALDOC', 'BNIP3', 'BPGM', 'CELSR2', 'CELSR3', 'CETP', 'CNP', 'DLAT', 'DLD', 'DMPK', 'DMWD', 'EFNA4', 'EFNA5', 'EIF2D', 'ENO1', 'ENO2', 'ENO3', 'FBP1', 'G6PC3', 'GALM', 'GAPDH', 'GPI', 'HK1', 'HK2', 'HK3', 'HKDC1', 'LDHA', 'LDHB', 'MEGF6', 'MEGF8', 'MINPP1', 'PCK2', 'PDHA1', 'PDHB', 'PFKL', 'PFKM', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGM1', 'PGM2', 'PKM', 'TPI1']

# myfilter3 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56dim, brain_normal') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')

# # Assuming 'subset' is the column with different groups
# group_labels = adata.obs['subset'].unique()

# for gene in glyc_marker_genes:
#     plt.figure(figsize=(30, 6))
#     for label in group_labels:
#         subset_data = adata[myfilter3 & (adata.obs['subset'] == label)]
#         expression_values = subset_data.obs_vector(gene)

#         plt.scatter([label] * len(expression_values), expression_values, label=label)

#     # Perform t-test between groups
#     p_value = stats.f_oneway(*[adata[myfilter3 & (adata.obs['subset'] == label)].obs_vector(gene) for label in group_labels]).pvalue

#     # Add significance marker
#     if p_value < 0.05:
#         plt.text(len(group_labels) // 2, max(expression_values) + 1, f'p-value: {p_value:.4f}', ha='center')

#     plt.title(f'Dotplot for {gene}')
#     plt.xlabel('Group')
#     plt.ylabel('Gene Expression')
#     #ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
#     #ax.figure
#     #plt.legend()
#     plt.show()

# ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
# ax.figure

# ###same but place all on same image
# ncols = 4  # Number of columns in the subplot grid
# nrows = math.ceil(len(glyc_marker_genes) / ncols)

# fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4 * nrows))

# for i, gene in enumerate(glyc_marker_genes):
#     ax = axes.flatten()[i] if len(glyc_marker_genes) > 1 else axes  # Adjust if there's only one gene

#     for label in group_labels:
#         subset_data = adata[myfilter3 & (adata.obs['subset'] == label)]
#         expression_values = subset_data.obs_vector(gene)

#         ax.scatter([label] * len(expression_values), expression_values, label=label)

#     # Perform t-test between groups
#     p_value = stats.f_oneway(*[adata[myfilter3 & (adata.obs['subset'] == label)].obs_vector(gene) for label in group_labels]).pvalue

#     # Add significance marker
#     if p_value < 0.05:
#         ax.text(len(group_labels) // 2, max(expression_values) + 1, f'p-value: {p_value:.4f}', ha='center')

#     ax.set_title(f'Dotplot for {gene}')
#     ax.set_xlabel('Group')
#     ax.set_ylabel('Gene Expression')
#     ax.legend()

# plt.tight_layout()
# plt.show()

# #from statsmodels.stats.multicomp import pairwise_tukeyhsd

# # Get all expression values for each group
# data = [adata[myfilter3 & (adata.obs['subset'] == label)].obs_vector(gene) for label in group_labels]

# # Flatten the data and create corresponding labels
# flat_data = np.concatenate(data)
# labels = np.concatenate([[label] * len(data_point) for label, data_point in zip(group_labels, data)])

# # Perform Tukey-Kramer post hoc test
# tukey_results = pairwise_tukeyhsd(flat_data, labels)

# # Print the summary of the test
# print(tukey_results.summary())
