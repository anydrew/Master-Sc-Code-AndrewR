# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 13:36:49 2024

@author: AROCHE6
"""

import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import minmax_scale
from scipy.stats import mannwhitneyu, shapiro, kstest, norm, probplot

# Load the AnnData object
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = ad.read_h5ad(file_path)

# Define the TGFB genes of interest
TGFB_genes = ["ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "AMH", "ASCL2", "ATF1", "ATF3", "ATP1B2", "ATP1B3", "BAMBI",
              "BMP2", "BMP7", "BMP8B", "BMPR1A", "BMPR1B", "BMPR2", "CALD1", "CD40", "CDK8", "CDKN1B", "CDKN2B",
              "CREBBP", "CTSC", "CUL1", "EIPR1", "E2F4", "E2F5", "EP300", "FST", "ID1", "ID2", "ID3", "ID4", "INHBA",
              "LRRC2", "LTBP1", "NBL1", "RBX1", "SKP1", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6",
              "SMAD7", "SMAD9", "SMURF1", "SMURF2", "SP1", "TFDP1", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2",
              "TGIF1", "TGIF2", "VPS53", "VPS54", "VPS52", "VPS51", "ZFYVE16", "ZFYVE9"]

# Filter for TGFB genes present in the dataset
genes_in_data = [gene for gene in TGFB_genes if gene in adata.var_names]
print(f"Filtered TGFB genes present in the data: {genes_in_data}")

# Subset the data for these genes
adata_tgfb = adata[:, genes_in_data]

# Calculate the mean expression values for each cell
adata_tgfb.obs['mean_expression'] = adata_tgfb.X.mean(axis=1)

# Normalize the scores
adata_tgfb.obs['mean_expression_scaled'] = minmax_scale(adata_tgfb.obs['mean_expression'])

# Extract the glioblastoma and PBMC cells
glioblastoma_cells = adata_tgfb[adata_tgfb.obs['source'] == 'glioblastoma']
pbmc_cells = adata_tgfb[adata_tgfb.obs['source'] == 'PBMC']

# Perform normality tests
glioblastoma_scores = glioblastoma_cells.obs['mean_expression_scaled']
pbmc_scores = pbmc_cells.obs['mean_expression_scaled']

shapiro_glioblastoma = shapiro(glioblastoma_scores)
shapiro_pbmc = shapiro(pbmc_scores)
ks_glioblastoma = kstest(glioblastoma_scores, 'norm', args=(glioblastoma_scores.mean(), glioblastoma_scores.std()))
ks_pbmc = kstest(pbmc_scores, 'norm', args=(pbmc_scores.mean(), pbmc_scores.std()))

print(f'Shapiro-Wilk test for glioblastoma: statistic={shapiro_glioblastoma.statistic}, p-value={shapiro_glioblastoma.pvalue}')
print(f'Shapiro-Wilk test for PBMC: statistic={shapiro_pbmc.statistic}, p-value={shapiro_pbmc.pvalue}')
print(f'Kolmogorov-Smirnov test for glioblastoma: statistic={ks_glioblastoma.statistic}, p-value={ks_glioblastoma.pvalue}')
print(f'Kolmogorov-Smirnov test for PBMC: statistic={ks_pbmc.statistic}, p-value={ks_pbmc.pvalue}')

# Plot the distribution of scores for glioblastoma and PBMC cells
plt.figure(figsize=(12, 6))
sns.kdeplot(glioblastoma_scores, label='Glioblastoma', fill=True)
sns.kdeplot(pbmc_scores, label='PBMC', fill=True)
plt.xlabel('Scaled Mean Expression')
plt.ylabel('Density')
plt.title('Distribution of TGFB Gene Scores')
plt.legend()
plt.show()

# QQ-plots for normality
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
probplot(glioblastoma_scores, dist="norm", plot=plt)
plt.title('QQ-plot for Glioblastoma')

plt.subplot(1, 2, 2)
probplot(pbmc_scores, dist="norm", plot=plt)
plt.title('QQ-plot for PBMC')

plt.show()
