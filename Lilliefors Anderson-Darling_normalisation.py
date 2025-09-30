# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 01:40:23 2025

@author: AROCHE6
"""

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
from scipy.stats import mannwhitneyu, shapiro, probplot, zscore, anderson
from statsmodels.stats.diagnostic import lilliefors  # <-- NEW

# Load the AnnData object
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = ad.read_h5ad(file_path)

# Define the TGFB genes of interest (unchanged) ...
TGFB_genes = ["ACVR1B","ACVR1C","ACVR2A","ACVR2B","AMH","ASCL2","ATF1","ATF3","ATP1B2","ATP1B3","BAMBI",
              "BMP2","BMP7","BMP8B","BMPR1A","BMPR1B","BMPR2","CALD1","CD40","CDK8","CDKN1B","CDKN2B",
              "CREBBP","CTSC","CUL1","EIPR1","E2F4","E2F5","EP300","FST","ID1","ID2","ID3","ID4","INHBA",
              "LRRC2","LTBP1","NBL1","RBX1","SKP1","SMAD1","SMAD2","SMAD3","SMAD4","SMAD5","SMAD6",
              "SMAD7","SMAD9","SMURF1","SMURF2","SP1","TFDP1","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2",
              "TGIF1","TGIF2","VPS53","VPS54","VPS52","VPS51","ZFYVE16","ZFYVE9"]

# Filter for TGFB genes present in the dataset
genes_in_data = [gene for gene in TGFB_genes if gene in adata.var_names]
print(f"Filtered TGFB genes present in the data: {genes_in_data}")

# Subset and compute mean expression
adata_tgfb = adata[:, genes_in_data]
adata_tgfb.obs['mean_expression'] = np.asarray(adata_tgfb.X.mean(axis=1)).ravel()

# Optional: keep your min-max scaled version for visualization if you like
adata_tgfb.obs['mean_expression_scaled'] = minmax_scale(adata_tgfb.obs['mean_expression'])

# >>> For normality testing, use an unbounded transform (log1p), then z-score:
adata_tgfb.obs['mean_expression_log1p'] = np.log1p(adata_tgfb.obs['mean_expression'])
adata_tgfb.obs['mean_expression_z'] = zscore(adata_tgfb.obs['mean_expression_log1p'], nan_policy='omit')

# Extract cohorts
glioblastoma_cells = adata_tgfb[adata_tgfb.obs['source'] == 'glioblastoma']
pbmc_cells        = adata_tgfb[adata_tgfb.obs['source'] == 'PBMC']

glioblastoma_scores_z = np.asarray(glioblastoma_cells.obs['mean_expression_z'].dropna())
pbmc_scores_z         = np.asarray(pbmc_cells.obs['mean_expression_z'].dropna())

# ---- Normality tests ----
def shapiro_safe(x, max_n=5000, random_state=0):
    """Shapiro is only validated up to ~5000; subsample if larger."""
    x = np.asarray(x)
    if x.size > max_n:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(x.size, size=max_n, replace=False)
        x = x[idx]
    return shapiro(x)

# Shapiro–Wilk
shapiro_glioblastoma = shapiro_safe(glioblastoma_scores_z)
shapiro_pbmc         = shapiro_safe(pbmc_scores_z)

# Lilliefors (KS with mu, sigma estimated from data; proper null)
lf_gliostat, lf_gliop = lilliefors(glioblastoma_scores_z, dist='norm')
lf_pbmcstat,  lf_pbmcp = lilliefors(pbmc_scores_z, dist='norm')

# Anderson–Darling (more sensitive in tails)
ad_gli = anderson(glioblastoma_scores_z, dist='norm')
ad_pbm = anderson(pbmc_scores_z, dist='norm')

print("\n=== Normality tests on log1p z-scored TGFB mean expression ===")
print(f"Shapiro–Wilk (glioblastoma): W={shapiro_glioblastoma.statistic:.4f}, p={shapiro_glioblastoma.pvalue:.3e}")
print(f"Shapiro–Wilk (PBMC)        : W={shapiro_pbmc.statistic:.4f}, p={shapiro_pbmc.pvalue:.3e}")

print(f"Lilliefors (glioblastoma)  : KS*={lf_gliostat:.4f}, p={lf_gliop:.3e}")
print(f"Lilliefors (PBMC)          : KS*={lf_pbmcstat:.4f}, p={lf_pbmcp:.3e}")

# Anderson returns a statistic and critical values for several alpha levels
alphas = [15, 10, 5, 2.5, 1]  # %
print(f"Anderson–Darling (glioblastoma): A2={ad_gli.statistic:.4f}; crit={list(zip(alphas, np.round(ad_gli.critical_values,4)))}")
print(f"Anderson–Darling (PBMC)        : A2={ad_pbm.statistic:.4f}; crit={list(zip(alphas, np.round(ad_pbm.critical_values,4)))}")

# ---- Visuals (use z-scores for diagnostics) ----
plt.figure(figsize=(12, 6))
sns.kdeplot(glioblastoma_scores_z, label='Glioblastoma', fill=True)
sns.kdeplot(pbmc_scores_z,         label='PBMC',        fill=True)
plt.xlabel('log1p mean TGFB (z-score)')
plt.ylabel('Density')
plt.title('Distribution of TGFB gene scores (log1p z-scored)')
plt.legend()
plt.show()

# QQ-plots
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
probplot(glioblastoma_scores_z, dist="norm", plot=plt)
plt.title('QQ-plot: Glioblastoma (z)')

plt.subplot(1, 2, 2)
probplot(pbmc_scores_z, dist="norm", plot=plt)
plt.title('QQ-plot: PBMC (z)')
plt.tight_layout()
plt.show()
