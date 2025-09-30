# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 19:24:36 2025

@author: AROCHE6
"""

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gseapy as gp

# -------------------------------
# Load data
# -------------------------------
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)

# -------------------------------
# Define TGF-beta genes
# -------------------------------
tgfb_genes = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]

# tgfb sig genes = [
#     "JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1", "CTNNB1", "IFNG",
#     "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300", "SMURF2", "TFRC", "SPTBN1", "ATF3",
#     "MAPK14", "MAPK13", "APC", "MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111",
#     "MAPK9", "ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12"
# ]
# ranked rig tgf b ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]

tgfb_genes = [g for g in tgfb_genes if g in adata.var_names]
print(f"Using {len(tgfb_genes)} genes")

# -------------------------------
# Subset data
# -------------------------------
nk_gbm = adata[(adata.obs["source"] == "glioblastoma")].copy()
nk_pbmc = adata[(adata.obs["source"] == "PBMC")].copy()
# New mask: count how many genes each cell exceeds PBMC mean for
exceed_matrix = []

for gene in tgfb_genes:
    pbmc_expr = nk_pbmc[:, gene].X
    pbmc_mean = pbmc_expr.toarray().mean() if sp.issparse(pbmc_expr) else np.mean(pbmc_expr)

    gbm_expr = nk_gbm[:, gene].X
    gbm_expr = gbm_expr.toarray().flatten() if sp.issparse(gbm_expr) else np.asarray(gbm_expr).flatten()

    exceed_matrix.append(gbm_expr > pbmc_mean)

# Combine into array and count per cell
exceed_array = np.stack(exceed_matrix, axis=1)  # shape: (n_cells, n_genes)
exceed_count = exceed_array.sum(axis=1)

# Label as High if >= 50% genes above PBMC mean
threshold = int(len(tgfb_genes) * 0.5)
nk_gbm.obs["TGFB_High_vs_PBMC"] = np.where(exceed_count >= threshold, "High", "Other")

# # -------------------------------
# # Build filter for "High" GBM NK cells
# # -------------------------------
# cell_mask = np.ones(nk_gbm.shape[0], dtype=bool)

# for gene in tgfb_genes:
#     pbmc_expr = nk_pbmc[:, gene].X
#     pbmc_mean = pbmc_expr.toarray().mean() if sp.issparse(pbmc_expr) else np.mean(pbmc_expr)

#     gbm_expr = nk_gbm[:, gene].X
#     gbm_expr = gbm_expr.toarray().flatten() if sp.issparse(gbm_expr) else np.asarray(gbm_expr).flatten()

#     cell_mask &= gbm_expr > pbmc_mean

# nk_gbm.obs["TGFB_High_vs_PBMC"] = np.where(cell_mask, "High", "Other")
# print(nk_gbm.obs["TGFB_High_vs_PBMC"].value_counts())

# # -------------------------------
# Differential Expression
# -------------------------------
sc.tl.rank_genes_groups(nk_gbm, groupby="TGFB_High_vs_PBMC", reference="Other", method="wilcoxon")
sc.pl.rank_genes_groups(nk_gbm, n_genes=20, sharey=False)

de_df = sc.get.rank_genes_groups_df(nk_gbm, group="High")

# -------------------------------
# GSEA on DE genes
# -------------------------------
rnk = de_df[["names", "logfoldchanges"]].sort_values("logfoldchanges", ascending=False)

pre_res = gp.prerank(
    rnk=rnk,
    gene_sets="KEGG_2016",
    outdir="gsea_results_tgfb_vs_pbmc",
    min_size=15,
    max_size=500,
    permutation_num=100,
    seed=42
)

# -------------------------------
# Save Results
# -------------------------------
output_path = r"C:\Users\aroche6\abc\12_08_vs_pbmc.xlsx"
with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
    nk_gbm.obs[["TGFB_High_vs_PBMC"]].to_excel(writer, sheet_name="Cell_Labels")
    de_df.to_excel(writer, sheet_name="DE_High_vs_Other", index=False)
    rnk.to_excel(writer, sheet_name="GSEA_Rank_List", index=False)
    pre_res.res2d.reset_index().rename(columns={"index": "Pathway"}).to_excel(writer, sheet_name="GSEA_Results", index=False)

print(f"✅ Analysis complete. Results saved to:\n{output_path}")

