# -*- coding: utf-8 -*-
"""
Created on Mon May 26 14:59:15 2025

@author: AROCHE6
"""

import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import gseapy as gp
import matplotlib.pyplot as plt
from pathlib import Path

# ---------------------------------------------------
# 1.  Load the data
# ---------------------------------------------------
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
# TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]
#TGFBSIGgene_list = [
#     "JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1",
#     "CTNNB1", "IFNG", "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300",
#     "SMURF2", "TFRC", "SPTBN1", "ATF3", "MAPK14", "MAPK13", "APC",
#     "MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111", "MAPK9",
#     "ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12"
# ]
# ---------------------------------------------------
# 2.  Create a High / Low label for TGFB1
# ---------------------------------------------------
gene = "TGFB1"

x = adata[:, gene].X            # (n_cells, 1) sparse or dense

# Convert to dense 1-D array
if sp.issparse(x):
    expr = x.toarray().ravel()  # <- use .toarray() instead of .A1
else:
    expr = np.asarray(x).ravel()

# Choose threshold
threshold = np.median(expr)     # or set manually, e.g. threshold = 1.5
adata.obs["TGFB1_status"] = np.where(expr > threshold, "High", "Low")

print(f"TGFB1 threshold = {threshold:.3f}")
print(adata.obs["TGFB1_status"].value_counts(), "\n")

# ---------------------------------------------------
# 3.  Violin plot of another gene split by TGFB1 status
# ---------------------------------------------------
sc.pl.violin(
    adata,
    keys="NFKB1",
    groupby="TGFB1_status",
    stripplot=True
)

# ---------------------------------------------------
# 4.  Differential expression: High vs Low (Wilcoxon)
# ---------------------------------------------------
sc.tl.rank_genes_groups(
    adata,
    groupby="TGFB1_status",
    reference="Low",
    method="wilcoxon"
)
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

de_df = sc.get.rank_genes_groups_df(adata, group="High")

# ---------------------------------------------------
# 5.  Pre-ranked GSEA
# ---------------------------------------------------
rnk = de_df[["names", "logfoldchanges"]].sort_values(
    "logfoldchanges",
    ascending=False
)

pre_res = gp.prerank(
    rnk=rnk,
    gene_sets="KEGG_2016",      # change to any library you like
    outdir="gsea_results",
    min_size=15,
    max_size=500,
    permutation_num=100,
    seed=42
)

print("\nTop enriched pathways:")
print(pre_res.res2d.head())

import os

# Set output path
output_path = r"C:\Users\aroche6\abc\tgfb1_analysis_results.xlsx"

# Create a writer to save multiple sheets
with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
    # Save TGFB1 labels + key metadata
    obs_cols = ["TGFB1_status"]
    if "cell_type" in adata.obs:
        obs_cols.append("cell_type")
    adata.obs[obs_cols].to_excel(writer, sheet_name="Cell_Labels")

    # Save differential expression results
    de_df.to_excel(writer, sheet_name="DE_High_vs_Low", index=False)

    # Save GSEA ranked gene list
    rnk.to_excel(writer, sheet_name="GSEA_Rank_List", index=False)

    # Save top enriched GSEA pathways
    top_gsea = pre_res.res2d.reset_index().rename(columns={"index": "Pathway"})
    top_gsea.to_excel(writer, sheet_name="GSEA_Results", index=False)

print(f"✅ Results exported to:\n{output_path}")
