# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 20:03:28 2025

@author: AROCHE6
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 16:39:25 2025

@author: AROCHE6
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import ranksums
import seaborn as sns
import matplotlib.pyplot as plt
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
adata_gbm = adata[adata.obs['source'] == 'glioblastoma']
adata_pbmc = adata[adata.obs['source'] == 'PBMC']
gene_list = ["JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1",
"CTNNB1", "IFNG", "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300",
"SMURF2", "TFRC", "SPTBN1", "ATF3", "MAPK14", "MAPK13", "APC",
"MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111", "MAPK9",
"ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12"
]# TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]
#TGFBSIGgene_list = [
    # "JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1",
    # "CTNNB1", "IFNG", "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300",
    # "SMURF2", "TFRC", "SPTBN1", "ATF3", "MAPK14", "MAPK13", "APC",
    # "MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111", "MAPK9",
    # "ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12"
# ]
# TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]


results = []
for gene in gene_list:
    if gene in adata.var_names:
        expr_gbm = adata_gbm[:, gene].X.toarray().flatten()
        expr_pbmc = adata_pbmc[:, gene].X.toarray().flatten()

        stat, p_value = ranksums(expr_gbm, expr_pbmc)
        log2fc = np.log2(np.mean(expr_gbm) + 1) - np.log2(np.mean(expr_pbmc) + 1)

        results.append({'Gene': gene, 'Log2FC': log2fc, 'p_value': p_value})

df_results = pd.DataFrame(results)
df_results['Significant'] = (df_results['p_value'] < 0.05) & (df_results['Log2FC'] > 1)
significant_genes = df_results[df_results['Significant']]['Gene'].tolist()
gbm_high_expr_cells = adata_gbm[:, significant_genes].to_df().sum(axis=1) > 0
#all
#gbm_high_expr_cells = (adata_gbm[:, significant_genes].to_df() > 0).all(axis=1)

# gbm_low_expr_cells = adata_gbm[:, significant_genes].to_df().sum(axis=1) < 0
#all
#gbm_low_expr_cells = (adata_gbm[:, significant_genes].to_df() < 0).all(axis=1)

adata_gbm_high = adata_gbm[gbm_high_expr_cells]

gbm_low_expr_cells = ~gbm_high_expr_cells
adata_gbm_low = adata_gbm[gbm_low_expr_cells]
interferon_gamma_genes = [
    "ABL1", "AXL", "BCL3", "C1QBP", "CD2", "CD3E", "CD14", "CD47", "CEBPG", 
    "CCR7", "CR1", "DDIT3", "F2RL1", "GAS6", "GATA3", "HLA-A", "HLA-DPA1",
    "HLA-DPB1", "HLA-DRB1", "HMGB1", "HRAS", "HSPD1", "IRF8", "IL1B",
    "IL1R1", "IL2", "IL10", "IL12A", "IL12B", "IL12RB1", "IL12RB2",
    "IL18", "INHA", "INHBA", "ISL1", "JAK2", "LGALS9", "LTA", "PDE4B",
    "PDE4D", "PRNP", "RARA", "TRIM27", "XCL1", "SLAMF1", "SLC11A1",
    "TLR3", "TLR4", "TNF", "TNFSF4", "TXK", "TYK2", "SCGB1A1", "WNT5A",
    "ZP3", "LAPTM5", "FZD5", "SLC7A5", "RIPK2", "FADD", "IL18R1",
    "PGLYRP1", "IL1RL1", "IL27RA", "ISG15", "NR1H4", "RASGRP1", "EBI3",
    "CD96", "CD226", "LILRB1", "ARID5A", "LILRB4", "RIPK3", "BTN3A2",
    "BTN3A1", "CD160", "KLRK1", "SCRIB", "PTPN22", "IL36RN", "PYCARD",
    "CD274", "FOXP3", "TLR7", "TLR8", "IL23A", "CYRIB", "CD244", "IL20RB",
    "TLR9", "SASH3", "CRTAM", "HMHB1", "IL21", "VSIR", "NOD2", "CLEC7A",
    "ZC3H12A", "PDCD1LG2", "CD276", "HAVCR2", "PGLYRP2", "SLAMF6",
    "SIRPA", "IL23R", "ZFPM1", "NLRP6", "IL27", "IFNL1", "LGALS9B",
    "MIR24-1", "LGALS7B", "LGALS9C", "CCR2", "MIR708", "KLRC4-KLRK1"
]

interferon_results = []
for gene in interferon_gamma_genes:
    if gene in adata.var_names:
        expr_high = adata_gbm_high[:, gene].X.toarray().flatten()
        expr_low = adata_gbm_low[:, gene].X.toarray().flatten()

        stat, p_value = ranksums(expr_high, expr_low)
        log2fc = np.log2(np.mean(expr_high) + 1) - np.log2(np.mean(expr_low) + 1)

        interferon_results.append({'Gene': gene, 'Log2FC': log2fc, 'p_value': p_value})

df_interferon = pd.DataFrame(interferon_results)
df_interferon['Significant'] = (df_interferon['p_value'] < 0.05) & (df_interferon['Log2FC'] > 1)

output_path = r"C:\Users\aroche6\abc\2.0interferon_volcano_data.xlsx"
df_interferon.to_excel(output_path, index=False)

plt.figure(figsize=(10, 6))
sns.scatterplot(data=df_interferon, x="Log2FC", y=-np.log10(df_interferon["p_value"]),
                hue="Significant", palette={True: "red", False: "gray"})
plt.axhline(-np.log10(0.05), linestyle="--", color="black")
plt.axvline(1, linestyle="--", color="black")
plt.axvline(x=-1, color='black', linestyle='--')
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 p-value")
plt.title("Volcano Plot of Interferon Gamma Production Geneset")
plt.xlim(-1.5, 1.5)  # show only x from -5 to 10
plt.ylim(-1.5, 20)  # show only y from -5 to 10

plt.show()

print(f"Interferon volcano plot data saved to: {output_path}")
