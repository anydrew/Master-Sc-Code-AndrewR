# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:17:27 2025

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
from adjustText import adjust_text  # Correct import

file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)
adata_gbm = adata[adata.obs['source'] == 'glioblastoma']
adata_pbmc = adata[adata.obs['source'] == 'PBMC']
gene_list = ["CSNK2A1", "ELK1", "FOS", "GRB2", "HRAS", "IL2", "IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "JUN", "LCK", "MAP2K1", "MAPK3", "MAPK8", "RAF1", "SHC1", "SOS1", "STAT5A", "STAT5B", "SYK"] 

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
gbm_high_expr_cells = (adata_gbm[:, significant_genes].to_df() > 0).all(axis=1)
adata_gbm_high = adata_gbm[gbm_high_expr_cells]
GLYCOL_GLUCONEOGEN_genes = ['ADORA2B', 'AGL', 'AGRN', 'AK3', 'AK4', 'AKR1A1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALG1', 'ANG', 'ANGPTL4', 'ANKZF1', 'ARPP19', 'AURKA', 'B3GALT6', 'B3GAT1', 'B3GAT3', 'B4GALT1', 'B4GALT2', 'B4GALT4', 'B4GALT7', 'BIK', 'BPNT1', 'CACNA1H', 'CAPN5', 'CASP6', 'CD44', 'CDK1', 'CENPA', 'CHPF', 'CHPF2', 'CHST1', 'CHST12', 'CHST2', 'CITED2', 'CLDN3', 'CLN6', 'COG2', 'COL5A1', 'COPB2', 'CTH', 'CXCR4', 'CYB5A', 'DCN', 'DDIT4', 'DEPDC1', 'DLD', 'DSC2', 'ECD', 'EGFR', 'EGLN3', 'ELF3', 'ENO1', 'ENO2', 'ERO1A', 'EXT1', 'EXT2', 'FAM162A', 'FKBP4', 'FUT8', 'G6PD', 'GALE', 'GALK1', 'GALK2', 'GCLC', 'GFPT1', 'GLCE', 'GLRX', 'GMPPA', 'GMPPB', 'GNE', 'GNPDA1', 'GOT1', 'GOT2', 'GPC1', 'GPC3', 'GPC4', 'GUSB', 'GYS1', 'HAX1', 'HDLBP', 'HK2', 'HMMR', 'HOMER1', 'HS2ST1', 'HS6ST2', 'HSPA5', 'IDH1', 'IDUA', 'IER3', 'IGFBP3', 'IL13RA1', 'IRS2', 'ISG20', 'KDELR3', 'KIF20A', 'KIF2A', 'LDHA', 'LHPP', 'MDH1', 'MDH2', 'ME2', 'MED24', 'MET', 'MIF', 'MPI', 'MXI1', 'NANP', 'NASP', 'NDUFV3', 'NOL3', 'NSDHL', 'NT5E', 'P4HA1', 'P4HA2', 'PAM', 'PAXIP1', 'PC', 'PDK3', 'PFKM', 'PFKL', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGLS', 'PGM2', 'PHKA2', 'PKM', 'PKP2', 'PLOD1', 'PLOD2', 'PMM2', 'POLR3K', 'PPIA', 'PPP2CB', 'PRPS1', 'PSMC4', 'PYGB', 'PYGL', 'QSOX1', 'RBCK1', 'RPE', 'RRAGD', 'SAP30', 'SDC1', 'SDC2', 'SDC3', 'SDHC', 'SLC16A3', 'SLC25A13', 'SLC35A3', 'SLC37A4', 'SOD1', 'SPAG4', 'SRD5A3', 'STMN1', 'TALDO1', 'TFF3', 'TGFA', 'TGFBI', 'TPBG', 'TPI1', 'TPST1', 'TXN', 'UGP2', 'VCAN', 'VEGFA', 'VLDLR', 'XYLT2', 'ZNF292']
interferon_results = []
for gene in GLYCOL_GLUCONEOGEN_genes:
    if gene in adata.var_names:
        expr_high = adata_gbm_high[:, gene].X.toarray().flatten()
        expr_rest = adata_gbm[:, gene].X.toarray().flatten()

        stat, p_value = ranksums(expr_high, expr_rest)
        log2fc = np.log2(np.mean(expr_high) + 1) - np.log2(np.mean(expr_rest) + 1)

        interferon_results.append({'Gene': gene, 'Log2FC': log2fc, 'p_value': p_value})

df_interferon = pd.DataFrame(interferon_results)
df_interferon['Significant'] = (df_interferon['p_value'] < 0.05) & (df_interferon['Log2FC'] > 0.58)
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df_interferon, x="Log2FC", y=-np.log10(df_interferon["p_value"]),
                hue="Significant", palette={True: "red", False: "gray"})
plt.xlim(-3, 3)  # Adjust based on actual min/max Log2FC values
plt.axhline(-np.log10(0.05), linestyle="--", color="black")
plt.axvline(0.58, linestyle="--", color="black")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 p-value")
plt.title("Volcano Plot of OTHER NADPH DEPENDENT Family sig  Geneset")
plt.show()
# # Save the interferon results to an Excel file in the same folder
# output_path = file_path.rsplit("\\", 1)[0] + "\\NADPH DEPENDENT ENZYMES Family _expression_results_inhighIL2cells.xlsx"
# df_interferon.to_excel(output_path, index=False)

# print(f"Exported results to: {output_path}")

# significant_genes = df_interferon[(df_interferon['p_value'] < 0.05) & (abs(df_interferon['Log2FC']) > 1)]
# # plt.scatter(significant_genes['Log2FC'], significant_genes['p_value'], color='red')

# # # Add labels to significant genes
# # texts = []
# # for i, row in significant_genes.iterrows():
# #     texts.append(plt.text(row['Log2 Fold Change'], row['-Log10 p-value'], row['gene'], fontsize=8))

# # # Adjust text to prevent overlap
# # adjust_text(texts)


# output_file_path = r"C:\Users\aroche6\abc\s123.xlsx"
# significant_genes.to_excel(output_file_path, index=False)

# Define output file path for Excel export
# output_file = r"C:\Users\aroche6\abc\890.xlsx"
# combined_data = pd.concat(
#     [mean_expression_df, fraction_of_cells_df], 
#     keys=["Mean Expression", "Fraction of Cells"], 
#     axis=1
# )

# # Save data to Excel
# combined_data.to_excel(output_file)

# print(f"Dot plot data with fractions of cells extracted and saved to {output_file}")
# # Export the results to an Excel file with two sheets
# with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
#     df_results.to_excel(writer, sheet_name='Gene_Comparison', index=False)
#     df_interferon.to_excel(writer, sheet_name='Interferon_Gamma_Geneset', index=False)