from statsmodels.stats.multitest import multipletests
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import ranksums
import gseapy as gp

# ---------------------------
# 1. Load data
# ---------------------------
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
adata = sc.read_h5ad(file_path)

# Subset GBM and PBMC NK cells
adata_gbm = adata[adata.obs['source'] == 'glioblastoma']
adata_pbmc = adata[adata.obs['source'] == 'PBMC']

# TGFB-related gene list
gene_list = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]

## TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]
#TGFSIG = [
#     "JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1",
#     "CTNNB1", "IFNG", "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300",
#     "SMURF2", "TFRC", "SPTBN1", "ATF3", "MAPK14", "MAPK13", "APC",
#     "MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111", "MAPK9",
#     "ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12"
# ]
# ---------------------------
# 2. GBM vs PBMC DE for TGFB-related genes
# ---------------------------
results = []
for gene in gene_list:
    if gene in adata.var_names:
        expr_gbm = adata_gbm[:, gene].X.toarray().flatten()
        expr_pbmc = adata_pbmc[:, gene].X.toarray().flatten()

        stat, p_value = ranksums(expr_gbm, expr_pbmc)
        log2fc = np.log2(np.mean(expr_gbm) + 1) - np.log2(np.mean(expr_pbmc) + 1)

        results.append({'Gene': gene, 'Log2FC': log2fc, 'p_value': p_value})

df_results = pd.DataFrame(results)

# Adjust p-values (FDR)
reject, adj_pvals, _, _ = multipletests(df_results['p_value'], method='fdr_bh')
df_results['adj_p_value'] = adj_pvals
df_results['Significant'] = (df_results['adj_p_value'] < 0.05) & (df_results['Log2FC'] > 1)

significant_genes = df_results[df_results['Significant']]['Gene'].tolist()

print("Significant TGFB-related genes:", significant_genes)

# ---------------------------
# 3. Label GBM cells as TGFB-high or TGFB-low
# ---------------------------
# expr_df = adata_gbm[:, significant_genes].to_df() > 0
# gbm_high_expr_cells = expr_df.sum(axis=1) >= len(significant_genes) / 2
# gbm_high_expr_cells = (adata_gbm[:, significant_genes].to_df() > 0).all(axis=1)
gbm_high_expr_cells = adata_gbm[:, significant_genes].to_df().sum(axis=1) > 0
adata_gbm.obs['TGFB_status'] = "Low"
adata_gbm.obs.loc[gbm_high_expr_cells, 'TGFB_status'] = "High"

# ---------------------------
# 4. Differential expression between TGFB-high vs TGFB-low GBM cells
# ---------------------------
sc.tl.rank_genes_groups(adata_gbm, groupby="TGFB_status", reference="Low", method="wilcoxon")

# Save DE results
de_df = sc.get.rank_genes_groups_df(adata_gbm, group="High")
de_df = de_df[['names', 'logfoldchanges', 'pvals_adj']]
de_df.to_csv("DELETE_DE_TGFB_High_vs_Low_GBM.csv", index=False)

# ---------------------------
# 5. Prepare ranked list for GSEA
# ---------------------------
rnk = de_df[['names', 'logfoldchanges']].sort_values('logfoldchanges', ascending=False)
rnk_path = "DELETE_high_vs_low_gbm.rnk"
rnk.to_csv(rnk_path, sep="\t", index=False, header=False)

# ---------------------------
# 6. Run GSEA
# ---------------------------
pre_res = gp.prerank(
    rnk=rnk_path,
    gene_sets="KEGG_2016",  # Can change to GO_Biological_Process_2021, Reactome_2016, etc.
    outdir="DELETEgsea_results_tgfb_high_vs_low",
    min_size=15,
    max_size=500,
    permutation_num=1000,
    seed=42
)
adata_gbm.obs['TGFB_status'] = "Low"
adata_gbm.obs.loc[gbm_high_expr_cells, 'TGFB_status'] = "High"
high_count = (adata_gbm.obs['TGFB_status'] == "High").sum()
print(f"Number of TGFB-high GBM cells: {high_count}")
# ---------------------------
# 7. Save GSEA summary table
# ---------------------------
gsea_df = pre_res.res2d
gsea_df.to_csv("DELETE_GSEA_TGFB_High_vs_Low_GBM.csv")
