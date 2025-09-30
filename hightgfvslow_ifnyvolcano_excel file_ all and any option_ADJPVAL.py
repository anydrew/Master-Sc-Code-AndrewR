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
from statsmodels.stats.multitest import multipletests

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
gene_list = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]
# TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]
# TGFBSIGgene_list = [
#     "JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1",
#     "CTNNB1", "IFNG", "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300",
#     "SMURF2", "TFRC", "SPTBN1", "ATF3", "MAPK14", "MAPK13", "APC",
#     "MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111", "MAPK9",
#     "ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12"
# ]
# TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]
#MITOcomp1 = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'SDHA', 'SDHB', 'SDHC', 'UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']
# gene_list = [
#     'ABCC1', 'ABCD1', 'ABL1', 'ADAM9', 'AGAP3', 'AIF1', 'AIFM1', 'AIFM2', 'AKT1', 'ALAD', 
#     'ALDH3B1', 'ALOX5', 'ALS2', 'ANKZF1', 'ANXA1', 'APOD', 'APOE', 'APP', 'AREG', 'ARL6IP5', 
#     'ARNT', 'ATF2', 'ATF4', 'ATM', 'ATOX1', 'ATP13A2', 'ATP2A2', 'ATP7A', 'AXL', 'BAK1', 
#     'BANF1', 'BCL2', 'BECN1', 'BMP7', 'BNIP3', 'BRF2', 'BTK', 'CAMKK2', 'CAPN2', 'CASP3', 
#     'CAT', 'CBX8', 'CCS', 'CD36', 'CD38', 'CDK1', 'CHCHD2', 'CHUK', 'COL1A1', 'CPEB2', 
#     'CRK', 'CRYAB', 'CYB5B', 'CYGB', 'CYP1B1', 'DAPK1', 'DDR2', 'DHCR24', 'DHFR', 'DUOX1', 
#     'ECT2', 'EDN1', 'EGFR', 'EIF2S1', 'EPAS1', 'ERCC1', 'ERCC2', 'ERCC3', 'ERCC6', 
#     'ERCC6L2', 'ERCC8', 'ERMP1', 'ERN1', 'ERO1A', 'ETFDH', 'ETV5', 'EZH2', 'FANCC', 
#     'FANCD2', 'FBLN5', 'FER', 'FKBP1B', 'FOS', 'FOSL1', 'FOXO1', 'FOXO3', 'FOXO4', 
#     'FOXP1', 'FUT8', 'FXN', 'FYN', 'G6PD', 'GCH1', 'GCLC', 'GCLM', 'GGT7', 'GJB2', 
#     'GLRX2', 'GPX2', 'GPX3', 'GPX4', 'GPX7', 'GPX8', 'GSR', 'GSS', 'GSTP1', 'HDAC2', 
#     'HDAC6', 'HIF1A', 'HM13', 'HMOX1', 'HMOX2', 'HSF1', 'HSPA1A', 'HSPA1B', 'HTRA2', 
#     'HYAL1', 'HYAL2', 'IDH1', 'IL18BP', 'IL18RAP', 'IL1A', 'IL6', 'IPCEF1', 'JAK2', 
#     'JUN', 'KAT2B', 'KDM6B', 'KEAP1', 'KLF2', 'LIAS', 'LONP1', 'LRRK2', 'MAP1LC3A', 
#     'MAP2K4', 'MAP3K5', 'MAPK1', 'MAPK13', 'MAPK3', 'MAPK7', 'MAPK8', 'MAPK9', 
#     'MAPKAP1', 'MAPT', 'MCTP1', 'MDM2', 'MET', 'MGAT3', 'MGST1', 'MICB', 'MITF', 
#     'MLH1', 'MMP7', 'MPO', 'MSRA', 'MSRB1', 'MTHFD1', 'MTHFD2', 'MTHFR', 'NADK', 
#     'NDUFA9', 'NDUFS1', 'NDUFS2', 'NFE2L1', 'NFE2L2', 'NLRP3', 'NOX4', 'NQO1', 
#     'NRG1', 'NRG2', 'OGG1', 'OLR1', 'OTUD1', 'OXSR1', 'P4HB', 'PA2G4', 'PARK7', 
#     'PDXK', 'PGD', 'PGK1', 'PGM3', 'PIK3R1', 'PIK3R5', 'PKM', 'PLA2G4A', 'PLA2G5', 
#     'PLCB1', 'PLCB3', 'PLCG1', 'PLCG2', 'PLK1', 'PML', 'PNP', 'POLB', 'POLG', 
#     'PPARGC1A', 'PRDX1', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6', 'PRKAA1', 
#     'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKCQ', 
#     'PXDN', 'RAB20', 'RAD23A', 'RAP1B', 'RB1', 'RBM38', 'RCAN1', 'RNF4', 'ROS1', 
#     'RUNX2', 'SIRT1', 'SIRT3', 'SIRT5', 'SLC11A2', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
#     'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A1', 
#     'SLC25A12', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
#     'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC25A4', 'SLC25A5', 'SLC25A6', 
#       'SLC25A12', 'SLC25A15', 'SLC25A19', 
#     'SLC25A22', 'SLC25A25', 'SLC25A29', 'SLC25A32', 'SLC25A38', 'SLC25A39', 
#     'SOD1', 'SOD2', 'SOD3', 'SRXN1', 'STK11', 'TBXAS1', 'TEK', 'TFAM', 'TGFB1', 
#     'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TGFBR3', 'TLR2', 'TLR4', 
#     'TNF', 'TNFAIP3', 'TNFRSF1A', 'TNFRSF1B', 'TP53', 'TXN', 'TXNDC5', 'TXNIP', 
#     'UQCRC1', 'UQCRC2', 'UQCRFS1', 'VHL', 'XPC', 'ZFP36'
# ]

results = []
for gene in gene_list:
    if gene in adata.var_names:
        expr_gbm = adata_gbm[:, gene].X.toarray().flatten()
        expr_pbmc = adata_pbmc[:, gene].X.toarray().flatten()

        stat, p_value = ranksums(expr_gbm, expr_pbmc)
        log2fc = np.log2(np.mean(expr_gbm) + 1) - np.log2(np.mean(expr_pbmc) + 1)

        results.append({'Gene': gene, 'Log2FC': log2fc, 'p_value': p_value})

df_results = pd.DataFrame(results)
# Adjust p-values for multiple testing
reject, adj_pvals, _, _ = multipletests(df_results['p_value'], method='fdr_bh')
df_results['adj_p_value'] = adj_pvals
df_results['Significant'] = (df_results['adj_p_value'] < 0.05) & (df_results['Log2FC'] > 1)
significant_genes = df_results[df_results['Significant']]['Gene'].tolist()
#any
gbm_high_expr_cells = adata_gbm[:, significant_genes].to_df().sum(axis=1) > 0
#all
#gbm_high_expr_cells = (adata_gbm[:, significant_genes].to_df() > 0).all(axis=1)

# gbm_low_expr_cells = adata_gbm[:, significant_genes].to_df().sum(axis=1) < 0
#all
#gbm_low_expr_cells = (adata_gbm[:, significant_genes].to_df() < 0).all(axis=1)

adata_gbm_high = adata_gbm[gbm_high_expr_cells]
# Low = expresses NONE of the significant genes (all 0)
gbm_low_expr_cells = (adata_gbm[:, significant_genes].to_df() == 0).all(axis=1)
# gbm_low_expr_cells = ~gbm_high_expr_cells
adata_gbm_low = adata_gbm[gbm_low_expr_cells]
interferon_gamma_genes = ["ACVR1", "BID", "CD244", "CD247", "CD48", "FCER1G", "FCGR3A", "FCGR3B", "GZMB", "HCST", "HLA-A", "HLA-C", "HLA-E", "HLA-G", "ICAM1", "ICAM2", "IFNAR1", "IFNAR2", "IFNG", "IFNGR1", "IFNGR2", "ITGAL", "ITGB2", "KIR2DL3", "KIR2DL4", "KLRC1", "KLRC2", "KLRD1", "LAT", "LCK", "LCP2", "MICB", "NCR1", "NCR3", "NFATC1", "NFATC2", "PRF1", "PTK2B", "RAC2", "RAC3", "SH2D1A", "SH2D1B", "SH3BP2", "SHC1", "SYK", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D", "TNFSF10", "TYROBP", "ULBP1", "ULBP2", "VAV1", "VAV3", "ZAP70", 'GZMB', 'GZMH', 'GZMK', 'GZMM']

# interferon_gamma_genes = [
#     "ABL1", "AXL", "BCL3", "C1QBP", "CD2", "CD3E", "CD14", "CD47", "CEBPG", 
#     "CCR7", "CR1", "DDIT3", "F2RL1", "GAS6", "GATA3", "HLA-A", "HLA-DPA1",
#     "HLA-DPB1", "HLA-DRB1", "HMGB1", "HRAS", "HSPD1", "IRF8", "IL1B",
#     "IL1R1", "IL2", "IL10", "IL12A", "IL12B", "IL12RB1", "IL12RB2",
#     "IL18", "INHA", "INHBA", "ISL1", "JAK2", "LGALS9", "LTA", "PDE4B",
#     "PDE4D", "PRNP", "RARA", "TRIM27", "XCL1", "SLAMF1", "SLC11A1",
#     "TLR3", "TLR4", "TNF", "TNFSF4", "TXK", "TYK2", "SCGB1A1", "WNT5A",
#     "ZP3", "LAPTM5", "FZD5", "SLC7A5", "RIPK2", "FADD", "IL18R1",
#     "PGLYRP1", "IL1RL1", "IL27RA", "ISG15", "NR1H4", "RASGRP1", "EBI3",
#     "CD96", "CD226", "LILRB1", "ARID5A", "LILRB4", "RIPK3", "BTN3A2",
#     "BTN3A1", "CD160", "KLRK1", "SCRIB", "PTPN22", "IL36RN", "PYCARD",
#     "CD274", "FOXP3", "TLR7", "TLR8", "IL23A", "CYRIB", "CD244", "IL20RB",
#     "TLR9", "SASH3", "CRTAM", "HMHB1", "IL21", "VSIR", "NOD2", "CLEC7A",
#     "ZC3H12A", "PDCD1LG2", "CD276", "HAVCR2", "PGLYRP2", "SLAMF6",
#     "SIRPA", "IL23R", "ZFPM1", "NLRP6", "IL27", "IFNL1", "LGALS9B",
#     "MIR24-1", "LGALS7B", "LGALS9C", "CCR2", "MIR708", "KLRC4-KLRK1"
# ]
# NK_CYTO_GENES = ["ACVR1", "BID", "CD244", "CD247", "CD48", "FCER1G", "FCGR3A", "FCGR3B", "GZMB", "HCST", "HLA-A", "HLA-C", "HLA-E", "HLA-G", "ICAM1", "ICAM2", "IFNAR1", "IFNAR2", "IFNG", "IFNGR1", "IFNGR2", "ITGAL", "ITGB2", "KIR2DL3", "KIR2DL4", "KLRC1", "KLRC2", "KLRD1", "LAT", "LCK", "LCP2", "MICB", "NCR1", "NCR3", "NFATC1", "NFATC2", "PRF1", "PTK2B", "RAC2", "RAC3", "SH2D1A", "SH2D1B", "SH3BP2", "SHC1", "SYK", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D", "TNFSF10", "TYROBP", "ULBP1", "ULBP2", "VAV1", "VAV3", "ZAP70", 'GZMB', 'GZMH', 'GZMK', 'GZMM']

# OXPHOS_genes = ['ATP1B1', 'COX10', 'UQCRH', 'ATP6V1B2', 'NDUFA5', 'NDUFA4', 'ATP6V1C1', 'NDUFC2', 'CLDN3', 'COX6C', 'LHPP', 'UQCRFS1', 'NDUFB2', 'NDUFA10', 'NDUFA4L2', 'NDUFS5', 'CD4', 'NDUFA13', 'RERE', 'UQCRQ', 'UQCRB', 'COX17', 'NDUFB7', 'COX11', 'NDUFS1', 'NDUFS4', 'ATP6V0B', 'COX7A2L', 'ATP6V0C', 'TPP1', 'NDUFB3', 'ATP6V1D', 'SURF1', 'CYC1', 'NDUFB5', 'NDUFV2', 'COX6A1', 'COX4I1', 'COX8A', 'NDUFB11', 'NDUFA12', 'COX7A2', 'NDUFA2', 'NDUFC1', 'PPA1', 'NDUFA9', 'NDUFB6', 'NDUFV1', 'NDUFB1', 'ATP6V0E2', 'COX7C', 'ATP6V1E1', 'COX5A', 'ATP6V1G1', 'COX5B', 'NDUFS6', 'ATP6V0E1', 'CPM', 'COX7B', 'ATP6V1A', 'CP', 'NDUFAB1', 'UQCRC2', 'NDUFA11', 'NDUFS8', 'NDUFS3', 'NDUFA6', 'ATP6V1E2', 'UQCR10', 'NDUFA3', 'COX15', 'NDUFB9', 'SDHA', 'NDUFA8', 'ATP6V0A2', 'SDHC', 'CPE', 'ATP6V1H', 'NDUFV3', 'ATP6V0D1', 'NDUFS7', 'NDUFA1', 'NDUFB10', 'TCIRG1', 'NDUFS2', 'UQCRC1', 'ATP6V1F', 'PPA2', 'NDUFB4', 'ATP6AP1', 'ATP6V0A1', 'SDHB', 'COX6B1']
#GLYCOL_GLUCONEOGEN_genes = ['ADORA2B', 'AGL', 'AGRN', 'AK3', 'AK4', 'AKR1A1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALG1', 'ANG', 'ANGPTL4', 'ANKZF1', 'ARPP19', 'AURKA', 'B3GALT6', 'B3GAT1', 'B3GAT3', 'B4GALT1', 'B4GALT2', 'B4GALT4', 'B4GALT7', 'BIK', 'BPNT1', 'CACNA1H', 'CAPN5', 'CASP6', 'CD44', 'CDK1', 'CENPA', 'CHPF', 'CHPF2', 'CHST1', 'CHST12', 'CHST2', 'CITED2', 'CLDN3', 'CLN6', 'COG2', 'COL5A1', 'COPB2', 'CTH', 'CXCR4', 'CYB5A', 'DCN', 'DDIT4', 'DEPDC1', 'DLD', 'DSC2', 'ECD', 'EGFR', 'EGLN3', 'ELF3', 'ENO1', 'ENO2', 'ERO1A', 'EXT1', 'EXT2', 'FAM162A', 'FKBP4', 'FUT8', 'G6PD', 'GALE', 'GALK1', 'GALK2', 'GCLC', 'GFPT1', 'GLCE', 'GLRX', 'GMPPA', 'GMPPB', 'GNE', 'GNPDA1', 'GOT1', 'GOT2', 'GPC1', 'GPC3', 'GPC4', 'GUSB', 'GYS1', 'HAX1', 'HDLBP', 'HK2', 'HMMR', 'HOMER1', 'HS2ST1', 'HS6ST2', 'HSPA5', 'IDH1', 'IDUA', 'IER3', 'IGFBP3', 'IL13RA1', 'IRS2', 'ISG20', 'KDELR3', 'KIF20A', 'KIF2A', 'LDHA', 'LHPP', 'MDH1', 'MDH2', 'ME2', 'MED24', 'MET', 'MIF', 'MPI', 'MXI1', 'NANP', 'NASP', 'NDUFV3', 'NOL3', 'NSDHL', 'NT5E', 'P4HA1', 'P4HA2', 'PAM', 'PAXIP1', 'PC', 'PDK3', 'PFKM', 'PFKL', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGLS', 'PGM2', 'PHKA2', 'PKM', 'PKP2', 'PLOD1', 'PLOD2', 'PMM2', 'POLR3K', 'PPIA', 'PPP2CB', 'PRPS1', 'PSMC4', 'PYGB', 'PYGL', 'QSOX1', 'RBCK1', 'RPE', 'RRAGD', 'SAP30', 'SDC1', 'SDC2', 'SDC3', 'SDHC', 'SLC16A3', 'SLC25A13', 'SLC35A3', 'SLC37A4', 'SOD1', 'SPAG4', 'SRD5A3', 'STMN1', 'TALDO1', 'TFF3', 'TGFA', 'TGFBI', 'TPBG', 'TPI1', 'TPST1', 'TXN', 'UGP2', 'VCAN', 'VEGFA', 'VLDLR', 'XYLT2', 'ZNF292']


# ros_genes = [
#     "ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC", 
#     "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB", 
#     "LAMTOR5", "MBP", "MGST1", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2", "NQO1", "OXSR1", 
#     "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "SBNO2", "SCAF4", 
#     "SOD1", "SOD2", "STK25", "TXN", "TXNRD1", "TXNRD2"
# ]

#MITOcomp1 = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'SDHA', 'SDHB', 'SDHC', 'UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']


# oxidative_stress_genes = [
#     'ABCC1', 'ABCD1', 'ABL1', 'ADAM9', 'AGAP3', 'AIF1', 'AIFM1', 'AIFM2', 'AKT1', 'ALAD', 
#     'ALDH3B1', 'ALOX5', 'ALS2', 'ANKZF1', 'ANXA1', 'APOD', 'APOE', 'APP', 'AREG', 'ARL6IP5', 
#     'ARNT', 'ATF2', 'ATF4', 'ATM', 'ATOX1', 'ATP13A2', 'ATP2A2', 'ATP7A', 'AXL', 'BAK1', 
#     'BANF1', 'BCL2', 'BECN1', 'BMP7', 'BNIP3', 'BRF2', 'BTK', 'CAMKK2', 'CAPN2', 'CASP3', 
#     'CAT', 'CBX8', 'CCS', 'CD36', 'CD38', 'CDK1', 'CHCHD2', 'CHUK', 'COL1A1', 'CPEB2', 
#     'CRK', 'CRYAB', 'CYB5B', 'CYGB', 'CYP1B1', 'DAPK1', 'DDR2', 'DHCR24', 'DHFR', 'DUOX1', 
#     'ECT2', 'EDN1', 'EGFR', 'EIF2S1', 'EPAS1', 'ERCC1', 'ERCC2', 'ERCC3', 'ERCC6', 
#     'ERCC6L2', 'ERCC8', 'ERMP1', 'ERN1', 'ERO1A', 'ETFDH', 'ETV5', 'EZH2', 'FANCC', 
#     'FANCD2', 'FBLN5', 'FER', 'FKBP1B', 'FOS', 'FOSL1', 'FOXO1', 'FOXO3', 'FOXO4', 
#     'FOXP1', 'FUT8', 'FXN', 'FYN', 'G6PD', 'GCH1', 'GCLC', 'GCLM', 'GGT7', 'GJB2', 
#     'GLRX2', 'GPX2', 'GPX3', 'GPX4', 'GPX7', 'GPX8', 'GSR', 'GSS', 'GSTP1', 'HDAC2', 
#     'HDAC6', 'HIF1A', 'HM13', 'HMOX1', 'HMOX2', 'HSF1', 'HSPA1A', 'HSPA1B', 'HTRA2', 
#     'HYAL1', 'HYAL2', 'IDH1', 'IL18BP', 'IL18RAP', 'IL1A', 'IL6', 'IPCEF1', 'JAK2', 
#     'JUN', 'KAT2B', 'KDM6B', 'KEAP1', 'KLF2', 'LIAS', 'LONP1', 'LRRK2', 'MAP1LC3A', 
#     'MAP2K4', 'MAP3K5', 'MAPK1', 'MAPK13', 'MAPK3', 'MAPK7', 'MAPK8', 'MAPK9', 
#     'MAPKAP1', 'MAPT', 'MCTP1', 'MDM2', 'MET', 'MGAT3', 'MGST1', 'MICB', 'MITF', 
#     'MLH1', 'MMP7', 'MPO', 'MSRA', 'MSRB1', 'MTHFD1', 'MTHFD2', 'MTHFR', 'NADK', 
#     'NDUFA9', 'NDUFS1', 'NDUFS2', 'NFE2L1', 'NFE2L2', 'NLRP3', 'NOX4', 'NQO1', 
#     'NRG1', 'NRG2', 'OGG1', 'OLR1', 'OTUD1', 'OXSR1', 'P4HB', 'PA2G4', 'PARK7', 
#     'PDXK', 'PGD', 'PGK1', 'PGM3', 'PIK3R1', 'PIK3R5', 'PKM', 'PLA2G4A', 'PLA2G5', 
#     'PLCB1', 'PLCB3', 'PLCG1', 'PLCG2', 'PLK1', 'PML', 'PNP', 'POLB', 'POLG', 
#     'PPARGC1A', 'PRDX1', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6', 'PRKAA1', 
#     'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKCQ', 
#     'PXDN', 'RAB20', 'RAD23A', 'RAP1B', 'RB1', 'RBM38', 'RCAN1', 'RNF4', 'ROS1', 
#     'RUNX2', 'SIRT1', 'SIRT3', 'SIRT5', 'SLC11A2', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
#     'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A1', 
#     'SLC25A12', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
#     'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC25A4', 'SLC25A5', 'SLC25A6', 
#       'SLC25A12', 'SLC25A15', 'SLC25A19', 
#     'SLC25A22', 'SLC25A25', 'SLC25A29', 'SLC25A32', 'SLC25A38', 'SLC25A39', 
#     'SOD1', 'SOD2', 'SOD3', 'SRXN1', 'STK11', 'TBXAS1', 'TEK', 'TFAM', 'TGFB1', 
#     'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TGFBR3', 'TLR2', 'TLR4', 
#     'TNF', 'TNFAIP3', 'TNFRSF1A', 'TNFRSF1B', 'TP53', 'TXN', 'TXNDC5', 'TXNIP', 
#     'UQCRC1', 'UQCRC2', 'UQCRFS1', 'VHL', 'XPC', 'ZFP36'
# ]
# mitofission_genes = ['VPS35', 'SLC25A46', 'MTFR1', 'COX10', 'TMEM135', 'DCN', 'AURKA', 'MTFR2', 'MUL1', 'MTFP1', 'AP3B1', 'DDHD1', 'MIEF1', 'MCU', 'FIS1', 'MFF', 'PINK1', 'MYO19', 'DDHD2', 'OPA1', 'RALBP1', 'RALA', 'INF2', 'BAX', 'CCAR2', 'UCP2', 'LRRK2', 'STAT2', 'SPIRE1', 'ATG3', 'PGAM5', 'LPIN1', 'MAPT', 'GDAP1', 'MIEF2', 'BNIP3', 'MTFR1L', 'DNM1L', 'PPARG']
# mitofussion_genes = ['MFN1', 'AFG3L2', 'MFN2', 'ADCK1', 'PID1', 'MUL1', 'STOML2', 'MIEF1', 'THG1L', 'ZDHHC6', 'FIS1', 'MFF', 'MTCH2', 'HUWE1', 'OPA1', 'BAK1', 'USP30', 'BAX', 'MCL1', 'GDAP1', 'VAT1', 'MIEF2', 'BNIP3', 'PLD6', 'CHCHD3', 'OMA1', 'PARL', 'TFRC', 'BCL2A1']

# cytochrome_c_release_apop_genes = ['AVEN', 'BBC3', 'BCL2L1', 'BMF', 'CASP3', 'PMAIP1', 'XIAP', 'PINK1', 'CARD8', 'MLLT11', 'NOL3', 'CASP9', 'BIRC2', 'AIFM1', 'CYCS', 'APIP', 'DNM1L', 'PARL', 'BCL2A1', 'BIRC3', 'UACA', 'ENDOG', 'BIK', 'DFFA', 'DIABLO', 'CIDEB', 'SOD2', 'FIS1', 'MFF', 'BAK1', 'BAX', 'GGCT', 'MMP9', 'FAM162A', 'BNIP3', 'TNFSF10', 'MAPK1', 'PPIF', 'APAF1', 'MOAP1', 'DFFB', 'TP53', 'BCL2', 'CASP8', 'IGF1', 'PDCD5', 'LMNA', 'CASP6', 'FXN', 'PSMD10', 'PYCARD', 'MAPK3', 'TIMM50', 'GHITM', 'BAD', 'BOK', 'PLAUR', 'PRELID1', 'SFN', 'OPA1', 'BCL2L2', 'MCL1', 'BID', 'CLU', 'AKT1', 'IFI6', 'CASP7', 'BCL2L11', 'TRIAP1']
# mitophagy_genes = ['ATG12', 'ATG5', 'CSNK2A1', 'CSNK2A2', 'CSNK2B', 'FUNDC1', 'MAP1LC3A', 'MAP1LC3B', 'MFN1', 'MFN2', 'MTERF3', 'PGAM5', 'PINK1', 'SQSTM1', 'SRC', 'TOMM20', 'TOMM22', 'TOMM40', 'TOMM5', 'TOMM7', 'UBB', 'UBC', 'ULK1', 'VDAC1']

interferon_results = []
for gene in interferon_gamma_genes:
    if gene in adata.var_names:
        expr_high = adata_gbm_high[:, gene].X.toarray().flatten()
        expr_low = adata_gbm_low[:, gene].X.toarray().flatten()

        stat, p_value = ranksums(expr_high, expr_low)
        log2fc = np.log2(np.mean(expr_high) + 1) - np.log2(np.mean(expr_low) + 1)

        interferon_results.append({'Gene': gene, 'Log2FC': log2fc, 'p_value': p_value})

df_interferon = pd.DataFrame(interferon_results)
# Adjust p-values (FDR-BH)
reject, adj_pvals, _, _ = multipletests(df_interferon['p_value'], method='fdr_bh')
df_interferon['adj_p_value'] = adj_pvals


df_interferon['Significant'] = (df_interferon['adj_p_value'] < 0.05) & (df_interferon['Log2FC'] > 1)
significant_genes = df_results[df_results['Significant']]['Gene'].tolist()
#output_path = r"C:\Users\aroche6\abc\all_sig_tgfb_interferon_volcano_data.xlsx"
#output_path = r"C:\Users\aroche6\abc\any_sig_tgfb_interferon_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\all_ranked_sig_tgfb_interferon_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\any_ranked_sig_cytoxic_interferon_volcano_data.xlsx"

#output_path = r"C:\Users\aroche6\abc\all_sig_tgfb_cytox_volcano_data.xlsx"
#output_path = r"C:\Users\aroche6\abc\any_sig_tgfb_cytox_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\all_ranked_sig_tgfb_cytox_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\any_ranked_sig_tgfb_cytox_volcano_data.xlsx"

#output_path = r"C:\Users\aroche6\abc\all_sig_tgfb_oxphos_volcano_data.xlsx"
#output_path = r"C:\Users\aroche6\abc\any_sig_tgfb_oxphos_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\all_ranked_sig_tgfb_oxphos_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\any_ranked_sig_tgfb_oxphos_volcano_data.xlsx"

#output_path = r"C:\Users\aroche6\abc\all_sig_tgfb_gluco_volcano_data.xlsx"
#output_path = r"C:\Users\aroche6\abc\any_sig_tgfb_gluco_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\all_ranked_sig_tgfb_gluco_volcano_data.xlsx"
# output_path = r"C:\Users\aroche6\abc\any_ranked_sig_tgfb_gluco_volcano_data.xlsx"


# output_path = r"C:\Users\aroche6\abc\any_ranked_sig_tgfb_mitocomp_volcano_data.xlsx"
# df_interferon.to_excel(output_path, index=False)

plt.figure(figsize=(10, 6))
sns.scatterplot(data=df_interferon, x="Log2FC", y=-np.log10(df_interferon["adj_p_value"]),
                hue="Significant", palette={True: "red", False: "gray"})
plt.axhline(-np.log10(0.05), linestyle="--", color="black")
plt.axvline(1, linestyle="--", color="black")
plt.axvline(x=-1, color='black', linestyle='--')
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 adjusted p-value")
plt.title("any_ranked_sig_cytochrome_C_volcano")
# plt.xlim(-2, 2)  # show only x from -5 to 10
# plt.ylim(-5, 50)  # show only y from -5 to 10

plt.show()

# print(f"Interferon volcano plot data saved to: {output_path}")
