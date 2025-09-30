# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 19:38:35 2024

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
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
sc.read_h5ad(file_path)
adata = sc.read_h5ad(file_path)
'adata_all_nk_after_mapping.h5ad'
adata

myfilter1 = (adata.obs['subset_source'] == 'CD56dim, brain_normal') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
myfilter2 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56dim, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC')
myfilter3 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56dim, brain_normal') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
myfilter4 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
myfilter5 = ((adata.obs['subset_source'] == 'CD56dim, lung_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, melanoma') | (adata.obs['subset_source'] == 'CD56dim, pancreas_tumor') | (adata.obs['subset_source'] == 'CD56dim, prostate_tumor') | (adata.obs['subset_source'] == 'CD56dim, sarcoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, breast_tumor'))
myfilter6 = ((adata.obs['subset_source'] == 'CD56bright, lung_tumor') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56bright, melanoma') | (adata.obs['subset_source'] == 'CD56bright, pancreas_tumor') | (adata.obs['subset_source'] == 'CD56bright, prostate_tumor') | (adata.obs['subset_source'] == 'CD56bright, sarcoma_tumor') | (adata.obs['subset_source'] == 'CD56bright, breast_tumor'))
myfilter7 = ((adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') |(adata.obs['subset_source'] == 'CD56bright, lung_tumor') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56bright, melanoma') | (adata.obs['subset_source'] == 'CD56bright, pancreas_tumor') | (adata.obs['subset_source'] == 'CD56bright, prostate_tumor') | (adata.obs['subset_source'] == 'CD56bright, sarcoma_tumor') | (adata.obs['subset_source'] == 'CD56bright, breast_tumor') | (adata.obs['subset_source'] == 'CD56dim, lung_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, melanoma') | (adata.obs['subset_source'] == 'CD56dim, pancreas_tumor') | (adata.obs['subset_source'] == 'CD56dim, prostate_tumor') | (adata.obs['subset_source'] == 'CD56dim, sarcoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, breast_tumor'))
myfilter8 = ((adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56bright, melanoma') | (adata.obs['subset_source'] == 'CD56dim, melanoma') | (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC'))
myfilter9 = (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
myfilter10 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor')
myfilter11 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')

myfilter0 = (adata.obs['bright_dim_subset'] == 'CD56dim') | (adata.obs['bright_dim_subset'] == 'CD56dim')
myfilter12 = (adata.obs['source'] == 'PBMC') | (adata.obs['source'] == 'glioblastoma')
myfilter13 = (adata.obs['source'] == 'PBMC') | (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
ax = sc.pl.dotplot(adata[myfilter12], nk_exc_genes, groupby='source', title="NK Cell Diagnostic Test: CD3, CD8, CD14 Expression")




Here they are in one line with quotes:

`"JUNB", "TXN", "PRDX2", "PRDX1", "IPCEF1"`
"PRELID1", "MCL1", "MAPK1", "TNFSF10", "BCL2", "CASP3", "GGCT", "BCL2L11"



"COX5A", "HIGD2A"

"COX5A", "NDUFS8", "NDUFB9", "RERE", "PPA1", "NDUFA10", "TCIRG1"


"MCL1", "HUWE1", "CHCHD3"

"DDHD1", "AP3B1"

"PPIA", "DDIT4", "TXN", "STMN1", "NASP", "TALDO1", "B4GALT1", "POLR3K", "ME2", "IER3", "IRS2"









sigtgfb = ["JUND", "JUNB", "ID2", "PPP1R15A", "FKBP1A", "CD9", "JUN", "MAPK1", "CTNNB1", "IFNG", "TFDP1", "ID3", "SKI", "MYC", "CUL1", "EP300", "SMURF2", "TFRC", "SPTBN1", "ATF3", "MAPK14", "MAPK13", "APC", "MAPK8", "NBL1", "NCOR2", "MAPK3", "HRAS", "RNF111", "MAPK9", "ACVR1B", "PPP2CB", "VPS54", "GDF11", "SMURF1", "PPP2R1B", "MAPK12" ]

ax = sc.pl.dotplot(filtered_adata, sigtgfb, groupby='source', title="Significantly Increased TGF-β genes GBM-NK vs HD-PBMC-NK")
sigoxstress = ["AREG", "FOS", "HSPA1A", "NR4A2", "GSTP1", "PPIA", "HSPA1B", "TXNIP", "TXN", "BCL2", "PRDX1", "CD38", "NCF1", "JUN", "MAPK1", "TNFAIP3", "NFE2L2", "NDUFS8", "PRDX2", "SETX", "HM13", "ERCC1", "CASP3", "PRDX3", "NUDT1", "SESN1", "KDM6B", "STK24", "EPAS1", "IPCEF1", "RHOB", "GCH1", "EZH2", "ATP2A2", "PPP1R15B", "OXR1", "MSRA", "FOXO1", "MAPK13", "AKT1", "KAT2B", "NUDT2", "GSR", "ABCC1", "MAPK8", "DHFR", "PYROXD1", "MAP2K4", "SLC25A24", "FUT8", "ABL1", "VRK2", "MAPK3", "MAP1LC3A", "SIRT1", "MAPK9", "CDK1", "SMPD3", "FER", "GCLC", "APOE", "AGAP3", "JAK2", "TRAP1", "AIF1", "PPP2CB", "PAWR", "APP", "MSRB2", "STX2", "AIFM2", "ERCC6", "FANCD2", "CPEB2", "PEX10", "DUOX1"]

ax = sc.pl.dotplot(filtered_adata, sigoxstress, groupby='source', title="Significantly Increased Oxidative Stress genes GBM-NK vs HD-PBMC-NK")

sigros = ["JUNB", "TXN", "PRDX1", "PRDX2", "IPCEF1", "TXNRD1", "MSRA", "GSR", "ABCC1", "GCLC", "SBNO2"
]

ax = sc.pl.dotplot(filtered_adata, sigros, groupby='source', title="Significantly Increased ROS genes GBM-NK vs HD-PBMC-NK")

sigcytoc = ["MCL1", "BCL2", "PRELID1", "MAPK1", "TNFSF10", "CASP3", "BCL2A1", "BCL2L11", "BCL2L1", "GGCT", "OPA1", "AKT1", "LMNA", "AVEN", "MAPK3", "CASP6", "ENDOG"
]

ax = sc.pl.dotplot(filtered_adata, sigcytoc, groupby='source', title="Significantly Increased Cytochrome C Release Induced Apoptosis genes GBM-NK vs HD-PBMC-NK")

"APOE", "AREG", "HSPA1B", "CDK1", "BCL2", "AIF1", "HSPA1A", "EPAS1", "SMPD3", "NCF1", "DUOX"

sigoxphos = ["COX5A", "NDUFS8", "NDUFB9", "PPA1", "NDUFA10", "TCIRG1", "ATP1B1", "RERE", "ATP6V0C", "ATP6V0A2", "COX10"
]
sigoxphosratelim = ["COX5A", "HIGD2A", "COX10"
]
sigmitocomp = ["COX5A", "NDUFS8", "NDUFA10"
]
sigmitofusion = ["MCL1", "HUWE1", "BCL2A1", "CHCHD3", "TFRC", "OPA1", "MFN1"
]
sigmitofission = ["DDHD1", "AP3B1", "OPA1", "TMEM135", "MTFR1", "SPIRE1", "COX10", "MCU", "AURKA"
]
sigmitophagy = ["CSNK2A2", "MAP1LC3A", "MFN1"
]
sigglycgluc = ["PPIA", "DDIT4", "TXN", "STMN1", "NASP", "TALDO1", "B4GALT1", "POLR3K", "IER3", "ME2", "B4GALT4", "PGM2", "IRS2", "ECD", "GMPPB", "GFPT1", "FUT8", "GALE", "RPE", "EXT1", "CDK1", "PMM2", "GCLC", "B3GALT6", "SAP30", "PPP2CB", "HS2ST1", "CASP6", "HMMR", "SLC25A13", "AURKA", "CENPA"]

ax = sc.pl.dotplot(filtered_adata, sigglycgluc, groupby='source', title="Significantly Increased Glycolysis and Gluconeogenic genes GBM-NK vs HD-PBMC-NK")


sigglycratelim = ["SLC16A1" ]

ax = sc.pl.dotplot(adata[myfilter4], LIST6, groupby='subset_source', title="TEST6")

nk_exc_genes = ['CD3D', 'CD3E', 'CD3G', 'CD14', 'CD8A']
ax = sc.pl.dotplot(adata[myfilter3], nk_exc_genes, groupby='subset_source', title="NK Cell Diagnostic Test: CD3, CD8, CD14 Expression")
nk_inc_gene = ['NCAM1']
ax = sc.pl.dotplot(adata[myfilter3], nk_inc_gene, groupby='subset_source', title="NK Cell Diagnostic Test: CD56 Expression")
nk_br_inc_genes = ['CCR7', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5']
ax = sc.pl.dotplot(adata[myfilter3], nk_br_inc_genes, groupby='subset_source', title="NK Brights inclusion Diagnostic Test: CCR7, HLA-DR")
myfilter5 = (adata.obs['subset_source'] == 'KIR+, PBMC') | (adata.obs['subset_source'] == 'CD57+, PBMC') | (adata.obs['subset_source'] == 'NKG2A+, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'CD56bright, PBMC') | (adata.obs['subset_source'] == 'Adaptive, PBMC')
nk_adaptive_genes = ['SELL', 'KLRC2', 'B3GAT1', 'KIR2DL3']
ax = sc.pl.dotplot(adata[myfilter5], nk_adaptive_genes, groupby='subset_source', title="NK Adaptives Diagnostic Test: CD62L-NKG2C+CD57+KIR+")
CD57_marker = ['B3GAT1']
ax = sc.pl.dotplot(adata[myfilter5], CD57_marker, groupby='subset_source', title="CD57+ NK Diagnostic Test")
KIR_MARKER = ['KIR2DL3']
ax = sc.pl.dotplot(adata[myfilter5], KIR_MARKER, groupby='subset_source', title="CDKIR+ NK Diagnostic Test")
NKG2A_MARKER = ['KLRC1']
ax = sc.pl.dotplot(adata[myfilter5], NKG2A_MARKER, groupby='subset_source', title="NKG2A+ NK Diagnostic Test")
conorgenes = ['SOD1', 'SOD2', 'CD38', 'HIF1A', 'GSR', 'GSS', 'TXN', 'TXNRD1', 'TXNIP', 'TXNDC2', 'P4HB', 'ADAM17', 'TXN2', 'PRDX1', 'NOX1', 'NOXA1', 'NOX2', 'NOX3', 'NOX4', 'NOX5', 'CYBA', 'CYBB']
nk_dm_inc_genes = ['FCGR3B', 'FCGR3A', 'CX3CR1', 'CD160', 'LILRB1','KIR2DL3']
ax = sc.pl.dotplot(adata[myfilter3], nk_dm_inc_genes, groupby='subset_source', title="NK Dims inclusion Diagnostic Test: ILT-2, KIR, 'CD16', 'CX3CR1', 'CD160")
nk_br_pref_genes = ['CD2', 'ITGAX', 'CD44', 'ITGA5', 'ICAM1', 'CD55', 'CD59', 'SELL', 'CD55', 'CD59', 'CXCR3']
ax = sc.pl.dotplot(adata[myfilter3], nk_br_pref_genes, groupby='subset_source', title="NK Brights Preferential Expression")
nk_dim_pref_genes = ['CD11A']
ax = sc.pl.dotplot(adata[myfilter3], nk_dim_pref_genes, groupby='subset_source', title="NK Dims Preferential Expression")
nk_peripheral_genes = ['ITGAM']
ax = sc.pl.dotplot(adata[myfilter3], nk_peripheral_genes, groupby='subset_source', title="Peripheral NK Diagnostic Test: CD11B")
nk_adaptive_genes = ['SELL', 'KLRC2', 'B3GAT1', 'KIR2DL3']
ax = sc.pl.dotplot(adata[myfilter3], nk_adaptive_genes, groupby='subset_source', title="NK Adaptives Diagnostic Test: CD62L-NKG2C+CD57+KIR+")
# adaptives =="CD62L-NKG2C+CD57+KIR+"
ER_stress_genes = ['ABCA7', 'ABL1', 'AFF4', 'AGR2', 'AIFM1', 'ALOX5', 'AMFR', 'ANKZF1', 'APAF1', 'AQP11', 'ATF3', 'ATF4', 'ATG10', 'ATP2A1', 'ATP2A2', 'ATP2A3', 'ATXN3', 'AUP1', 'BAK1', 'BAX', 'BBC3', 'BCAP31', 'BCL2', 'BCL2L1', 'BCL2L11', 'BFAR', 'BOK', 'BRSK2', 'CALR', 'CANX', 'CASP4', 'CAV1', 'CCDC47', 'CCND1', 'CDK5RAP3', 'CEBPB', 'CFTR', 'CHAC1', 'CLU', 'COPS5', 'CREB3', 'CREB3L2', 'CREBZF', 'CTH', 'CXCL8', 'DAB2IP', 'DDIT3', 'DDRGK1', 'DDX3X', 'DERL1', 'DERL2', 'DERL3', 'DNAJB12', 'DNAJB2', 'DNAJB9', 'DNAJC10', 'DNAJC3', 'EDEM1', 'EDEM2', 'EDEM3', 'EIF2AK2', 'EIF2AK3', 'EIF2AK4', 'EIF2B5', 'EIF2S1', 'EIF4G1', 'ERLEC1', 'ERLIN1', 'ERLIN2', 'ERMP1', 'ERN1', 'ERO1A', 'ERP27', 'ERP29', 'ERP44', 'FAF1', 'FAF2', 'FAM8A1', 'FBXO2', 'FBXO27', 'FBXO44', 'FBXO6', 'FCGR2B', 'FICD', 'FLOT1', 'FOXRED2', 'GORASP2', 'GRINA', 'GSK3B', 'HERPUD1', 'HERPUD2', 'HM13', 'HSP90B1', 'HSPA1A', 'HSPA5', 'HYOU1', 'ITPR1', 'JKAMP', 'JUN', 'LPCAT3', 'LRRK2', 'MAN1A1', 'MAN1A2', 'MAN1B1', 'MAN1C1', 'MANF', 'MAP3K5', 'MARCKS', 'MBTPS1', 'MBTPS2', 'NCK1', 'NCK2', 'NFE2L2', 'NGLY1', 'NOD1', 'NOD2', 'NPLOC4', 'NR1H2', 'NR1H3', 'NRBF2', 'NUPR1', 'OPA1', 'OS9', 'P4HB', 'PARK7', 'PARP16', 'PARP6', 'PARP8', 'PDIA3', 'PDIA4', 'PDIA6', 'PIGBOS1', 'PIK3R1', 'PMAIP1', 'PML', 'PPP1R15A', 'PPP1R15B', 'PPP2CB', 'PSMC6', 'PTPN1', 'PTPN2', 'QRICH1', 'RASGRF1', 'RASGRF2', 'RCN3', 'RHBDD1', 'RHBDD2', 'RNF103', 'RNF121', 'RNF139', 'RNF185', 'RNFT1', 'RNFT2', 'SDF2L1', 'SEC16A', 'SEC61B', 'SEL1L', 'SELENOK', 'SERINC3', 'SERP1', 'SERP2', 'SESN2', 'SGTA', 'SIRT1', 'SRPX', 'STT3B', 'STUB1', 'SVIP', 'SYVN1', 'TARDBP', 'TBL2', 'THBS1', 'TMBIM6', 'TMCO1', 'TMED2', 'TMEM117', 'TMEM129', 'TMEM258', 'TMEM259', 'TMEM33', 'TMEM67', 'TMTC3', 'TMTC4', 'TMUB1', 'TMUB2', 'TMX1', 'TNFRSF10B', 'TOR1A', 'TP53', 'TRAF2', 'TRIB3', 'TRIM13', 'TRIM25', 'TXNDC12', 'UBA5', 'UBAC2', 'UBE2G2', 'UBE2J1', 'UBE2J2', 'UBE4A', 'UBE4B', 'UBQLN1', 'UBQLN2', 'UBXN1', 'UBXN4', 'UBXN6', 'UFC1', 'UFL1', 'UFM1', 'UGGT1', 'UGGT2', 'USP13', 'USP14', 'USP19', 'USP25', 'VAPB', 'VCP', 'WFS1', 'XBP1', 'YOD1']
ax = sc.pl.dotplot(adata[myfilter3], ER_stress_genes, groupby='subset_source', title="ER Stress Geneset Expression: Glioblastoma TINKS vs HC NKs")
HIF_GENES = ['AKT1', 'AKT2', 'AKT3', 'ALDOA', 'ANGPT2', 'ARNT', 'BCL2', 'CAMK2B', 'CAMK2D', 'CAMK2G', 'CDKN1A', 'CDKN1B', 'CREBBP', 'CSNK2A1', 'CUL2', 'CYBB', 'EDN1', 'EGFR', 'EGLN1', 'EGLN2', 'EGLN3', 'EIF4E', 'EIF4E2', 'EIF4EBP1', 'ENO1', 'ENO2', 'EP300', 'ERBB2', 'FLT1', 'GAPDH', 'HIF1A', 'HK1', 'HK2', 'HK3', 'HKDC1', 'HMOX1', 'IFNG', 'IFNGR1', 'IFNGR2', 'IGF1', 'IGF1R', 'IL6', 'IL6R', 'INSR', 'LDHA', 'LTBR', 'MAP2K1', 'MAP2K2', 'MAPK1', 'MAPK3', 'MKNK1', 'MKNK2', 'MTOR', 'NFKB1', 'NOS3', 'PDHA1', 'PDHB', 'PDK1', 'PFKFB3', 'PFKL', 'PGK1', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3R1', 'PIK3R3', 'PLCG1', 'PLCG2', 'PRKCA', 'PRKCB', 'RBX1', 'RELA', 'SERPINE1', 'SLC2A1', 'STAT3', 'TFRC', 'TIMP1', 'TLR4', 'VEGFA', 'VHL']
ax = sc.pl.dotplot(adata[myfilter3], HIF_GENES, groupby='subset_source', title="HIF Geneset Expression: Glioblastoma TINKS vs HC NKs")
nitrogen_metab_genes = ['AP1G1', 'CA13', 'CA2', 'CA5B', 'CA8', 'CD6', 'CPS1', 'GLUD1', 'GLUD2', 'GLUL', 'PTPRC', 'NCAM1']
ax = sc.pl.dotplot(adata[myfilter3], nitrogen_metab_genes, groupby='subset_source', title="Nitrogen Metabolism Geneset Expression: Glioblastoma TINKS vs HC NKs")
hypoxia_genes = ['ACKR3', 'ADM', 'ADORA2B', 'AK4', 'AKAP12', 'ALDOA', 'AMPD3', 'ANGPTL4', 'ANKZF1', 'ANXA2', 'ATF3', 'ATP7A', 'B3GALT6', 'BCL2', 'BGN', 'BHLHE40', 'BTG1', 'CASP6', 'CAV1', 'CCNG2', 'CDKN1A', 'CDKN1B', 'CDKN1C', 'CHST2', 'CITED2', 'COL5A1', 'CP', 'CSRP2', 'CXCR4', 'DCN', 'DDIT3', 'DDIT4', 'DTNA', 'DUSP1', 'EFNA1', 'EGFR', 'ENO1', 'ENO2', 'ENO3', 'ERO1A', 'ERRFI1', 'ETS1', 'EXT1', 'F3', 'FAM162A', 'FBP1', 'FOS', 'FOSL2', 'FOXO3', 'GAA', 'GALK1', 'GAPDH', 'GBE1', 'GCNT2', 'GLRX', 'GPC1', 'GPC3', 'GPC4', 'GPI', 'GRHPR', 'GYS1', 'HDLBP', 'HEXA', 'HK1', 'HK2', 'HMOX1', 'HS3ST1', 'HSPA5', 'IER3', 'IGFBP3', 'IL6', 'ILVBL', 'IRS2', 'ISG20', 'JMJD6', 'JUN', 'KDELR3', 'KDM3A', 'KLF6', 'KLF7', 'KLHL24', 'MAFF', 'MAP3K1', 'MIF', 'MT1E', 'MT2A', 'MXI1', 'MYH9', 'NAGK', 'NDRG1', 'NDST1', 'NDST2', 'NEDD4L', 'NFIL3', 'NOCT', 'NR3C1', 'P4HA1', 'P4HA2', 'PAM', 'PDGFB', 'PDK1', 'PDK3', 'PFKFB3', 'PFKL', 'PFKP', 'PGAM2', 'PGF', 'PGK1', 'PGM1', 'PGM2', 'PHKG1', 'PIM1', 'PLAC8', 'PLAUR', 'PLIN2', 'PNRC1', 'PPP1R15A', 'PRDX5', 'PRKCA', 'PYGM', 'RBPJ', 'RORA', 'RRAGD', 'S100A4', 'SAP30', 'SCARB1', 'SDC2', 'SDC3', 'SDC4', 'SELENBP1', 'SERPINE1', 'SIAH2', 'SLC25A1', 'SLC2A1', 'SLC2A3', 'SLC2A5', 'SLC37A4', 'SLC6A6', 'SRPX', 'TES', 'TGFB3', 'TGFBI', 'TGM2', 'TIPARP', 'TMEM45A', 'TNFAIP3', 'TPBG', 'TPD52', 'TPI1', 'TPST2', 'UGP2', 'VEGFA', 'VHL', 'VLDLR', 'WSB1', 'XPNPEP1', 'ZFP36', 'ZNF292']
ax = sc.pl.dotplot(adata[myfilter3], hypoxia_genes, groupby='subset_source', title="Hypoxia Geneset Expression: Glioblastoma TINKS vs HC NKs")
iNOS_genes = ['ACTB', 'ACTN4', 'ALOX5AP', 'AP1B1', 'APP', 'ATP1A1', 'ATP6V0C', 'AZGP1', 'CACNB3', 'CAPNS1', 'CD151', 'CD36', 'CD81', 'CFH', 'CP', 'CTSZ', 'DDB1', 'DYNC1H1', 'EEF1A1', 'EEF2', 'EIF2S2', 'EIF3B', 'EIF4G3', 'EIF5A', 'F12', 'FDX1', 'FKBP3', 'FTL', 'FZD4', 'G3BP2', 'GALK1', 'GCLM', 'GNB2', 'HPN', 'HSPB8', 'IFNAR2', 'JUP', 'LONP1', 'MAN2A1', 'MAP2K2', 'MGAT1', 'MYH9', 'NFIA', 'PCYT1A', 'PDIA3', 'PDIA4', 'PEA15', 'PPP1CA', 'PPP2R1A', 'PRPSAP1', 'RABAC1', 'RIPK1', 'RRAS', 'SCN1B', 'SERPING1', 'SLC6A8', 'SORD', 'SOX4', 'STARD10', 'STAT3', 'TSPAN4', 'UBC', 'VDAC1', 'VWF', 'YWHAQ', 'ZFP36L1']
ax = sc.pl.dotplot(adata[myfilter3], iNOS_genes, groupby='subset_source', title="iNOS Geneset Expression: Glioblastoma TINKS vs HC NKs")
mitocomp1_genes = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6']
ax = sc.pl.dotplot(adata[myfilter3], mitocomp1_genes, groupby='subset_source', title="iNOS Geneset Expression: Glioblastoma TINKS vs HC NKs")
mitocomp1_genes = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6']
ax = sc.pl.dotplot(adata[myfilter3], mitocomp2_genes, groupby='subset_source', title="iNOS Geneset Expression: Glioblastoma TINKS vs HC NKs")
mitocomp3_genes = ['UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10']
ax = sc.pl.dotplot(adata[myfilter3], mitocomp3_genes, groupby='subset_source', title="Mitochondrial Complex III Geneset Expression: Glioblastoma TINKS vs HC NKs")
mitocomp4_genes = ['COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']
ax = sc.pl.dotplot(adata[myfilter3], mitocomp4_genes, groupby='subset_source', title="Mitochondrial Complex IV Geneset Expression: Glioblastoma TINKS vs HC NKs")
mitocomp3_genes = ['UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10']
mitocomp2_genes = ['SDHA', 'SDHB', 'SDHC']
mitocomp_genes = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'SDHA', 'SDHB', 'SDHC', 'UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']

mitocomp5_genes = ['ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5F1E', 'ATP5MC1', 'ATP5MC2', 'ATP5MC3', 'ATP5ME', 'ATP5MF', 'ATP5MG', 'ATP5MJ', 'ATP5MK', 'MT-ATP6', 'MT-ATP8', 'ATP5PB', 'ATP5PD', 'ATP5PF', 'ATP5PO', 'ATP5IF1']
ax = sc.pl.dotplot(adata[myfilter3], mitocomp5_gene, groupby='subset_source', title="Mitochondrial Complex IV Geneset Expression: Glioblastoma TINKS vs HC NKs")
mitocomp1_def_genes = ["ACAD9", "COA1", "ECSIT", "FOXRED1", "NDUFA1", "NDUFA10", "NDUFA12", "NDUFA13", "NDUFA2", "NDUFA3", "NDUFA5", "NDUFA6", "NDUFA8", "NDUFAB1", "NDUFAF1", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF6", "NDUFAF7", "NDUFB1", "NDUFB10", "NDUFB11", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFV1", "NDUFV2", "NDUFV3", "NUBPL", "TIMMDC1", "TMEM126B", "TMEM186", "TMEM70"]
ax = sc.pl.dotplot(adata[myfilter3], mitocomp1_def_genes, groupby='subset_source', title="Mitochondrial Complex I Deficiency Geneset Expression: Glioblastoma TINKS vs HC NKs")

MITOFISS_genes = ["AP3B1", "AURKA", "BNIP3", "C11orf65", "COX10", "CYRIB", "DCN", "DDHD1", "DDHD2", "DNM1L", "FIS1", "GDAP1", "INF2", "KDR", "LRRK2", "MAPT", "MARCHF5", "MCU", "MFF", "MIEF1", "MIEF2", "MTFP1", "MTFR1", "MTFR1L", "MTFR2", "MUL1", "MYO19", "OPA1", "PGAM5", "PINK1", "PPARG", "PRKN", "RALA", "RALBP1", "SLC25A46", "SPIRE1", "STAT2", "TMEM135", "UCP2", "VPS35"]
PosRegMITOFISS_genes = ["AURKA", "BNIP3", "DCN", "DDHD1", "DDHD2", "DNM1L", "FIS1", "KDR", "MARCHF5", "MCU", "MFF", "MIEF1", "MIEF2", "MUL1", "PGAM5", "PINK1", "PRKN", "RALA", "RALBP1", "SPIRE1", "VPS35"]
NegMitoFiss = ["CALR", "MAPT", "MFN1", "MFN2", "PINK1", "PPARG"]

mitoFUSS= ["ADCK1", "AFG3L2", "BAK1", "BAX", "BCL2A1", "BNIP3", "CHCHD3", "FIS1", "GDAP1", "HUWE1", "MCL1", "MFF", "MFN1", "MFN2", "MIEF1", "MIEF2", "MIGA1", "MIGA2", "MTCH2", "MUL1", "OMA1", "OPA1", "PARL", "PID1", "PLD6", "PRKN", "RCC1L", "STOML2", "TFRC", "THG1L", "USP30", "VAT1", "ZDHHC6"]
NEG_REG_MITOFUSS = ["ADCK1", "BNIP3", "HUWE1", "MUL1", "OMA1", "PRKN", "TFRC", "VAT1"]
POS_REG_MITOFUSS = ["MFN1", "OPA1", "PLD6", "ZDHHC6"]

RELEASECYTOC = ["AKT1", "AVP", "BAD", "BAK1", "BAX", "BBC3", "BCL2", "BCL2A1", "BCL2L1", "BCL2L11", "BCL2L2", "BID", "BIK", "BMF", "BNIP3", "BOK", "CIDEB", "CLU", "DNM1L", "FAM162A", "FIS1", "FXN", "FZD9", "GGCT", "GHITM", "GPER1", "GPX1", "HGF", "HRK", "IFI6", "IGF1", "LMNA", "MCL1", "MFF", "MLLT11", "MMP9", "MOAP1", "NOL3", "OPA1", "PARL", "PDCD5", "PINK1", "PLAUR", "PLSCR3", "PMAIP1", "PPIF", "PRELID1", "PRKN", "PSMD10", "PYCARD", "SFN", "SOD2", "TIMM50", "TNFSF10", "TP53", "TRIAP1"]
POSREGCYTORELEASE = ["BAD", "BAK1", "BAX", "BBC3", "BCL2L11", "BID", "BIK", "BMF", "BNIP3", "CIDEB", "DNM1L", "FAM162A", "GPER1", "HRK", "MFF", "MLLT11", "MMP9", "MOAP1", "PDCD5", "PINK1", "PLAUR", "PMAIP1", "PYCARD", "TNFSF10", "TP53"]
NEGREGCYTOCRELEASE = ["AKT1", "AVP", "BAK1", "BCL2L1", "BCL2L2", "CLU", "FXN", "GHITM", "GPX1", "HGF", "IGF1", "LMNA", "NOL3", "OPA1", "PARL", "PPIF", "PRELID1", "PRKN", "PSMD10", "TRIAP1"]

POSREGDEPOL = ["ADCY10", "DCN", "KDR", "MLLT11", "MYOC", "P2RX7", "PARP1", "RACK1", "TSPO"]
POSREG_MITOPHAGY_DEPO = ["CDC37", "HDAC6", "HUWE1", "MFN2", "OPTN", "PINK1", "PRKN", "TOMM7", "VDAC1", "VPS13C"]

TGFB_genes = ["ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "AMH", "ASCL2", "ATF1", "ATF3", "ATP1B2", "ATP1B3", "BAMBI", "BMP2", "BMP7", "BMP8B", "BMPR1A", "BMPR1B", "BMPR2", "CALD1", "CD40", "CDK8", "CDKN1B", "CDKN2B", "CREBBP", "CTSC", "CUL1", "EIPR1", "E2F4", "E2F5", "EP300", "FST", "ID1", "ID2", "ID3", "ID4", "INHBA", "LRRC2", "LTBP1", "NBL1", "RBX1", "SKP1", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9", "SMURF1", "SMURF2", "SP1", "TFDP1", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGIF1", "TGIF2", "VPS53", "VPS54", "VPS52", "VPS51", "ZFYVE16", "ZFYVE9"]
ax = sc.pl.dotplot(adata[myfilter3], TGFB_genes, groupby='subset_source', title="Mitochondrial Complex I Deficiency Geneset Expression: Glioblastoma TINKS vs HC NKs")
nutrient_transporter_genes = ["SLC2A1", "SLC2A3", "SLC2A5", "SLC2A6", "SLC2A8", "SLC27A1", "SLC27A2", "SLC30A1", "SLC30A4", "SLC30A5", "SLC36A1", "SLC36A4", "SLC38A1", "SLC38A2", "SLC38A5", "SLC39A1", "SLC39A3", "SLC39A4", "SLC39A6", "SLC39A7", "SLC39A8", "SLC11A1", "SLC11A2", "SLC40A1", "SLC31A1", "SLC31A2", "ATP7A", "ABCA1", "ABCG1", "AQP3", "AQP5", "DHCR7", "CYP2R1", "CYP27B1"]
ax = sc.pl.dotplot(adata[myfilter3], nutrient_transporter_genes, groupby='subset_source', title="Mitochondrial Complex I Deficiency Geneset Expression: Glioblastoma TINKS vs HC NKs")
ROS_genes = ["NFE2L2", "SOD1", "SOD2", "SOD3", "CAT", "GPX2", "GPX3", "GPX4", "GSR", "HMOX1", "KEAP1", "TP53", "BCL2", "BAX", "HSPA1A", "HSPB1"]
ax = sc.pl.dotplot(adata[myfilter3], ROS_genes, groupby='subset_source', title="Key ROS Geneset Expression: Glioblastoma TINKS vs HC NKs")
ER_GENES = ["ATF6", "XBP1", "CALR", "CANX", "SYVN1", "SEL1L"]
ax = sc.pl.dotplot(adata[myfilter3], ER_GENES, groupby='subset_source', title="Key ER Stress Geneset Expression: Glioblastoma TINKS vs HC NKs")
NK_chemokine_funct_genes = ['LAMP', 'TNF', 'IFNG', 'NFKB1', 'NFKB2']
ax = sc.pl.dotplot(adata[myfilter3], NK_chemokine_funct_genes, groupby='subset_source', title="NK Cell Chemokine and Granularity Function Gene Expression: Glioblastoma TINKS vs HC NKs")
NK_CYTO_GENES = ["ACVR1", "BID", "CD244", "CD247", "CD48", "FCER1G", "FCGR3A", "FCGR3B", "GZMB", "HCST", "HLA-A", "HLA-C", "HLA-E", "HLA-G", "ICAM1", "ICAM2", "IFNAR1", "IFNAR2", "IFNG", "IFNGR1", "IFNGR2", "ITGAL", "ITGB2", "KIR2DL3", "KIR2DL4", "KLRC1", "KLRC2", "KLRD1", "LAT", "LCK", "LCP2", "MICB", "NCR1", "NCR3", "NFATC1", "NFATC2", "PRF1", "PTK2B", "RAC2", "RAC3", "SH2D1A", "SH2D1B", "SH3BP2", "SHC1", "SYK", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D", "TNFSF10", "TYROBP", "ULBP1", "ULBP2", "VAV1", "VAV3", "ZAP70"]
CYTOTOXICITY =  ["FASLG", "FAS", "FADD", "TNFRSF10A/TRAILR1", "TNFRSF10B/TRAILR2", "TNFRSF10C/TRAILR3", "TNFRSF10D/TRAILR4", "TNFRSF10", "TNFRSF11B/OPG", "MAPK8/JNK", "CASP8", "CFLAR3", "CASP3", "NCR3", "NCR2", "NCR1", "CD244", "CD226", "KLRF1", "CD69", "TNF", "TNFRSF1A", "TNFAIP3", "TNFRSF1B", "TNFSF9", "TNFRSF9", "TNFRSF8", "TNFRSF4", "TNIK", "TRADD", "TRAF1", "TRAF2", "TNIK", "IFNG", "IFNGR1", "IRF1", "IRF9", "JAK1", "JAK2", "MAP2K1", "IFNA1", "IDO1", "DAPK1", "MTOR", "LAMP1", "PRF1", "STAT4", "GNLY", "Gzmba", "Gzmb", "Gzmh", "Gzmm", "Gzmk", "NKG2D", "NKG2C", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5"]


ax = sc.pl.dotplot(adata[myfilter3], NK_CYTO_GENES, groupby='subset_source', title="NK Cell Cytoxicity Function Gene Expression: Glioblastoma TINKS vs HC NKs")
BRIGHTS_MARKER_GENES = ['XCL1', 'CD44', 'LMNA', 'CD52', 'SELL', 'CCR7', 'CXCR4', 'CCL3', 'CCL4L1', 'PIM1', 'NFKBIA', 'NR4A1']
ax = sc.pl.dotplot(adata[myfilter3], NK_CYTO_GENES, groupby='subset_source', title="Brights Criteria Diagnostic: Glioblastoma TINKS vs HC NKs")
INTERLEUKINS = ["IL1A" "IL1B" "IL1F10" "IL1RN" "IL2" "IL3" "IL4" "IL5" "IL6" "IL7" "CXCL8" "IL9" "IL10" "IL11" "IL12A" "IL12B" "IL13" "IL15" "IL16" "IL17A" "IL17B" "IL17C" "IL17D" "IL17F" "IL18" "IL19" "IL20" "IL21" "IL22" "IL23A" "IL24" "IL25" "IL26" "IL27" "IL31" "IL32" "IL33" "IL34" "IL36A" "IL36B" "IL36G" "IL36RN" "IL37"]
#XCLI L SELECTIN (sell), 
TGF = ["ID2", "SKP1", "CTSC", "TGFB1", "RBX1", "SMAD2", "SMAD3", "SMAD7", "SMURF1", "SMURF2", "LTBP1", "LTBP2", "LTBP3", "TGFBR1", "TGFBR2", "TGFB2", "TGFB3", "LRRC32", "MYC", "BMPR1", "BMPR2", "BMP2", "BMP7", "BMP4", "TGIF1", "TGIF2", "HDAC1", "SKI", "SKIL", "RNF111", "ARK2N", "ARK2C", "MMP2", "MMP3", "MMP9", "TIMP1", "ATF3", "JUNB", "FOXP3", "CDH1", "ITGAE", "CD9", "ITGAV", "ITGB3", "ITGB8", "ILK", "VPS51", "VPS52", "VPS53", "VPS54"]
ACTRHLA = ["HLA-A", "HLA-C", "HLA-DMA", "HLA-DOA", "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-E",  "HLA-G", "NCR3",  "NCR1", "CD244", "CD226", "KLRF1"]
inhibreceptors = ["KIR2DL4", "KIR2DL3", "KLRD1", "LILRA2", "LILRB1", "LILRB1", "LAIR1"]
CYTOTOXICITY =  ["FASLG", "FAS", "FADD", "TNFRSF10A/TRAILR1", "TNFRSF10B/TRAILR2", "TNFRSF10C/TRAILR3", "TNFRSF10D/TRAILR4", "TNFRSF10", "TNFRSF11B/OPG", "MAPK8/JNK", "CASP8", "CFLAR3", "CASP3", "NCR3", "NCR2", "NCR1", "CD244", "CD226", "KLRF1", "CD69", "TNF", "TNFRSF1A", "TNFAIP3", "TNFRSF1B", "TNFSF9", "TNFRSF9", "TNFRSF8", "TNFRSF4", "TNIK", "TRADD", "TRAF1", "TRAF2", "TNIK", "IFNG", "IFNGR1", "IRF1", "IRF9", "JAK1", "JAK2", "MAP2K1", "IFNA1", "IDO1", "DAPK1", "MTOR", "LAMP1", "PRF1", "STAT4", "GNLY", "Gzmba", "Gzmb", "Gzmh", "Gzmm", "Gzmk", "NKG2D", "NKG2C", "KIR2DS1", "KIR2DS2", "KIR2DS4", "KIR2DS5"]
ax = sc.pl.dotplot(adata[myfilter4], CYTOTOXICITY, groupby='subset_source', title="Brights Criteria Diagnostic: Glioblastoma TINKS vs HC NKs")
chemo_cytogenes = ["ADCY1", "ADCY2", "ADCY3", "ADCY4", "ADCY5", "ADCY6", "ADCY7", "ADCY8", "ADCY9", "AKT1", "AKT2", "AKT3", "ARRB1", "ARRB2", "BCAR1", "BRAF", "CCL1", "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL17", "CCL18", "CCL19", "CCL2", "CCL20", "CCL21", "CCL22", "CCL23", "CCL24", "CCL25", "CCL26", "CCL27", "CCL28", "CCL3", "CCL3L1", "CCL3L3", "CCL4", "CCL4L2", "CCL5", "CCL7", "CCL8", "CCR1", "CCR10", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CDC42", "CHUK", "CRK", "CRKL", "CSK", "CX3CL1", "CX3CR1", "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9", "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "DOCK2", "ELMO1", "FGR", "FOXO3", "GNAI1", "GNAI2", "GNAI3", "GNB1", "GNB2", "GNB3", "GNB4", "GNB5", "GNG10", "GNG11", "GNG12", "GNG13", "GNG2", "GNG3", "GNG4", "GNG5", "GNG7", "GNG8", "GNGT1", "GNGT2", "GRB2", "GRK1", "GRK2", "GRK3", "GRK4", "GRK5", "GRK6", "GRK7", "GSK3A", "GSK3B", "HCK", "HRAS", "IKBKB", "IKBKG", "ITK", "JAK2", "JAK3", "KRAS", "LYN", "MAP2K1", "MAPK1", "MAPK3", "NCF1", "NFKB1", "NFKBIA", "NFKBIB", "NRAS", "PAK1", "PARD3", "PF4", "PF4V1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R5", "PLCB1", "PLCB2", "PLCB3", "PLCB4", "PPBP", "PPBPP1", "PREX1", "PRKACA", "PRKACB", "PRKACG", "PRKCB", "PRKCD", "PRKCZ", "PRKX", "PTK2", "PTK2B", "PXN", "RAC1", "RAC2", "RAF1", "RAP1A", "RAP1B", "RASGRP2", "RELA", "RHOA", "ROCK1", "ROCK2", "SHC1", "SHC2", "SHC3", "SHC4", "SOS1", "SOS2", "STAT1", "STAT2", "STAT3", "STAT5B", "TIAM1", "TIAM2", "VAV1", "VAV2", "VAV3", "WAS", "WASL", "XCL1", "XCL2", "XCR1", "CCL2", "CCL3", "CCL4", "CCL5", "CXCL8", "CXCL10", "IL1B", "IL6", "IL7", "IL10", "IL12B", "IL2", "IL4", "IL5", "IL13", "IL15", "IL17", "CCL26", "IL2", "IL4", "IL5", "IL13", "IL-15", "IL-18"]
rosgene = ["ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC", "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB", "MBP", "MGST1", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2", "NQO1", "OXSR1", "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "SBNO2", "SCAF4", "SOD1", "SOD2", "STK25", "TXN", "TXNRD1", "TXNRD2"]
tca = ["ACLY", "ACO1", "ACO2", "CDK7", "CS", "DLST", "FH", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", "MDH1", "MDH2", "OGDH", "PC", "SUCLA2", "SUCLG1", "SUCLG2"]
OXPHOS = ["ATP6AP1", "ATP6V0A1", "ATP6V0A2", "ATP6V0B", "ATP6V0C", "ATP6V0D1", "ATP6V0E1", "ATP6V0E2", "ATP6V1A", "ATP6V1B2", "ATP6V1C1", "ATP6V1D", "ATP6V1E1", "ATP6V1E2", "ATP6V1F", "ATP6V1G1", "ATP6V1H", "CD4", "CLDN3", "COX10", "COX11", "COX15", "COX17", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7A2L", "COX7B", "COX7C", "COX8A", "CP", "CPE", "CPM", "CYC1", "LHPP", "NDUFA1", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA4L2", "NDUFA5", "NDUFA6", "NDUFA8", "NDUFA9", "NDUFAB1", "NDUFB1", "NDUFB10", "NDUFB11", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB9", "NDUFC1", "NDUFC2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "NDUFV3", "PPA1", "PPA2", "RERE", "SDHA", "SDHB", "SDHC", "TCIRG1", "TPP1", "UQCR10", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ"]
ppp = ["AHR", "CHEK1", "DERA", "G6PD", "GLYCTK", "H6PD", "IDNK", "PGD", "PGLS", "PRPS1", "PRPS2", "RBKS", "RPE", "RPIA", "TALDO1", "TKT"]
gly_ser_thre_genes = ["AHCY", "ALAS1", "AOC3", "CDH5", "CHDH", "CTH", "DMGDH", "GAMT", "GATM", "GCAT", "GNMT", "JAG1", "MAOB", "PHGDH", "PIPOX", "PSAT1", "PSPH", "SARDH", "SDS", "SDSL", "SRR"]
ala_asp_glut_genes = ["ABAT", "ALDH4A1", "ALDH5A1", "AP1G1", "ASL", "ASNS", "ASPA", "ASS1", "ATF2", "CAD", "CASP2", "CD6", "CDKN1C", "CPS1", "DDO", "GFPT1", "GFPT2", "GLS", "GLS2", "GLUD1", "GLUD2", "GOT1", "GOT2", "GPT2", "IL4I1", "NIT2", "PPAT", "RIMKLA", "RIMKLB"]
oxyst_genes = ["ACOT11", "ACOT13", "ACOT2", "ACOT4", "ACOT7", "ACOT8", "ACOT9", "AMACR", "CH25H", "CYP27A1", "CYP46A1", "CYP7B1", "DBP", "DHCR7", "EBP", "EPHX2", "ESR2", "GPR183", "HSD11B1", "HSD11B2", "HSD3B7", "INSIG1", "NR1H2", "RORC", "SCP2", "SLC27A2", "SLC27A5", "SULT2B1"]
protgly_genes = ["AKT2", "AKT3", "ANK1", "ANK3", "ARAF", "ARHGEF1", "ARHGEF12", "BRAF", "CAMK2B", "CAMK2D", "CAMK2G", "CASP3", "CAV1", "CAV2", "CBL", "CBLB", "CCND1", "CD44", "CD63", "CDC42", "CDKN1A", "CTNNB1", "CTSL", "CTTN", "CYCS", "DCN", "DDX5", "DROSHA", "EGFR", "EIF4B", "ELK1", "ERBB2", "ERBB3", "ESR1", "EZR", "FAS", "FASLG", "FGFR1", "FLNA", "FLNB", "FN1", "FRS2", "FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "GAB1", "GPC1", "GPC3", "GRB2", "HBEGF", "HCLS1", "HIF1A", "HPSE", "HRAS", "HSPG2", "IGF1", "IGF1R", "IQGAP1", "ITGA2", "ITGA5", "ITGAV", "ITGB1", "ITPR1", "ITPR2", "ITPR3", "KRAS", "LUM", "MAP2K1", "MAP2K2", "MAPK1", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "MAPK3", "MDM2", "MET", "MMP2", "MMP9", "MRAS", "MSN", "MTOR", "MYC", "NRAS", "NUDT16L1", "PAK1", "PDCD4", "PDPK1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R3", "PLAU", "PLAUR", "PLCE1", "PLCG1", "PLCG2", "PPP1CA", "PPP1CB", "PPP1CC", "PPP1R12A", "PPP1R12B", "PPP1R12C", "PRKACA", "PRKACB", "PRKCA", "PRKCB", "PTCH1", "PTK2", "PTPN11", "PTPN6", "PXN", "RAC1", "RAF1", "RDX", "RHOA", "ROCK1", "ROCK2", "RRAS", "RRAS2", "SDC1", "SDC2", "SDC4", "SEMA4F", "SMAD2", "SOS1", "SOS2", "SRC", "STAT3", "TFAP4", "TGFB1", "TGFB2", "THBS1", "TIAM1", "TIMP3", "TLR2", "TLR4", "TP53", "VAV2", "VEGFA", "WNT10A", "WNT2B", "WNT3", "WNT5B", "WNT7B", "WNT9A"]
ketone_genes = ["BDH1", "BDH2", "CHM", "CHN1", "HMGCL", "HMGCS1", "OXCT1"]
glyx_dicarbxyl = ["AFMID", "AMT", "CAT", "GCSH", "GLUL", "HYI", "MCEE", "PCCA", "PCCB", "PGP", "SHMT1", "SHMT2"]
fruct_mann_genes = ["ADD1", "AKR1B1", "ENOSF1", "FPGT", "GMDS", "GMPPA", "GMPPB", "KHK", "MPI", "PFKFB2", "PFKFB3", "PFKFB4", "PMM1", "PMM2", "SORD", "TIGAR", "TKFC"]
glyclip_genes = ["ADRB2", "AGK", "AGPAT1", "AGPAT2", "AGPAT3", "AGPAT4", "AGPAT5", "CHI3L2", "COL5A2", "DGAT1", "DGAT2", "DGKA", "DGKD", "DGKE", "DGKG", "DGKH", "DGKQ", "DGKZ", "GK", "GLA", "GPAM", "GPAT2", "GPAT3", "GPAT4", "LCLAT1", "LPIN1", "LPIN2", "LPL", "MBOAT1", "MBOAT2", "MGLL", "MOGAT1", "PLPP1", "PLPP2", "PLPP3", "PLPP5", "PNPLA2"]
glycophos_genes = ["ACHE", "ADPRM", "ADSL", "AKT1", "BMP6", "CA11", "CDIPT", "CDS1", "CDS2", "CEPT1", "CHKA", "CHKB", "CHPT1", "CRLS1", "DSP", "ETNK1", "ETNK2", "GNPAT", "GPCPD1", "GPD1L", "GPD2", "LCAT", "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "LPGAT1", "LYPLA1", "LYPLA2", "MBOAT7", "PCYT1A", "PCYT2", "PEMT", "PGS1", "PISD", "PLA2G10", "PLA2G12A", "PLA2G15", "PLA2G4A", "PLA2G4C", "PLA2G6", "PLB1", "PLD1", "PLD2", "PLD3", "PLD4", "PNPLA6", "PNPLA7", "PTDSS1", "PTDSS2"]
MAT_NK = ['TCF7']
ax = sc.pl.dotplot(adata[myfilter2], MAT_NK, groupby='subset_source', title="HC-NK Cell Diagnostic Test: CD57, CD56, CD158b1/b2 or KIR2DL1, NKG2C, NKG2A Expression")
nitrogen_metabolism = ["AMT", "ASNS", "CA2", "CA5B", "CA8", "CPS1", "CTH", "GLS", "GLS2", "GLUD1", "GLUD2", "GLUL"]

import scanpy as sc
import pandas as pd
import os

# Load the data
file_path = r"C:\Users\aroche6\abc\all_nk_cells.h5ad"
print("Loading data from:", file_path)
adata = sc.read_h5ad(file_path)
print("Data loaded successfully")

# Define the filter for specific cell subsets
myfilter4 = (
    (adata.obs['subset_source'] == 'KIR+, PBMC') | 
    (adata.obs['subset_source'] == 'CD57+, PBMC') | 
    (adata.obs['subset_source'] == 'NKG2A+, PBMC') | 
    (adata.obs['subset_source'] == 'CD56bright, PBMC') | 
    (adata.obs['subset_source'] == 'Adaptive, PBMC') | 
    (adata.obs['subset_source'] == 'CD56bright, glioblastoma_tumor') | 
    (adata.obs['subset_source'] == 'CD56dim, glioblastoma_tumor')
)

# Apply the filter
filtered_adata = adata[myfilter4]
print("Filtering done")

# Verify the filtering
print("Number of cells after filtering:", filtered_adata.shape[0])

# Define the genes of interest
nk_exc_genes = ['CD3D', 'CD3E', 'CD3G', 'CD14', 'CD8A']

# Initialize an empty list to store the results
results = []

# Calculate mean expression values for each gene in each subset
for subset in filtered_adata.obs['subset_source'].unique():
    subset_data = filtered_adata[filtered_adata.obs['subset_source'] == subset]
    mean_expression = subset_data[:, nk_exc_genes].X.mean(axis=0).tolist()
    results.append([subset] + mean_expression)
    print(f"Processed subset: {subset}")

# Create a pandas DataFrame with the results
columns = ['Subset'] + nk_exc_genes
df = pd.DataFrame(results, columns=columns)
print("DataFrame created")
print(df)

# Define the output file path (change the directory to ensure write permissions)
output_file = r"C:\Users\aroche6\abc\mean_expression_values.xlsx"
alternative_output_file = r"C:\mean_expression_values.xlsx"

# Ensure the directory exists
os.makedirs(os.path.dirname(output_file), exist_ok=True)
print("Directory checked/created")

# Save the DataFrame to an Excel file (try both paths)
try:
    df.to_excel(output_file, index=False)
    print(f"Mean expression values saved to {output_file}")
    # Check if the file exists
    if os.path.isfile(output_file):
        print("File saved successfully and found at the specified location.")
    else:
        print("File was not found at the specified location after attempting to save.")
except Exception as e:
    print(f"An error occurred while saving the file to {output_file}: {e}")

# Try saving to an alternative path
try:
    df.to_excel(alternative_output_file, index=False)
    print(f"Mean expression values saved to {alternative_output_file}")
    # Check if the file exists
    if os.path.isfile(alternative_output_file):
        print("File saved successfully and found at the alternative location.")
    else:
        print("File was not found at the alternative location after attempting to save.")
except Exception as e:
    print(f"An error occurred while saving the file to {alternative_output_file}: {e}")

import pandas as pd
df.to_csv('data.csv', index=False)

df.to_excel('data.xlsx', index=False)

output_file = os.path.join('C:/My Documents/Python/Dataset/ny_station.xlsx')
writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
writer._save()

ax = sc.pl.dotplot(adata[myfilter14], trans, groupby='subset_source', title="Glut1 and GLut3"