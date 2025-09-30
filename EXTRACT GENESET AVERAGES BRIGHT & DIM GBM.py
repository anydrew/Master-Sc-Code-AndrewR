# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:37:02 2025

@author: AROCHE6
"""

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

# Define the file path and gene lists
file_path = r"C:\Users\aroche6\abc\adata_all_nk_after_mapping.h5ad"
ros_genes = [
    "ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC", 
    "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB", 
    "LAMTOR5", "MBP", "MGST1", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2", "NQO1", "OXSR1", 
    "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "SBNO2", "SCAF4", 
    "SOD1", "SOD2", "STK25", "TXN", "TXNRD1", "TXNRD2"
]
tgfbeta_genes = [
    "CD9", "ATF3", "SIN3A", "MAPK11", "DCN", "ENG", "MAPK9", "HFE", "CDK8", "PPP1R15A", 
    "BMP8B", "FNTA", "ATP1B2", "CITED2", "IFNGR2", "LTBP2", "SMAD1", "BMPR1B", "NCOR1", 
    "TAB1", "E2F5", "LTBP3", "PPM1A", "XIAP", "SMAD3", "MAPK12", "ID1", "TRIM33", 
    "SERPINE1", "GDF15", "BMP2", "THBS1", "CREBBP", "CDK9", "CTSC", "IFNG", "BCAR3", 
    "RGMB", "EMP3", "GDF11", "FURIN", "WWTR1", "SKIL", "MAPK14", "RGMA", "PPP2CA", 
    "EP300", "FKBP1A", "VPS51", "ACVR2A", "ACVR1B", "E2F4", "CDKN1C", "SLC20A1", 
    "BMPR1A", "LTBP1", "TGFBR2", "RHOA", "ACVR1", "SMAD4", "INHBA", "ID3", "PPP1CA", 
    "MMP9", "CALD1", "JUNB", "SMURF1", "SP1", "TGFB3", "ZFYVE9", "NRAS", "NBL1", 
    "MAPK13", "NRROS", "SNIP1", "BAMBI", "TFRC", "MAPK1", "MYC", "NCOR2", "GREM1", 
    "VPS54", "SMAD9", "UBE2D3", "HDAC1", "SMURF2", "ITGB8", "ID2", "APC", "ACVRL1", 
    "TIMP1", "PPP2CB", "CTNNB1", "BMP6", "RRAS", "BMPR2", "SPTBN1", "KLF10", "AMH", 
    "ATF2", "FOSL1", "MAPK8", "HRAS", "SKP1", "JUN", "TFDP1", "TGFB1", "RAB31", 
    "ATP1B3", "ITGAV", "RNF111", "LRRC32", "SMAD6", "ILK", "MAPK3", "TFR2", "PPP2R1A", 
    "SMAD2", "BMP7", "TGIF2", "MAP3K7", "VPS53", "ZFYVE16", "TGFB2", "ARID4B", 
    "ITGAE", "JUND", "DCP1B", "TJP1", "MMP2", "ACVR1C", "ROCK1", "MAPK10", "CDKN2B", 
    "MAP3K7CL", "FBN1", "BMP1", "NEO1", "HDAC2", "TGIF1", "PMEPA1", "CDKN1B", "RBX1", 
    "SKI", "FST", "ATF1", "LRRC2", "CUL1", "SMAD7", "ASCL2", "ID4", "TGFBR1", "ACVR2B", 
    "CD40", "PPP2R1B", "SMAD5", "HIPK2", "THSD4"
]

# Define Oxidative Stress gene list
oxidative_stress_genes = [
    'ABCC1', 'ABCD1', 'ABL1', 'ADAM9', 'AGAP3', 'AIF1', 'AIFM1', 'AIFM2', 'AKT1', 'ALAD', 
    'ALDH3B1', 'ALOX5', 'ALS2', 'ANKZF1', 'ANXA1', 'APOD', 'APOE', 'APP', 'AREG', 'ARL6IP5', 
    'ARNT', 'ATF2', 'ATF4', 'ATM', 'ATOX1', 'ATP13A2', 'ATP2A2', 'ATP7A', 'AXL', 'BAK1', 
    'BANF1', 'BCL2', 'BECN1', 'BMP7', 'BNIP3', 'BRF2', 'BTK', 'CAMKK2', 'CAPN2', 'CASP3', 
    'CAT', 'CBX8', 'CCS', 'CD36', 'CD38', 'CDK1', 'CHCHD2', 'CHUK', 'COL1A1', 'CPEB2', 
    'CRK', 'CRYAB', 'CYB5B', 'CYGB', 'CYP1B1', 'DAPK1', 'DDR2', 'DHCR24', 'DHFR', 'DUOX1', 
    'ECT2', 'EDN1', 'EGFR', 'EIF2S1', 'EPAS1', 'ERCC1', 'ERCC2', 'ERCC3', 'ERCC6', 
    'ERCC6L2', 'ERCC8', 'ERMP1', 'ERN1', 'ERO1A', 'ETFDH', 'ETV5', 'EZH2', 'FANCC', 
    'FANCD2', 'FBLN5', 'FER', 'FKBP1B', 'FOS', 'FOSL1', 'FOXO1', 'FOXO3', 'FOXO4', 
    'FOXP1', 'FUT8', 'FXN', 'FYN', 'G6PD', 'GCH1', 'GCLC', 'GCLM', 'GGT7', 'GJB2', 
    'GLRX2', 'GPX2', 'GPX3', 'GPX4', 'GPX7', 'GPX8', 'GSR', 'GSS', 'GSTP1', 'HDAC2', 
    'HDAC6', 'HIF1A', 'HM13', 'HMOX1', 'HMOX2', 'HSF1', 'HSPA1A', 'HSPA1B', 'HTRA2', 
    'HYAL1', 'HYAL2', 'IDH1', 'IL18BP', 'IL18RAP', 'IL1A', 'IL6', 'IPCEF1', 'JAK2', 
    'JUN', 'KAT2B', 'KDM6B', 'KEAP1', 'KLF2', 'LIAS', 'LONP1', 'LRRK2', 'MAP1LC3A', 
    'MAP2K4', 'MAP3K5', 'MAPK1', 'MAPK13', 'MAPK3', 'MAPK7', 'MAPK8', 'MAPK9', 
    'MAPKAP1', 'MAPT', 'MCTP1', 'MDM2', 'MET', 'MGAT3', 'MGST1', 'MICB', 'MITF', 
    'MLH1', 'MMP7', 'MPO', 'MSRA', 'MSRB1', 'MTHFD1', 'MTHFD2', 'MTHFR', 'NADK', 
    'NDUFA9', 'NDUFS1', 'NDUFS2', 'NFE2L1', 'NFE2L2', 'NLRP3', 'NOX4', 'NQO1', 
    'NRG1', 'NRG2', 'OGG1', 'OLR1', 'OTUD1', 'OXSR1', 'P4HB', 'PA2G4', 'PARK7', 
    'PDXK', 'PGD', 'PGK1', 'PGM3', 'PIK3R1', 'PIK3R5', 'PKM', 'PLA2G4A', 'PLA2G5', 
    'PLCB1', 'PLCB3', 'PLCG1', 'PLCG2', 'PLK1', 'PML', 'PNP', 'POLB', 'POLG', 
    'PPARGC1A', 'PRDX1', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6', 'PRKAA1', 
    'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKCQ', 
    'PXDN', 'RAB20', 'RAD23A', 'RAP1B', 'RB1', 'RBM38', 'RCAN1', 'RNF4', 'ROS1', 
    'RUNX2', 'SIRT1', 'SIRT3', 'SIRT5', 'SLC11A2', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
    'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A1', 
    'SLC25A12', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
    'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC25A4', 'SLC25A5', 'SLC25A6', 
     'SLC25A12', 'SLC25A15', 'SLC25A19', 
    'SLC25A22', 'SLC25A25', 'SLC25A29', 'SLC25A32', 'SLC25A38', 'SLC25A39', 
    'SOD1', 'SOD2', 'SOD3', 'SRXN1', 'STK11', 'TBXAS1', 'TEK', 'TFAM', 'TGFB1', 
    'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TGFBR3', 'TLR2', 'TLR4', 
    'TNF', 'TNFAIP3', 'TNFRSF1A', 'TNFRSF1B', 'TP53', 'TXN', 'TXNDC5', 'TXNIP', 
    'UQCRC1', 'UQCRC2', 'UQCRFS1', 'VHL', 'XPC', 'ZFP36'
]
OXPHOS_genes = ['ATP1B1', 'COX10', 'UQCRH', 'ATP6V1B2', 'NDUFA5', 'NDUFA4', 'ATP6V1C1', 'NDUFC2', 'CLDN3', 'COX6C', 'LHPP', 'UQCRFS1', 'NDUFB2', 'NDUFA10', 'NDUFA4L2', 'NDUFS5', 'CD4', 'NDUFA13', 'RERE', 'UQCRQ', 'UQCRB', 'COX17', 'NDUFB7', 'COX11', 'NDUFS1', 'NDUFS4', 'ATP6V0B', 'COX7A2L', 'ATP6V0C', 'TPP1', 'NDUFB3', 'ATP6V1D', 'SURF1', 'CYC1', 'NDUFB5', 'NDUFV2', 'COX6A1', 'COX4I1', 'COX8A', 'NDUFB11', 'NDUFA12', 'COX7A2', 'NDUFA2', 'NDUFC1', 'PPA1', 'NDUFA9', 'NDUFB6', 'NDUFV1', 'NDUFB1', 'ATP6V0E2', 'COX7C', 'ATP6V1E1', 'COX5A', 'ATP6V1G1', 'COX5B', 'NDUFS6', 'ATP6V0E1', 'CPM', 'COX7B', 'ATP6V1A', 'CP', 'NDUFAB1', 'UQCRC2', 'NDUFA11', 'NDUFS8', 'NDUFS3', 'NDUFA6', 'ATP6V1E2', 'UQCR10', 'NDUFA3', 'COX15', 'NDUFB9', 'SDHA', 'NDUFA8', 'ATP6V0A2', 'SDHC', 'CPE', 'ATP6V1H', 'NDUFV3', 'ATP6V0D1', 'NDUFS7', 'NDUFA1', 'NDUFB10', 'TCIRG1', 'NDUFS2', 'UQCRC1', 'ATP6V1F', 'PPA2', 'NDUFB4', 'ATP6AP1', 'ATP6V0A1', 'SDHB', 'COX6B1']
oxphos_ratelim = ["COX10", "COX11", "COX14", "COX15", "COX16", "COX17", "COX18", "COX19", "COX20", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7A2L", "COX7B", "COX7C", "COX8A", "SURF1", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "HIGD1A", "HIGD2A"]

GLYCOL_GLUCONEOGEN_genes = ['ADORA2B', 'AGL', 'AGRN', 'AK3', 'AK4', 'AKR1A1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALG1', 'ANG', 'ANGPTL4', 'ANKZF1', 'ARPP19', 'AURKA', 'B3GALT6', 'B3GAT1', 'B3GAT3', 'B4GALT1', 'B4GALT2', 'B4GALT4', 'B4GALT7', 'BIK', 'BPNT1', 'CACNA1H', 'CAPN5', 'CASP6', 'CD44', 'CDK1', 'CENPA', 'CHPF', 'CHPF2', 'CHST1', 'CHST12', 'CHST2', 'CITED2', 'CLDN3', 'CLN6', 'COG2', 'COL5A1', 'COPB2', 'CTH', 'CXCR4', 'CYB5A', 'DCN', 'DDIT4', 'DEPDC1', 'DLD', 'DSC2', 'ECD', 'EGFR', 'EGLN3', 'ELF3', 'ENO1', 'ENO2', 'ERO1A', 'EXT1', 'EXT2', 'FAM162A', 'FKBP4', 'FUT8', 'G6PD', 'GALE', 'GALK1', 'GALK2', 'GCLC', 'GFPT1', 'GLCE', 'GLRX', 'GMPPA', 'GMPPB', 'GNE', 'GNPDA1', 'GOT1', 'GOT2', 'GPC1', 'GPC3', 'GPC4', 'GUSB', 'GYS1', 'HAX1', 'HDLBP', 'HK2', 'HMMR', 'HOMER1', 'HS2ST1', 'HS6ST2', 'HSPA5', 'IDH1', 'IDUA', 'IER3', 'IGFBP3', 'IL13RA1', 'IRS2', 'ISG20', 'KDELR3', 'KIF20A', 'KIF2A', 'LDHA', 'LHPP', 'MDH1', 'MDH2', 'ME2', 'MED24', 'MET', 'MIF', 'MPI', 'MXI1', 'NANP', 'NASP', 'NDUFV3', 'NOL3', 'NSDHL', 'NT5E', 'P4HA1', 'P4HA2', 'PAM', 'PAXIP1', 'PC', 'PDK3', 'PFKM', 'PFKL', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGLS', 'PGM2', 'PHKA2', 'PKM', 'PKP2', 'PLOD1', 'PLOD2', 'PMM2', 'POLR3K', 'PPIA', 'PPP2CB', 'PRPS1', 'PSMC4', 'PYGB', 'PYGL', 'QSOX1', 'RBCK1', 'RPE', 'RRAGD', 'SAP30', 'SDC1', 'SDC2', 'SDC3', 'SDHC', 'SLC16A3', 'SLC25A13', 'SLC35A3', 'SLC37A4', 'SOD1', 'SPAG4', 'SRD5A3', 'STMN1', 'TALDO1', 'TFF3', 'TGFA', 'TGFBI', 'TPBG', 'TPI1', 'TPST1', 'TXN', 'UGP2', 'VCAN', 'VEGFA', 'VLDLR', 'XYLT2', 'ZNF292']
Glycolysis_ratelim = ["SLC2A1", "HK1", "HK2", "HK3", "PFKM", "PFKL", "PFKP", "SLC16A4", "SLC16A1"]

comp1 = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6']
comp2 = ['SDHA', 'SDHB', 'SDHC']
comp3 = ['UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10']
comp4 = ['COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']

mitofission_genes = ['VPS35', 'SLC25A46', 'MTFR1', 'COX10', 'TMEM135', 'DCN', 'AURKA', 'MTFR2', 'MUL1', 'MTFP1', 'AP3B1', 'DDHD1', 'MIEF1', 'MCU', 'FIS1', 'MFF', 'PINK1', 'MYO19', 'DDHD2', 'OPA1', 'RALBP1', 'RALA', 'INF2', 'BAX', 'CCAR2', 'UCP2', 'LRRK2', 'STAT2', 'SPIRE1', 'ATG3', 'PGAM5', 'LPIN1', 'MAPT', 'GDAP1', 'MIEF2', 'BNIP3', 'MTFR1L', 'DNM1L', 'PPARG']
mitofussion_genes = ['MFN1', 'AFG3L2', 'MFN2', 'ADCK1', 'PID1', 'MUL1', 'STOML2', 'MIEF1', 'THG1L', 'ZDHHC6', 'FIS1', 'MFF', 'MTCH2', 'HUWE1', 'OPA1', 'BAK1', 'USP30', 'BAX', 'MCL1', 'GDAP1', 'VAT1', 'MIEF2', 'BNIP3', 'PLD6', 'CHCHD3', 'OMA1', 'PARL', 'TFRC', 'BCL2A1']

cytochrome_c_release_apop_genes = ['AVEN', 'BBC3', 'BCL2L1', 'BMF', 'CASP3', 'PMAIP1', 'XIAP', 'PINK1', 'CARD8', 'MLLT11', 'NOL3', 'CASP9', 'BIRC2', 'AIFM1', 'CYCS', 'APIP', 'DNM1L', 'PARL', 'BCL2A1', 'BIRC3', 'UACA', 'ENDOG', 'BIK', 'DFFA', 'DIABLO', 'CIDEB', 'SOD2', 'FIS1', 'MFF', 'BAK1', 'BAX', 'GGCT', 'MMP9', 'FAM162A', 'BNIP3', 'TNFSF10', 'MAPK1', 'PPIF', 'APAF1', 'MOAP1', 'DFFB', 'TP53', 'BCL2', 'CASP8', 'IGF1', 'PDCD5', 'LMNA', 'CASP6', 'FXN', 'PSMD10', 'PYCARD', 'MAPK3', 'TIMM50', 'GHITM', 'BAD', 'BOK', 'PLAUR', 'PRELID1', 'SFN', 'OPA1', 'BCL2L2', 'MCL1', 'BID', 'CLU', 'AKT1', 'IFI6', 'CASP7', 'BCL2L11', 'TRIAP1']
mitophagy_genes = ['ATG12', 'ATG5', 'CSNK2A1', 'CSNK2A2', 'CSNK2B', 'FUNDC1', 'MAP1LC3A', 'MAP1LC3B', 'MFN1', 'MFN2', 'MTERF3', 'PGAM5', 'PINK1', 'SQSTM1', 'SRC', 'TOMM20', 'TOMM22', 'TOMM40', 'TOMM5', 'TOMM7', 'UBB', 'UBC', 'ULK1', 'VDAC1']

# Load the AnnData object
adata = sc.read_h5ad(file_path)


# Filter the data for the specified subsets
myfilter12 = (adata.obs['subset_source'] == 'TiCD56bright_glioblastoma') | (adata.obs['subset_source'] == 'TiCD56dim_glioblastoma')
filtered_adata = adata[myfilter12]

# Check if the genes are in the dataset
available_ros_genes = [gene for gene in ros_genes if gene in filtered_adata.var_names]
available_tgfbeta_genes = [gene for gene in tgfbeta_genes if gene in filtered_adata.var_names]
available_oxidative_stress_genes = [gene for gene in oxidative_stress_genes if gene in filtered_adata.var_names]
available_OXPHOS_genes = [gene for gene in OXPHOS_genes if gene in filtered_adata.var_names]
available_oxphos_ratelim = [gene for gene in oxphos_ratelim if gene in filtered_adata.var_names]
available_GLYCOL_GLUCONEOGEN_genes = [gene for gene in GLYCOL_GLUCONEOGEN_genes if gene in filtered_adata.var_names]
available_Glycolysis_ratelim = [gene for gene in Glycolysis_ratelim if gene in filtered_adata.var_names]
available_comp1 = [gene for gene in comp1 if gene in filtered_adata.var_names]
available_comp2 = [gene for gene in comp2 if gene in filtered_adata.var_names]
available_comp3 = [gene for gene in comp3 if gene in filtered_adata.var_names]
available_comp4 = [gene for gene in comp4 if gene in filtered_adata.var_names]
available_mitofission_genes = [gene for gene in mitofission_genes if gene in filtered_adata.var_names]
available_mitofussion_genes = [gene for gene in mitofussion_genes if gene in filtered_adata.var_names]
available_cytochrome_c_release_apop_genes = [gene for gene in cytochrome_c_release_apop_genes if gene in filtered_adata.var_names]
available_mitophagy_genes = [gene for gene in mitophagy_genes if gene in filtered_adata.var_names]

if len(available_ros_genes) == 0 or len(available_tgfbeta_genes) == 0 or len(available_oxidative_stress_genes) == 0 or len(available_OXPHOS_genes) == 0 or len(available_oxphos_ratelim) == 0 or len(available_GLYCOL_GLUCONEOGEN_genes) == 0 or len(available_Glycolysis_ratelim) == 0 or len(available_comp1) == 0 or len(available_comp2) == 0 or len(available_comp3) == 0 or len(available_comp4) == 0 or len(available_mitofission_genes) == 0 or len(available_mitofussion_genes) == 0 or len(available_cytochrome_c_release_apop_genes) == 0 or len(available_mitophagy_genes) == 0:
    raise ValueError("None of the specified genes are available in the dataset.")

# Subset the data for the available ROS and TGF-beta genes
# Subset and analyze ROS genes (including overlaps)
adata_ros = filtered_adata[:, available_ros_genes]
adata_oxidative_stress = filtered_adata[:, available_oxidative_stress_genes]

# Subset and analyze TGF-beta genes (including overlaps)
adata_tgfbeta = filtered_adata[:, available_tgfbeta_genes]


adata_OXPHOS = filtered_adata[:, available_OXPHOS_genes]

# Subset and analyze TGF-beta genes (including overlaps)
adata_oxphos_ratelim = filtered_adata[:, available_oxphos_ratelim]

adata_GLYCOL_GLUCONEOGEN = filtered_adata[:, available_GLYCOL_GLUCONEOGEN_genes]

adata_Glycolysis_ratelim = filtered_adata[:, available_Glycolysis_ratelim]

adata_comp1 = filtered_adata[:, available_comp1]

# Subset and analyze TGF-beta genes (including overlaps)
adata_comp2 = filtered_adata[:, available_comp2]

adata_comp3 = filtered_adata[:, available_comp3]

adata_comp4 = filtered_adata[:, available_comp4]

adata_mitofission = filtered_adata[:, available_mitofission_genes]

# Subset and analyze TGF-beta genes (including overlaps)
adata_mitofussion = filtered_adata[:, available_mitofussion_genes]

adata_cytochrome_c_release_apop = filtered_adata[:, available_cytochrome_c_release_apop_genes]

# Subset and analyze TGF-beta genes (including overlaps)
adata_mitophagy = filtered_adata[:, available_mitophagy_genes]


# Calculate cumulative expression and perform statistical tests separately
mean_expression_ros_gbm = adata_ros[adata_ros.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_ros_pbmc = adata_ros[adata_ros.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_tgfbeta_gbm = adata_tgfbeta[adata_tgfbeta.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_tgfbeta_pbmc = adata_tgfbeta[adata_tgfbeta.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_oxidative_stress_gbm = adata_oxidative_stress[adata_oxidative_stress.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_oxidative_stress_pbmc = adata_oxidative_stress[adata_oxidative_stress.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

# Calculate cumulative expression and perform statistical tests separately
mean_expression_OXPHOS_gbm = adata_OXPHOS[adata_OXPHOS.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_OXPHOS_pbmc = adata_OXPHOS[adata_OXPHOS.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_oxphos_ratelim_gbm = adata_oxphos_ratelim[adata_oxphos_ratelim.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_oxphos_ratelim_pbmc = adata_oxphos_ratelim[adata_oxphos_ratelim.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_GLYCOL_GLUCONEOGEN_gbm = adata_GLYCOL_GLUCONEOGEN[adata_GLYCOL_GLUCONEOGEN.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_GLYCOL_GLUCONEOGEN_pbmc = adata_GLYCOL_GLUCONEOGEN[adata_GLYCOL_GLUCONEOGEN.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_Glycolysis_ratelim_gbm = adata_Glycolysis_ratelim[adata_Glycolysis_ratelim.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_Glycolysis_ratelim_pbmc = adata_Glycolysis_ratelim[adata_Glycolysis_ratelim.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_comp1_gbm = adata_comp1[adata_comp1.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_comp1_pbmc = adata_comp1[adata_comp1.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_comp2_gbm = adata_comp2[adata_comp2.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_comp2_pbmc = adata_comp2[adata_comp2.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_comp3_gbm = adata_comp3[adata_comp3.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_comp3_pbmc = adata_comp3[adata_comp3.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_comp4_gbm = adata_comp4[adata_comp4.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_comp4_pbmc = adata_comp4[adata_comp4.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

# Calculate cumulative expression and perform statistical tests separately
mean_expression_mitofission_gbm = adata_mitofission[adata_mitofission.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_mitofission_pbmc = adata_mitofission[adata_mitofission.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_mitofussion_gbm = adata_mitofussion[adata_mitofussion.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_mitofussion_pbmc = adata_mitofussion[adata_mitofussion.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_cytochrome_c_release_apop_gbm = adata_cytochrome_c_release_apop[adata_cytochrome_c_release_apop.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_cytochrome_c_release_apop_pbmc = adata_cytochrome_c_release_apop[adata_cytochrome_c_release_apop.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)

mean_expression_mitophagy_gbm = adata_mitophagy[adata_mitophagy.obs['subset_source'] == 'TiCD56bright_glioblastoma'].X.mean(axis=1)
mean_expression_mitophagy_pbmc = adata_mitophagy[adata_mitophagy.obs['subset_source'] == 'TiCD56dim_glioblastoma'].X.mean(axis=1)


# Continue with your analysis as before

# Perform Mann-Whitney U test
u_stat_ros, p_val_ros = mannwhitneyu(mean_expression_ros_gbm, mean_expression_ros_pbmc, alternative='two-sided')
u_stat_tgfbeta, p_val_tgfbeta = mannwhitneyu(mean_expression_tgfbeta_gbm, mean_expression_tgfbeta_pbmc, alternative='two-sided')
u_stat_oxidative_stress, p_val_oxidative_stress = mannwhitneyu(mean_expression_oxidative_stress_gbm, mean_expression_oxidative_stress_pbmc, alternative='two-sided')
# Perform Mann-Whitney U test
u_stat_OXPHOS, p_val_OXPHOS = mannwhitneyu(mean_expression_OXPHOS_gbm, mean_expression_OXPHOS_pbmc, alternative='two-sided')
u_stat_oxphos_ratelim, p_val_oxphos_ratelim = mannwhitneyu(mean_expression_oxphos_ratelim_gbm, mean_expression_oxphos_ratelim_pbmc, alternative='two-sided')
u_stat_GLYCOL_GLUCONEOGEN, p_val_GLYCOL_GLUCONEOGEN = mannwhitneyu(mean_expression_GLYCOL_GLUCONEOGEN_gbm, mean_expression_GLYCOL_GLUCONEOGEN_pbmc, alternative='two-sided')
u_stat_Glycolysis_ratelim, p_val_Glycolysis_ratelim = mannwhitneyu(mean_expression_Glycolysis_ratelim_gbm, mean_expression_Glycolysis_ratelim_pbmc, alternative='two-sided')
u_stat_comp1, p_val_comp1 = mannwhitneyu(mean_expression_comp1_gbm, mean_expression_comp1_pbmc, alternative='two-sided')
u_stat_comp2, p_val_comp2 = mannwhitneyu(mean_expression_comp2_gbm, mean_expression_comp2_pbmc, alternative='two-sided')
u_stat_comp3, p_val_comp3 = mannwhitneyu(mean_expression_comp3_gbm, mean_expression_comp3_pbmc, alternative='two-sided')
u_stat_comp4, p_val_comp4 = mannwhitneyu(mean_expression_comp4_gbm, mean_expression_comp4_pbmc, alternative='two-sided')
u_stat_mitofission, p_val_mitofission = mannwhitneyu(mean_expression_mitofission_gbm, mean_expression_mitofission_pbmc, alternative='two-sided')
u_stat_mitofussion, p_val_mitofussion = mannwhitneyu(mean_expression_mitofussion_gbm, mean_expression_mitofussion_pbmc, alternative='two-sided')

u_stat_cytochrome_c_release_apop, p_val_cytochrome_c_release_apop = mannwhitneyu(mean_expression_cytochrome_c_release_apop_gbm, mean_expression_cytochrome_c_release_apop_pbmc, alternative='two-sided')
u_stat_mitophagy, p_val_mitophagy = mannwhitneyu(mean_expression_mitophagy_gbm, mean_expression_mitophagy_pbmc, alternative='two-sided')

# Create DataFrame to store results
results_df = pd.DataFrame({
    'Gene Set': ['ROS', 'Oxidative Stress', 'TGF-beta', 'Oxidative Phosphorylation', 'Ox Phos Rate Limiting Enzymes', 'Glycolysis and Gluconeogenesis', 'Glycolysis Rate Limiting Enzymes', 'Mitochondrial Complex I', 'Mitochondrial Complex II', 'Mitochondrial Complex III', 'Mitochondrial Complex IV', 'Mitochondrial Fission', 'Mitochondrial Fussion', 'Cytochrome C Release Induced Apoptosis', 'Mitophagy'],
    'U Statistic': [u_stat_ros, u_stat_oxidative_stress, u_stat_tgfbeta, u_stat_OXPHOS, u_stat_oxphos_ratelim, u_stat_GLYCOL_GLUCONEOGEN, u_stat_Glycolysis_ratelim, u_stat_comp1, u_stat_comp2, u_stat_comp3, u_stat_comp4, u_stat_mitofission, u_stat_mitofussion, u_stat_cytochrome_c_release_apop, u_stat_mitophagy],
    'p-value': [p_val_ros, p_val_oxidative_stress, p_val_tgfbeta, p_val_OXPHOS, p_val_oxphos_ratelim, p_val_GLYCOL_GLUCONEOGEN,  p_val_Glycolysis_ratelim, p_val_comp1, p_val_comp2, p_val_comp3, p_val_comp4, p_val_mitofission, p_val_mitofussion, p_val_cytochrome_c_release_apop, p_val_mitophagy],
    'TiCD56bright_glioblastoma Mean Expression': [mean_expression_ros_gbm.mean(), mean_expression_oxidative_stress_gbm.mean(), mean_expression_tgfbeta_gbm.mean(), mean_expression_OXPHOS_gbm.mean(), mean_expression_oxphos_ratelim_gbm.mean(), mean_expression_GLYCOL_GLUCONEOGEN_gbm.mean(), mean_expression_Glycolysis_ratelim_gbm.mean(), mean_expression_comp1_gbm.mean(), mean_expression_comp2_gbm.mean(), mean_expression_comp3_gbm.mean(), mean_expression_comp4_gbm.mean(), mean_expression_mitofission_gbm.mean(), mean_expression_mitofussion_gbm.mean(), mean_expression_cytochrome_c_release_apop_gbm.mean(), mean_expression_mitophagy_gbm.mean()],
    'TiCD56dim_glioblastoma Mean Expression': [mean_expression_ros_pbmc.mean(), mean_expression_oxidative_stress_pbmc.mean(), mean_expression_tgfbeta_pbmc.mean(), mean_expression_OXPHOS_pbmc.mean(), mean_expression_oxphos_ratelim_pbmc.mean(), mean_expression_GLYCOL_GLUCONEOGEN_pbmc.mean(), mean_expression_Glycolysis_ratelim_pbmc.mean(), mean_expression_comp1_pbmc.mean(), mean_expression_comp2_pbmc.mean(), mean_expression_comp3_pbmc.mean(), mean_expression_comp4_pbmc.mean(), mean_expression_mitofission_pbmc.mean(), mean_expression_mitofussion_pbmc.mean(), mean_expression_cytochrome_c_release_apop_pbmc.mean(), mean_expression_mitophagy_pbmc.mean()]
})

# Determine significance and which group to annotate
def get_annotation(row):
    if row['p-value'] < 0.05:
        if row['TiCD56bright_glioblastoma Mean Expression'] > row['TiCD56dim_glioblastoma Mean Expression']:
            return '*'
        else:
            return 'non-sig'
    else:
        return 'non-sig'

results_df['Significance'] = results_df.apply(get_annotation, axis=1)

# Prepare data for heatmap plotting
heatmap_data = pd.DataFrame({
    '(ROS)': [mean_expression_ros_gbm.mean(), mean_expression_ros_pbmc.mean()],
    '(Oxidative Stress)': [mean_expression_oxidative_stress_gbm.mean(), mean_expression_oxidative_stress_pbmc.mean()],
    '(TGF-beta)': [mean_expression_tgfbeta_gbm.mean(), mean_expression_tgfbeta_pbmc.mean()],
    '(Oxidative Phosphorylation)': [mean_expression_OXPHOS_gbm.mean(), mean_expression_OXPHOS_pbmc.mean()],
    '(Ox Phos Rate Limiting Enzymes)': [mean_expression_oxphos_ratelim_gbm.mean(), mean_expression_oxphos_ratelim_pbmc.mean()],
    '(Glycolysis and Gluconeogenesis)': [mean_expression_GLYCOL_GLUCONEOGEN_gbm.mean(), mean_expression_GLYCOL_GLUCONEOGEN_pbmc.mean()],
    '(Glycolysis Rate Limiting Enzymes)': [mean_expression_Glycolysis_ratelim_gbm.mean(), mean_expression_Glycolysis_ratelim_pbmc.mean()],
    '(Mitochondrial Complex I)': [mean_expression_comp1_gbm.mean(), mean_expression_comp1_pbmc.mean()],
    '(Mitochondrial Complex II)': [mean_expression_comp2_gbm.mean(), mean_expression_comp2_pbmc.mean()],
    '(Mitochondrial Complex III)': [mean_expression_comp3_gbm.mean(), mean_expression_comp3_pbmc.mean()],
    '(Mitochondrial Complex IV)': [mean_expression_comp4_gbm.mean(), mean_expression_comp4_pbmc.mean()],
    '(Mitochondrial Fission)': [mean_expression_mitofission_gbm.mean(), mean_expression_mitofission_pbmc.mean()],
    '(Mitochondrial Fussion)': [mean_expression_mitofussion_gbm.mean(), mean_expression_mitofussion_pbmc.mean()],
    '(Cytochrome C Release Induced Apoptosis)': [mean_expression_cytochrome_c_release_apop_gbm.mean(), mean_expression_cytochrome_c_release_apop_pbmc.mean()],
    '(Mitophagy)': [mean_expression_mitophagy_gbm.mean(), mean_expression_mitophagy_pbmc.mean()],
}, index=['TiCD56bright_glioblastoma', 'TiCD56dim_glioblastoma'])

# Add significance markers to heatmap annotations
def significance_annotation(gene_set):
    subset = results_df[results_df['Gene Set'] == gene_set]
    if not subset.empty:
        significance = subset['Significance'].values[0]
        return ['*', 'non-sig'] if significance == '*' else ['non-sig', 'non-sig']
    else:
        return ['non-sig', 'non-sig']

heatmap_annotations = pd.DataFrame({
    'ROS': significance_annotation('ROS'),
    'Oxidative Stress': significance_annotation('Oxidative Stress'),
    'TGF-beta': significance_annotation('TGF-beta'),
    'Oxidative Phosphorylation': significance_annotation('Oxidative Phosphorylation'),
    'Ox Phos Rate Limiting Enzymes': significance_annotation('Ox Phos Rate Limiting Enzymes'),
    'Glycolysis and Gluconeogenesis': significance_annotation('Glycolysis and Gluconeogenesis'),
    'Glycolysis Rate Limiting Enzymes': significance_annotation('Glycolysis Rate Limiting Enzymes'),
    'Mitochondrial Complex I': significance_annotation('Mitochondrial Complex I'),
    'Mitochondrial Complex II': significance_annotation('Mitochondrial Complex II'),
    'Mitochondrial Complex III': significance_annotation('Mitochondrial Complex III'),
    'Mitochondrial Complex IV': significance_annotation('Mitochondrial Complex IV'),
    'Mitochondrial Fission': significance_annotation('Mitochondrial Fission'),
    'Mitochondrial Fussion': significance_annotation('Mitochondrial Fussion'),
    'Cytochrome C Release Induced Apoptosis': significance_annotation('Cytochrome C Release Induced Apoptosis'),
    'Mitophagy': significance_annotation('Mitophagy'),
}, index=['Glioblastoma', 'TiCD56dim_glioblastoma'])

# Plot the heatmap with expression levels and significance markers
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data, cmap='RdBu_r', annot=heatmap_annotations, fmt='', cbar_kws={'label': 'Cumulative Expression Level'})

# Annotate the heatmap with significance markers
plt.title('Heatmap for Differential Gene Family Expression: TINKS and Metabolic Signatures')
plt.xlabel('Gene Set')
plt.ylabel('Subsets')

output_path = r"C:\Users\aroche6\abc\heatmap_gene_expression.svg"
plt.savefig(output_path, format="svg")
plt.close()  # Close the plot to free memory
# Save the plot as an SVG file

# Save the cumulative expression values and significance to an Excel file
excel_file_path = r"C:\Users\aroche6\abc\TRIAL20.xlsx"
with pd.ExcelWriter(excel_file_path) as writer:
    heatmap_data.to_excel(writer, sheet_name='Cumulative Expression')
    results_df.to_excel(writer, sheet_name='Significance')
