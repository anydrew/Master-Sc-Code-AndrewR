# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 00:37:24 2025

@author: AROCHE6
"""
# Define color mapping for genes based on significance and log2_fold_change
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text  # Correct import

# Load the Excel file
file_path = r"C:\Users\aroche6\abc\LOG DFS PBMC VS TUMOR.xlsx"
deg_results = pd.read_excel(file_path)

# Ensure the data has the necessary columns
required_columns = ['gene', 'log2_fold_change', 'pvals', 'pvals_adj']
if not all(col in deg_results.columns for col in required_columns):
    raise ValueError(f"Input file must contain the following columns: {required_columns}")

# Calculate -log10 p-values
deg_results['-log10_pvals_adj'] = -np.log10(deg_results['pvals_adj'])

# Define the TGFB genes of interest
TGFB_genes = ['TGFBR1', 'SLC25A12', 'PRKAB2', 'VHL', 'SMPD3', 'MAPK9', 'MAPK3', 'TGFBR3', 'PRKCQ', 'EPAS1', 'PGD', 'SLC25A1', 'FOS', 'NUDT1', 'TNFRSF1A', 'ZFP36', 'SRXN1', 'CYB5B', 'ANXA1', 'ABCC1', 'KLF2', 'NDUFS2', 'PLK1', 'HDAC2', 'ERO1A', 'PPIA', 'PRDX6', 'PARK7', 'TNF', 'POLG', 'ERCC8', 'PRKAG2', 'EIF2S1', 'MSRB1', 'PRDX3', 'LONP1', 'TXN', 'FER', 'FOXO1', 'CDK1', 'MTHFD1', 'SOD1', 'TNFAIP3', 'PLCG2', 'GSS', 'NDUFS1', 'MAP1LC3A', 'MAPKAP1', 'BNIP3', 'HSF1', 'MCTP1', 'OLR1', 'PLA2G4A', 'POLB', 'SLC2A4', 'SLC25A22', 'MAPK13', 'P4HB', 'MET', 'NRG1', 'DDR2', 'FANCC', 'CAT', 'STK11', 'UQCRC1', 'FOSL1', 'BRF2', 'CYGB', 'HM13', 'ETFDH', 'APP', 'NDUFA9', 'CRYAB', 'GPX4', 'GCLC', 'IL18RAP', 'FOXO4', 'KDM6B', 'MITF', 'HDAC6', 'ARNT', 'RB1', 'NFE2L1', 'ATP7A', 'CD38', 'GPX2', 'IL18BP', 'G6PD', 'GCH1', 'TXNIP', 'GPX7', 'DUOX1', 'LRRK2', 'IL1A', 'IL6', 'MICB', 'OGG1', 'HMOX1', 'MAP2K4', 'PLCG1', 'ERCC2', 'ABL1', 'PRKAA2', 'ROS1', 'CBX8', 'NR4A2', 'PRDX1', 'PRDX2', 'AXL', 'DHCR24', 'KAT2B', 'TGFBR2', 'TGFB1', 'CASP3', 'SLC25A4', 'PLA2G5', 'SLC25A25', 'AREG', 'PDXK', 'TEK', 'CHUK', 'GGT7', 'DAPK1', 'MAPK7', 'JUN', 'FYN', 'CAMKK2', 'PGM3', 'SLC25A5', 'ALAD', 'PML', 'HSPA1A', 'BTK', 'HYAL2', 'PLCB1', 'AIFM2', 'MTHFR', 'PIK3R5', 'AKT1', 'FBLN5', 'SOD2', 'ETV5', 'RNF4', 'RAP1B', 'SESN1', 'SIRT3', 'MSRA', 'SLC25A39', 'TLR2', 'NRG2', 'ERCC3', 'EZH2', 'MDM2', 'SLC2A1', 'KEAP1', 'BAK1', 'EGFR', 'HYAL1', 'ERMP1', 'FOXP1', 'NLRP3', 'RCAN1', 'RAD23A', 'SETX', 'CAPN2', 'BCL2', 'GSR', 'LIAS', 'TGFB2', 'MPO', 'ATF2', 'PRKAG1', 'TFAM', 'ATF4', 'PA2G4', 'TNFRSF1B', 'PKM', 'NOX4', 'SIRT1', 'CPEB2', 'ANKZF1', 'MLH1', 'ALOX5', 'TXNDC5', 'XPC', 'CCS', 'GSTP1', 'SLC25A32', 'DHFR', 'OTUD1', 'GJB2', 'UQCRC2', 'ATP13A2', 'NADK', 'APOD', 'ATOX1', 'ERCC6L2', 'MTHFD2', 'HSPA1B', 'MAP3K5', 'PRDX4', 'SLC25A29', 'HIF1A', 'SLC25A19', 'PRDX5', 'SIRT5', 'FKBP1B', 'SLC25A6', 'SLC25A38', 'CD36', 'ERN1', 'PLCB3', 'COL1A1', 'ARL6IP5', 'SLC11A2', 'SLC25A15', 'PPARGC1A', 'BMP7', 'ATP2A2', 'GLRX2', 'JAK2', 'FOXO3', 'CRK', 'RBM38', 'EDN1', 'APOE', 'PXDN', 'NQO1', 'RAB20', 'HMOX2', 'ERCC1', 'FUT8', 'GPX8', 'AIF1', 'ERCC6', 'GPX3', 'ALS2', 'BECN1', 'CHCHD2', 'MMP7', 'ALDH3B1', 'FANCD2', 'PRKAB1', 'TP53', 'NFE2L2', 'PGK1', 'AIFM1', 'ECT2', 'OXSR1', 'HTRA2', 'MGST1', 'NDUFS8', 'ADAM9', 'SLC2A3', 'MAPK8', 'TBXAS1', 'IDH1', 'MAPT', 'RUNX2', 'CYP1B1', 'FXN', 'MGAT3', 'PRKAA1', 'TLR4', 'SOD3', 'PIK3R1', 'MAPK1', 'TGFB3', 'BANF1', 'UQCRFS1', 'IPCEF1', 'NCF1', 'AGAP3', 'GCLM', 'ABCD1', 'PNP', 'ATM']

# [
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
#     'MLH1', 'MMP7', 'MPO', 'MSRA', 'MSRB1', 'MTHFD1', 'MTHFD2', 'MTHFR', 'NADK', 'NCF1', 
#     'NDUFA9', 'NDUFS1', 'NDUFS2', 'NDUFS8', 'NFE2L1', 'NFE2L2', 'NLRP3', 'NOX4', 'NQO1', 
#     'NRG1', 'NRG2', 'NR4A2', 'NUDT1', 'OGG1', 'OLR1', 'OTUD1', 'OXSR1', 'P4HB', 'PA2G4', 'PARK7', 
#     'PDXK', 'PGD', 'PGK1', 'PGM3', 'PIK3R1', 'PIK3R5', 'PKM', 'PLA2G4A', 'PLA2G5', 
#     'PLCB1', 'PLCB3', 'PLCG1', 'PLCG2', 'PLK1', 'PML', 'PNP', 'POLB', 'POLG', 
#     'PPARGC1A', 'PPIA', 'PRDX1', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6', 'PRKAA1', 
#     'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKCQ', 
#     'PXDN', 'RAB20', 'RAD23A', 'RAP1B', 'RB1', 'RBM38', 'RCAN1', 'RNF4', 'ROS1', 
#     'RUNX2', 'SIRT1', 'SIRT3', 'SIRT5', 'SLC11A2', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
#     'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A1', 
#     'SLC25A12', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
#     'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC25A4', 'SLC25A5', 'SLC25A6', 
#      'SLC25A12', 'SLC25A15', 'SLC25A19', 'SETX', 'SESN1',
#     'SLC25A22', 'SLC25A25', 'SLC25A29', 'SLC25A32', 'SLC25A38', 'SLC25A39', 
#     'SOD1', 'SOD2', 'SOD3', 'SRXN1', 'STK11', 'TBXAS1', 'TEK', 'TFAM', 'TGFB1', 
#     'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TGFBR3', 'TLR2', 'TLR4', 
#     'TNF', 'TNFAIP3', 'TNFRSF1A', 'TNFRSF1B', 'TP53', 'TXN', 'TXNDC5', 'TXNIP', 
#     'UQCRC1', 'UQCRC2', 'UQCRFS1', 'VHL', 'XPC', 'ZFP36'
# ]
# Filter the results to include only the TGFB genes
tgfb_results = deg_results[deg_results['gene'].isin(TGFB_genes)]

# TGF_SIG_RANK_GENES = ["JUNB", "ID2", "JUND", "JUN", "PPP1R15A", "MAPK1", "FKBP1A", "CTNNB1", "CD9", "IFNG", "SKI", "TFDP1", "CUL1", "EP300", "SMURF2", "SPTBN1", "ID3", "MYC", "MAPK8"]

def gene_color(row):
    if row['pvals_adj'] < 0.05:
        if row['log2_fold_change'] <= -1:
            return 'blue'
        elif row['log2_fold_change'] > 1:
            return 'red'
    return 'grey'

# Apply the color mapping to the data
tgfb_results['color'] = tgfb_results.apply(gene_color, axis=1)

# Filter significant genes for saving to Excel
significant_genes = tgfb_results[(tgfb_results['color'] != 'grey')]

# # Save significant genes to an Excel file
output_file_path = r"C:\Users\aroche6\abc\REAL_ADJpvalsignificant_OXSTRESS_genes.xlsx"
significant_genes.to_excel(output_file_path, index=False)

# Plot the volcano plot for TGFB genes using p-values with specified colors
plt.figure(figsize=(15, 10))
sns.scatterplot(
    x='log2_fold_change', 
    y='-log10_pvals_adj', 
    data=tgfb_results, 
    hue='color', 
    palette={'blue': 'blue', 'red': 'red', 'grey': 'grey'}, 
    legend=False,
    edgecolor=None
)

# Add labels to significant genes
texts = []
for i, row in tgfb_results[(tgfb_results['color'] != 'grey')].iterrows():
    texts.append(plt.text(row['log2_fold_change'], row['-log10_pvals_adj'], row['gene'], fontsize=8))

# Adjust text to prevent overlap
adjust_text(texts)

# Add plot annotations and styling
plt.title('Volcano Plot of Differential Gene Expression for ROS Genes: Glioblastoma vs PBMC')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-Value')
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')
plt.axvline(x=1, color='black', linestyle='--')
plt.axvline(x=-1, color='black', linestyle='--')

# Set x-axis limit to [-10, 10]
plt.xlim(-5, 10)

plt.grid(True)
plt.show()

# ros_genes = [
#     "ABCC1", "ATOX1", "CAT", "CDKN2D", "EGLN2", "ERCC2", "FES", "FTL", "G6PD", "GCLC", 
#     "GCLM", "GLRX", "GLRX2", "GPX3", "GPX4", "GSR", "HHEX", "HMOX2", "IPCEF1", "JUNB", 
#     "LAMTOR5", "MBP", "MGST1", "MSRA", "NDUFA6", "NDUFB4", "NDUFS2", "NQO1", "OXSR1", 
#     "PDLIM1", "PFKP", "PRDX1", "PRDX2", "PRDX4", "PRDX6", "PRNP", "SBNO2", "SCAF4", 
#     "SOD1", "SOD2", "STK25", "TXN", "TXNRD1", "TXNRD2"
# ]
# tgfbeta_genes = [
#     "CD9", "ATF3", "SIN3A", "MAPK11", "DCN", "ENG", "MAPK9", "HFE", "CDK8", "PPP1R15A", 
#     "BMP8B", "FNTA", "ATP1B2", "CITED2", "IFNGR2", "LTBP2", "SMAD1", "BMPR1B", "NCOR1", 
#     "TAB1", "E2F5", "LTBP3", "PPM1A", "XIAP", "SMAD3", "MAPK12", "ID1", "TRIM33", 
#     "SERPINE1", "GDF15", "BMP2", "THBS1", "CREBBP", "CDK9", "CTSC", "IFNG", "BCAR3", 
#     "RGMB", "EMP3", "GDF11", "FURIN", "WWTR1", "SKIL", "MAPK14", "RGMA", "PPP2CA", 
#     "EP300", "FKBP1A", "VPS51", "ACVR2A", "ACVR1B", "E2F4", "CDKN1C", "SLC20A1", 
#     "BMPR1A", "LTBP1", "TGFBR2", "RHOA", "ACVR1", "SMAD4", "INHBA", "ID3", "PPP1CA", 
#     "MMP9", "CALD1", "JUNB", "SMURF1", "SP1", "TGFB3", "ZFYVE9", "NRAS", "NBL1", 
#     "MAPK13", "NRROS", "SNIP1", "BAMBI", "TFRC", "MAPK1", "MYC", "NCOR2", "GREM1", 
#     "VPS54", "SMAD9", "UBE2D3", "HDAC1", "SMURF2", "ITGB8", "ID2", "APC", "ACVRL1", 
#     "TIMP1", "PPP2CB", "CTNNB1", "BMP6", "RRAS", "BMPR2", "SPTBN1", "KLF10", "AMH", 
#     "ATF2", "FOSL1", "MAPK8", "HRAS", "SKP1", "JUN", "TFDP1", "TGFB1", "RAB31", 
#     "ATP1B3", "ITGAV", "RNF111", "LRRC32", "SMAD6", "ILK", "MAPK3", "TFR2", "PPP2R1A", 
#     "SMAD2", "BMP7", "TGIF2", "MAP3K7", "VPS53", "ZFYVE16", "TGFB2", "ARID4B", 
#     "ITGAE", "JUND", "DCP1B", "TJP1", "MMP2", "ACVR1C", "ROCK1", "MAPK10", "CDKN2B", 
#     "MAP3K7CL", "FBN1", "BMP1", "NEO1", "HDAC2", "TGIF1", "PMEPA1", "CDKN1B", "RBX1", 
#     "SKI", "FST", "ATF1", "LRRC2", "CUL1", "SMAD7", "ASCL2", "ID4", "TGFBR1", "ACVR2B", 
#     "CD40", "PPP2R1B", "SMAD5", "HIPK2", "THSD4"
# ]

# # Define Oxidative Stress gene list ros_genes, tgfbeta_genes, oxidative_stress_genes, OXPHOS_genes, oxphos_ratelim, GLYCOL_GLUCONEOGEN_genes, Glycolysis_ratelim

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
# OXPHOS_genes = ['ATP1B1', 'COX10', 'UQCRH', 'ATP6V1B2', 'NDUFA5', 'NDUFA4', 'ATP6V1C1', 'NDUFC2', 'CLDN3', 'COX6C', 'LHPP', 'UQCRFS1', 'NDUFB2', 'NDUFA10', 'NDUFA4L2', 'NDUFS5', 'CD4', 'NDUFA13', 'RERE', 'UQCRQ', 'UQCRB', 'COX17', 'NDUFB7', 'COX11', 'NDUFS1', 'NDUFS4', 'ATP6V0B', 'COX7A2L', 'ATP6V0C', 'TPP1', 'NDUFB3', 'ATP6V1D', 'SURF1', 'CYC1', 'NDUFB5', 'NDUFV2', 'COX6A1', 'COX4I1', 'COX8A', 'NDUFB11', 'NDUFA12', 'COX7A2', 'NDUFA2', 'NDUFC1', 'PPA1', 'NDUFA9', 'NDUFB6', 'NDUFV1', 'NDUFB1', 'ATP6V0E2', 'COX7C', 'ATP6V1E1', 'COX5A', 'ATP6V1G1', 'COX5B', 'NDUFS6', 'ATP6V0E1', 'CPM', 'COX7B', 'ATP6V1A', 'CP', 'NDUFAB1', 'UQCRC2', 'NDUFA11', 'NDUFS8', 'NDUFS3', 'NDUFA6', 'ATP6V1E2', 'UQCR10', 'NDUFA3', 'COX15', 'NDUFB9', 'SDHA', 'NDUFA8', 'ATP6V0A2', 'SDHC', 'CPE', 'ATP6V1H', 'NDUFV3', 'ATP6V0D1', 'NDUFS7', 'NDUFA1', 'NDUFB10', 'TCIRG1', 'NDUFS2', 'UQCRC1', 'ATP6V1F', 'PPA2', 'NDUFB4', 'ATP6AP1', 'ATP6V0A1', 'SDHB', 'COX6B1']
# oxphos_ratelim = ["COX10", "COX11", "COX14", "COX15", "COX16", "COX17", "COX18", "COX19", "COX20", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", "COX7A2L", "COX7B", "COX7C", "COX8A", "SURF1", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "HIGD1A", "HIGD2A"]

# GLYCOL_GLUCONEOGEN_genes = ['ADORA2B', 'AGL', 'AGRN', 'AK3', 'AK4', 'AKR1A1', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALG1', 'ANG', 'ANGPTL4', 'ANKZF1', 'ARPP19', 'AURKA', 'B3GALT6', 'B3GAT1', 'B3GAT3', 'B4GALT1', 'B4GALT2', 'B4GALT4', 'B4GALT7', 'BIK', 'BPNT1', 'CACNA1H', 'CAPN5', 'CASP6', 'CD44', 'CDK1', 'CENPA', 'CHPF', 'CHPF2', 'CHST1', 'CHST12', 'CHST2', 'CITED2', 'CLDN3', 'CLN6', 'COG2', 'COL5A1', 'COPB2', 'CTH', 'CXCR4', 'CYB5A', 'DCN', 'DDIT4', 'DEPDC1', 'DLD', 'DSC2', 'ECD', 'EGFR', 'EGLN3', 'ELF3', 'ENO1', 'ENO2', 'ERO1A', 'EXT1', 'EXT2', 'FAM162A', 'FKBP4', 'FUT8', 'G6PD', 'GALE', 'GALK1', 'GALK2', 'GCLC', 'GFPT1', 'GLCE', 'GLRX', 'GMPPA', 'GMPPB', 'GNE', 'GNPDA1', 'GOT1', 'GOT2', 'GPC1', 'GPC3', 'GPC4', 'GUSB', 'GYS1', 'HAX1', 'HDLBP', 'HK2', 'HMMR', 'HOMER1', 'HS2ST1', 'HS6ST2', 'HSPA5', 'IDH1', 'IDUA', 'IER3', 'IGFBP3', 'IL13RA1', 'IRS2', 'ISG20', 'KDELR3', 'KIF20A', 'KIF2A', 'LDHA', 'LHPP', 'MDH1', 'MDH2', 'ME2', 'MED24', 'MET', 'MIF', 'MPI', 'MXI1', 'NANP', 'NASP', 'NDUFV3', 'NOL3', 'NSDHL', 'NT5E', 'P4HA1', 'P4HA2', 'PAM', 'PAXIP1', 'PC', 'PDK3', 'PFKM', 'PFKL', 'PFKP', 'PGAM1', 'PGAM2', 'PGK1', 'PGLS', 'PGM2', 'PHKA2', 'PKM', 'PKP2', 'PLOD1', 'PLOD2', 'PMM2', 'POLR3K', 'PPIA', 'PPP2CB', 'PRPS1', 'PSMC4', 'PYGB', 'PYGL', 'QSOX1', 'RBCK1', 'RPE', 'RRAGD', 'SAP30', 'SDC1', 'SDC2', 'SDC3', 'SDHC', 'SLC16A3', 'SLC25A13', 'SLC35A3', 'SLC37A4', 'SOD1', 'SPAG4', 'SRD5A3', 'STMN1', 'TALDO1', 'TFF3', 'TGFA', 'TGFBI', 'TPBG', 'TPI1', 'TPST1', 'TXN', 'UGP2', 'VCAN', 'VEGFA', 'VLDLR', 'XYLT2', 'ZNF292']
# Glycolysis_ratelim = ["SLC2A1", "HK1", "HK2", "HK3", "PFKM", "PFKL", "PFKP", "SLC16A4", "SLC16A1"]
# #ros_genes, tgfbeta_genes, oxidative_stress_genes, OXPHOS_genes, oxphos_ratelim, GLYCOL_GLUCONEOGEN_genes, Glycolysis_ratelim, comp1, comp2, comp3, comp4, mitofission_genes, mitofussion_genes, cytochrome_c_release_apop_genes, mitophagy_genes

# # comp1 = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6']
# # comp2 = ['SDHA', 'SDHB', 'SDHC']
# # comp3 = ['UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10']
# # comp4 = ['COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']
# mitocomps = ['NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFC1', 'NDUFC2', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'SDHA', 'SDHB', 'SDHC', 'UQCRB', 'UQCRQ', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRFS1', 'UQCRH', 'UQCR10', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7B', 'COX7C', 'COX8A']
# mitofission_genes = ['VPS35', 'SLC25A46', 'MTFR1', 'COX10', 'TMEM135', 'DCN', 'AURKA', 'MTFR2', 'MUL1', 'MTFP1', 'AP3B1', 'DDHD1', 'MIEF1', 'MCU', 'FIS1', 'MFF', 'PINK1', 'MYO19', 'DDHD2', 'OPA1', 'RALBP1', 'RALA', 'INF2', 'BAX', 'CCAR2', 'UCP2', 'LRRK2', 'STAT2', 'SPIRE1', 'ATG3', 'PGAM5', 'LPIN1', 'MAPT', 'GDAP1', 'MIEF2', 'BNIP3', 'MTFR1L', 'DNM1L', 'PPARG']
# mitofussion_genes = ['MFN1', 'AFG3L2', 'MFN2', 'ADCK1', 'PID1', 'MUL1', 'STOML2', 'MIEF1', 'THG1L', 'ZDHHC6', 'FIS1', 'MFF', 'MTCH2', 'HUWE1', 'OPA1', 'BAK1', 'USP30', 'BAX', 'MCL1', 'GDAP1', 'VAT1', 'MIEF2', 'BNIP3', 'PLD6', 'CHCHD3', 'OMA1', 'PARL', 'TFRC', 'BCL2A1']

# cytochrome_c_release_apop_genes = ['AVEN', 'BBC3', 'BCL2L1', 'BMF', 'CASP3', 'PMAIP1', 'XIAP', 'PINK1', 'CARD8', 'MLLT11', 'NOL3', 'CASP9', 'BIRC2', 'AIFM1', 'CYCS', 'APIP', 'DNM1L', 'PARL', 'BCL2A1', 'BIRC3', 'UACA', 'ENDOG', 'BIK', 'DFFA', 'DIABLO', 'CIDEB', 'SOD2', 'FIS1', 'MFF', 'BAK1', 'BAX', 'GGCT', 'MMP9', 'FAM162A', 'BNIP3', 'TNFSF10', 'MAPK1', 'PPIF', 'APAF1', 'MOAP1', 'DFFB', 'TP53', 'BCL2', 'CASP8', 'IGF1', 'PDCD5', 'LMNA', 'CASP6', 'FXN', 'PSMD10', 'PYCARD', 'MAPK3', 'TIMM50', 'GHITM', 'BAD', 'BOK', 'PLAUR', 'PRELID1', 'SFN', 'OPA1', 'BCL2L2', 'MCL1', 'BID', 'CLU', 'AKT1', 'IFI6', 'CASP7', 'BCL2L11', 'TRIAP1']
# mitophagy_genes = ['ATG12', 'ATG5', 'CSNK2A1', 'CSNK2A2', 'CSNK2B', 'FUNDC1', 'MAP1LC3A', 'MAP1LC3B', 'MFN1', 'MFN2', 'MTERF3', 'PGAM5', 'PINK1', 'SQSTM1', 'SRC', 'TOMM20', 'TOMM22', 'TOMM40', 'TOMM5', 'TOMM7', 'UBB', 'UBC', 'ULK1', 'VDAC1']