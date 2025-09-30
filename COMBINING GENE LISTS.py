# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 15:32:43 2024

@author: AROCHE6
"""

def combine_and_remove_duplicates(*lists):
    combined_list = []
    for lst in lists:
        combined_list.extend(lst)  # Add each list to the combined list
    unique_list = list(set(combined_list))  # Remove duplicates by converting to a set and back to a list
    return unique_list
# Example usage
# list1 = ["ACVR1", "ACVR1C", "ACVR2A", "ACVR2B", "ACVRL1", "AMH", "AMHR2", "BMP2", 
#     "BMP4", "BMP5", "BMP6", "BMP7", "BMP8A", "BMP8B", "BMPR1A", "BMPR1B", 
#     "BMPR2", "CDKN2B", "CHRD", "COMP", "CREBBP", "CUL1", "DCN", "E2F4", 
#     "E2F5", "EP300", "FST", "GDF5", "GDF6", "GDF7", "ID1", "ID2", "ID3", 
#     "ID4", "IFNG", "INHBA", "INHBB", "INHBC", "INHBE", "LEFTY1", "LEFTY2", 
#     "LTBP1", "MAPK1", "MAPK3", "MYC", "NODAL", "NOG", "PITX2", "PPP2CA", 
#     "PPP2CB", "PPP2R1A", "PPP2R1B", "RBL1", "RBL2", "RBX1", "RHOA", "ROCK1", 
#     "ROCK2", "RPS6KB1", "RPS6KB2", "SKP1", "SKP1P2", "SMAD1", "SMAD2", 
#     "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9", "SMURF1", "SMURF2", 
#     "SP1", "TFDP1", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "THBS1", 
#     "THBS2", "THBS3", "THBS4", "TNF", "ZFYVE16", "ZFYVE9"]
# list2 = ["ACVR1B", "ACVR2A", "ARRB2", "ASCL2", "ATF1", "ATF3", "ATF4", "ATP1B3",
#     "ATP6V0B", "BMP7", "BMP8B", "BMPR1A", "BMPR2", "CALM1", "CBLB", "CD40",
#     "CDK8", "CDKN1B", "CREB1", "CREBBP", "CTSC", "CUL1", "E2F4", "E2F5",
#     "EP300", "ID1", "ID2", "ID3", "IFNG", "MAPK1", "MAPK3", "MYC", "NBL1",
#     "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B", "RBX1", "RHOA", "ROCK1",
#     "RPS6KB1", "RPS6KB2", "SKP1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", 
#     "SMAD7", "SMURF1", "SMURF2", "SP1", "TFDP1", "TGFB1", "TGFBR1", 
#     "TGFBR2", "TGIF1", "TGIF2", "THBS1", "TNF", "ZFYVE16"]
# list3 = ['CHRD', 'NOG', 'NBL1', 'MICOS10-NBL1', 'GREM1', 'GREM2', 'THBS1', 'DCN', 'FMOD', 'LEFTY1', 'LEFTY2', 'FST', 'BMP2', 'BMP4', 'BMP6', 'INHBB', 'TF', 'BMP5', 'BMP7', 'BMP8B', 'BMP8A', 'GDF5', 'GDF6', 'GDF7', 'AMH', 'IFNG', 'TNF', 'THSD4', 'FBN1', 'LTBP1', 'TGFB1', 'TGFB2', 'TGFB3', 'INHBA', 'INHBC', 'INHBE', 'NODAL', 'BMPR1A', 'BMPR1B', 'ACVR1', 'BMPR2', 'ACVR2A', 'HJV', 'NEO1', 'HFE', 'TFR2', 'TFRC', 'RGMA', 'RGMB', 'AMHR2', 'TGFBR1', 'TGFBR2', 'EMP3', 'LRRC32', 'NRROS', 'ACVR1B', 'ACVR2B', 'ACVR1C', 'IGSF1', 'BAMBI', 'SMAD1', 'SMAD5', 'SMAD9', 'SMAD2', 'SMAD3', 'SMAD4', 'SMAD6', 'SMAD7', 'SMURF1', 'SMURF2', 'ZFYVE9', 'ZFYVE16', 'RBX1', 'CUL1', 'SKP1', 'MAPK1', 'MAPK3', 'RHOA', 'ROCK1', 'PPP2R1B', 'PPP2R1A', 'PPP2CA', 'PPP2CB', 'RPS6KB1', 'RPS6KB2', 'HAMP', 'ID1', 'ID2', 'ID3', 'ID4', 'RBL1', 'E2F4', 'E2F5', 'TFDP1', 'CREBBP', 'EP300', 'SP1', 'SKIL', 'SKI', 'SIN3A', 'HDAC1', 'HDAC2', 'NCOR1', 'TGIF1', 'TGIF2', 'MYC', 'CDKN2B', 'PITX2']
# list4 = ["ID2", "SKP1", "CTSC", "TGFB1", "RBX1", "SMAD2", "SMAD3", "SMAD7", "SMURF1", "SMURF2", "LTBP1", "LTBP2", "LTBP3", "TGFBR1", "TGFBR2", "TGFB2", "TGFB3", "LRRC32", "MYC", "BMPR1", "BMPR2", "BMP2", "BMP7", "BMP4", "TGIF1", "TGIF2", "HDAC1", "SKI", "SKIL", "RNF111", "ARK2N", "ARK2C", "MMP2", "MMP3", "MMP9", "TIMP1", "ATF3", "JUNB", "FOXP3", "CDH1", "ITGAE", "CD9", "ITGAV", "ITGB3", "ITGB8", "ILK", "VPS51", "VPS52", "VPS53", "VPS54"]
# list5 = ["ACVR1B", "ACVR1C", "ACVR2A", "ACVR2B", "AMH", "ASCL2", "ATF1", "ATF3", "ATP1B2", "ATP1B3", "BAMBI", "BMP2", "BMP7", "BMP8B", "BMPR1A", "BMPR1B", "BMPR2", "CALD1", "CD40", "CDK8", "CDKN1B", "CDKN2B", "CREBBP", "CTSC", "CUL1", "EIPR1", "E2F4", "E2F5", "EP300", "FST", "ID1", "ID2", "ID3", "ID4", "INHBA", "LRRC2", "LTBP1", "NBL1", "RBX1", "SKP1", "SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7", "SMAD9", "SMURF1", "SMURF2", "SP1", "TFDP1", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGIF1", "TGIF2", "VPS53", "VPS54", "VPS52", "VPS51", "ZFYVE16", "ZFYVE9"]
list1 = [
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
    'MLH1', 'MMP7', 'MPO', 'MSRA', 'MSRB1', 'MTHFD1', 'MTHFD2', 'MTHFR', 'NADK', 'NCF1', 
    'NDUFA9', 'NDUFS1', 'NDUFS2', 'NDUFS8', 'NFE2L1', 'NFE2L2', 'NLRP3', 'NOX4', 'NQO1', 
    'NRG1', 'NRG2', 'NR4A2', 'NUDT1', 'OGG1', 'OLR1', 'OTUD1', 'OXSR1', 'P4HB', 'PA2G4', 'PARK7', 
    'PDXK', 'PGD', 'PGK1', 'PGM3', 'PIK3R1', 'PIK3R5', 'PKM', 'PLA2G4A', 'PLA2G5', 
    'PLCB1', 'PLCB3', 'PLCG1', 'PLCG2', 'PLK1', 'PML', 'PNP', 'POLB', 'POLG', 
    'PPARGC1A', 'PPIA', 'PRDX1', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6', 'PRKAA1', 
    'PRKAA2', 'PRKAB1', 'PRKAB2', 'PRKAG1', 'PRKAG2', 'PRKCQ', 
    'PXDN', 'RAB20', 'RAD23A', 'RAP1B', 'RB1', 'RBM38', 'RCAN1', 'RNF4', 'ROS1', 
    'RUNX2', 'SIRT1', 'SIRT3', 'SIRT5', 'SLC11A2', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
    'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC2A1', 'SLC2A3', 'SLC2A4', 'SLC25A1', 
    'SLC25A12', 'SLC25A15', 'SLC25A19', 'SLC25A22', 'SLC25A25', 'SLC25A29', 
    'SLC25A32', 'SLC25A38', 'SLC25A39', 'SLC25A4', 'SLC25A5', 'SLC25A6', 
     'SLC25A12', 'SLC25A15', 'SLC25A19', 'SETX', 'SESN1',
    'SLC25A22', 'SLC25A25', 'SLC25A29', 'SLC25A32', 'SLC25A38', 'SLC25A39', 
    'SOD1', 'SOD2', 'SOD3', 'SRXN1', 'STK11', 'TBXAS1', 'TEK', 'TFAM', 'TGFB1', 
    'TGFB2', 'TGFB3', 'TGFBR1', 'TGFBR2', 'TGFBR3', 'TLR2', 'TLR4', 
    'TNF', 'TNFAIP3', 'TNFRSF1A', 'TNFRSF1B', 'TP53', 'TXN', 'TXNDC5', 'TXNIP', 
    'UQCRC1', 'UQCRC2', 'UQCRFS1', 'VHL', 'XPC', 'ZFP36'
]
list2 = [
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
result = combine_and_remove_duplicates(list1, list2)
print(result)
