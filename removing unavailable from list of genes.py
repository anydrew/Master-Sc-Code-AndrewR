# -*- coding: utf-8 -*-
"""
Created on Sun Aug 11 14:53:39 2024

@author: AROCHE6
"""

# Define the lists
list1 = ["ATG12", "ATG5", "CSNK2A1", "CSNK2A2", "CSNK2B", "FUNDC1", "MAP1LC3A", "MAP1LC3B", "MFN1", "MFN2", "MTERF3", "PGAM5", "PINK1", "PRKN", "RPS27A", "SQSTM1", "SRC", "TOMM20", "TOMM22", "TOMM40", "TOMM5", "TOMM6", "TOMM7", "TOMM70", "UBA52", "UBB", "UBC", "ULK1", "VDAC1"]
list2 = ['PRKN', 'RPS27A', 'TOMM6', 'TOMM70', 'UBA52']
# Remove items in list2 from list1
list1 = [item for item in list1 if item not in list2]

# Print the result
print(list1)