# Benjamin Adam Catching
# 2024-05-23
# NIH-NIAID-LVD-QVEU
# EV-D68 GxE Project

"""

"""

# Import packages
import pandas as pd 
import numpy as np

dNdS_list = []
gene_starts = np.array([69, 317, 564, 861, 1008, 1107, 1437, 1526, 1548, 1731, 2188])
genes = ['VP4', 'VP2', 'VP3', 'VP1', '2A', '2B', '2C', '3A', '3B', '3C', '3D']

P1_dNdS_list = list(sorted(glob.glob('../data/sequences/2023-02-06_passage_1/dNdS_per_gene/*.csv')))
P2_dNdS_list = list(sorted(glob.glob('../data/sequences/2023-02-27_passage_2/dNdS_per_gene/*.csv')))
P3_dNdS_list = list(sorted(glob.glob('../data/sequences/2023-03-14_passage_3/dNdS_per_gene/*.csv')))
P4_dNdS_list = list(sorted(glob.glob('../data/sequences/2023-04-06_passage_4/dNdS_per_gene/*.csv')))
P5_dNdS_list = list(sorted(glob.glob('../data/sequences/2023-05-30_passage_5/dNdS_per_gene/*.csv')))
P6_dNdS_list = list(sorted(glob.glob('../data/sequences/2023-06-18_passage_6/dNdS_per_gene/*.csv')))


for passage in range(0, 6, 1):
    print(passage)
    dNdS_files = [P1_dNdS_list, P2_dNdS_list, P3_dNdS_list, P4_dNdS_list, P5_dNdS_list, P6_dNdS_list][passage]
    count = 0
    for i, cell in enumerate(['SH-SY5Y', 'RD', 'A549']):
        for j, temperature in enumerate([33, 37]):
            for k, strain in enumerate(['US/MO/14-18947', 'Fermon']):
                for l, rep in enumerate('ABCDEFGH'):
                    temp_len = len(open(dNdS_files[count], 'r').readlines())
                    if temp_len != 0:
                        dNdS_temp_df = pd.read_csv(dNdS_files[count]).merge(dummy_df, on='#gene', how='right')
                        for kind in genes:
                            type_df = dNdS_temp_df[dNdS_temp_df['#gene'] == kind]
                            dNdS_list.append([passage+1, cell, temperature, strain, kind, rep, type_df['dn/ds'].iloc[0], count+1])
                    count += 1

dNdS_df = pd.DataFrame(dNdS_list)
dNdS_df.columns = ['passage', 'cell', 'temperature', 'strain', 'gene', 'replicate', 'dn/ds', 'number']

dNdS_df.to_csv