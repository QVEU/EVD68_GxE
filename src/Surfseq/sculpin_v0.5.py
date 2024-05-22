# Benjamin Adam Catching
# 2023-03-02
# NIH-NIAID-LVD
# Quantitative Virology and Evolution Unit

"""
Slimy sculpin (Cottus cognatus), a nocturnal fish, v0.5
From aligned, phred-score = 20 cutoff, Q20.txt files, call each consensus.
mutation in each replicate of a 8-replicate, 12-condition passaging experiment
Returns a single .csv file with information for each nucleotide-position.
"""

# +++ TO DO +++ #
# - Parse Sam files for more accurate codon information
# - Fix melted_df being 4X too large with NaN in 3 out of 4 nucleotide counts

# Import modules
import numpy as np
import pandas as pd
import glob
from Bio import Seq

# Find locations of the sequencing from the initial virus (US/MO/2014-18947=MO, 
# Fermon=Fermon)

fermon_consensus = open('data/sequences/full_genome/fermon.fa').readlines()[1]
mo_consensus = open('data/sequences/full_genome/18947.fa').readlines()[1][1:]

# Locations of Nanopore Q20 files
print('Reading in files')
new_P1_dir = 'data/sequences/2023-04-17_passage_1/trimmed_Q20/*.txt'
new_P2_dir = 'data/sequences/2023-05-01_passage_2/trimmed_Q20/*.txt'
new_P3_dir = 'data/sequences/2023-05-05_passage_3/trimmed_Q20/*.txt'
new_P4_dir = 'data/sequences/2023-05-08_passage_4/trimmed_Q20/*.txt'
directories = [new_P1_dir, new_P2_dir, new_P3_dir, new_P4_dir]

# Locations of Miseq Q20 files
"""new_P1_dir = 'data/sequences/2023-04-17_passage_1/trimmed_Q20/*.txt'
new_P2_dir = 'data/sequences/2023-05-01_passage_2/trimmed_Q20/*.txt'
new_P3_dir = 'data/sequences/2023-05-05_passage_3/trimmed_Q20/*.txt'
new_P4_dir = 'data/sequences/2023-05-08_passage_4/trimmed_Q20/*.txt'"""


# Iterate through Q20 directories
passage_files = [sorted(glob.glob(directory)) for directory in directories]

# List comprehension of the data
Q20_df =[[pd.read_csv(x, delimiter='\t') for x in y] for y in passage_files]

# Define conditions
strains = ['US/MO/14-18947', 'Fermon']
cell_types = ['SH-SY5Y', 'RD', 'A549']
temps = ['33', '37']
print('Converting to DataFrame')
for j, passage in enumerate(Q20_df):
    for i, replicate in enumerate(passage):

        batch = i // 8
        replicate_num = i % 8
        replicate_chr = chr(ord('@') + replicate_num+1)
        strain_num = batch % 2
        strain = strains[strain_num]
        cell_num = batch // 4
        cell_type = cell_types[cell_num]
        temperature_num = batch // 2 % 2
        temperature = temps[temperature_num]
        #print(j, i, strain, cell_num)
        
        replicate['passage'] = str(j+1)
        replicate['cell'] = cell_type
        replicate['strain'] = strain
        replicate['temperature'] = temperature
        replicate['replicate'] = replicate_chr
        replicate['nucleotide position'] = [x + 1 for x in list(replicate.index)]

# Concatanate Q20 file
Q20_df_list = pd.concat([pd.concat(Q20_df[i]) for i in range(len(Q20_df))])
Q20_df_list = Q20_df_list.rename(columns={'0': 'A', '0.1': 'C', '0.2': 'G', 
    '0.3': 'T', '0.4': 'N'})

Q20_df_list = Q20_df_list.drop(['1','8','2','4','3'], axis=1)
Q20_df_list['coverage'] = Q20_df_list['A'] + Q20_df_list['C'] +  Q20_df_list['G'] + Q20_df_list['T']

#Q20_df_list.to_csv('data/EV-D68_all_passage_data_v3.csv')
print('Melting DataFrame into nucleotide-position')
melted_df = pd.concat([Q20_df_list[[x, 'passage', 'cell', 'strain', 'temperature','replicate', 'nucleotide position', 'coverage']] for x in ['A', 'C', 'G', 'T']])
melted_df = pd.melt(melted_df, ['passage', 'cell', 'strain', 'temperature','replicate', 'nucleotide position', 'coverage'])
melted_df['percent'] = melted_df['value'] / melted_df['coverage']

# Need a list of Fermon protein start and stops
fermon_protein_pos = {"5' UTR":1, 'VP4':731, 'VP2':938, 'VP3':1682, 'VP1':2387,
                     '2A':3314, '2B':3755, '2C':4052,
                     '3A':5042, '3B':5309, '3C':5375, '3D':5924,
                     "3' UTR":7298}

# Need a list of US/MO/14-18947 protein start and stops
MO_protein_pos = {"5' UTR":1, 'VP4':697, 'VP2':904, 'VP3':1648, 'VP1':2353,
                     '2A':3280, '2B':3721, '2C':4018,
                     '3A':5008, '3B':5275, '3C':5341, '3D':5890,
                     "3' UTR":7264}
#melted_df.to_csv('data/EV-D68_all_passage_data_v5.csv')

print('Determining codon for each nucleotide')
temp_file = open('data/EV-D68_all_passage_data_v8.csv', 'w')
# Write header
temp_file.write('# Miseq data (passage 1-4)')
temp_file.write('passage,cell,strain,temperature,replicate,nucleotide position\
,coverage,nucleotide,nucleotide count,nucleotide percent,WT nucleotide,protein,protein \
position,polyprotein position,WT codon,mutant codon,WT aa,mutant aa,non-synonymous\n')

# Iterate through DataFrame (bad I know)
for i in range(len(melted_df)):
    # Define current row
    temp_row = melted_df.iloc[i+1]

    # Use strain from row to determine the protein and consensus data
    if temp_row['strain'] == 'Fermon':
        proteins = list(fermon_protein_pos.keys())
        start_list = list(fermon_protein_pos.values())
        consensus = fermon_consensus

    elif temp_row['strain'] == 'US/MO/14-18947':
        proteins = list(MO_protein_pos.keys())
        start_list = list(MO_protein_pos.values())
        consensus = mo_consensus

    # Split the input mutation into original nucleotide, position, and new nucleotide
    genome_pos = temp_row['nucleotide position']-1
    new_nuc = temp_row['variable']

    # Find the gene this position falls under
    nth_protein = sum(genome_pos >= start_list)-1
    # Find the position within the gene
    protein_pos = genome_pos - start_list[nth_protein]-1
    polyprotein_pos = genome_pos - start_list[1]-1
    
    if nth_protein not in [0, 12]:

        # Find the position within the codon of this position
        aa_pos = (protein_pos+1) % 3 


        # Return the amino acid start position in the consensus
        aa_num = ((protein_pos+1)//3)*3 + start_list[nth_protein]
        # Return the amino acid number
        aa_value = (protein_pos+1)//3 + 1
        # Return the aa position within the polyprotein
        pp_aa_value = (polyprotein_pos+1)//3 + 1

        # Find the old and new codon from this mutation
        old_codon = consensus[aa_num:aa_num+3]
        new_codon = old_codon[:aa_pos] + new_nuc + old_codon[aa_pos+1:]
        #print(aa_num, aa_pos, old_codon, new_codon, consensus[genome_pos], proteins[nth_protein], (protein_pos+1),aa_value,pp_aa_value,nth_protein)

        # Find the old and new amino acid
        old_aa = str(Seq.Seq(old_codon).translate())
        new_aa = str(Seq.Seq(new_codon).translate())
        #print(old_aa, new_aa, old_codon, new_codon)
        # Separate the synonymous from non-synonymous mutation
        if old_aa != new_aa:
            non_synon = 'non-synonymous'
        else:
            non_synon = 'synonymous'
        
        if genome_pos < len(consensus):
            temp_file.write(','.join([str(x) for x in list(temp_row[['passage', 'cell',
                'strain', 'temperature','replicate', 'nucleotide position', 'coverage', 
                'variable', 'value', 'percent']])]))
            temp_file.write(','.join(['',consensus[genome_pos], proteins[nth_protein], str(aa_value), 
                str(pp_aa_value), old_codon.upper(), new_codon.upper(), old_aa, new_aa, non_synon]))
            temp_file.write('\n')
    else:
        if nth_protein == 0:
            temp_protein = "3' UTR"
        elif nth_protein == 12:
            temp_protein = "5' UTR"
        #print(consensus[genome_pos], temp_protein, (protein_pos+1),nth_protein)

        if genome_pos < len(consensus):
            temp_file.write(','.join([str(x) for x in list(temp_row[['passage', 'cell',
                'strain', 'temperature','replicate', 'nucleotide position', 'coverage', 
                'variable', 'value', 'percent']])]))
            new_line = ','.join(['',consensus[genome_pos], temp_protein, 'NA', 
                'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
            temp_file.write(new_line)
            temp_file.write('\n')
temp_file.close()
print('Done')
