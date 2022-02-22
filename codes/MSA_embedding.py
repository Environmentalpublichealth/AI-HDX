"""
Created on July 21, 2021
Author: Jiali
Usage: Take the HHBlits output hhm format, and encode the protein sequence
python3 environment
Run:
python MSA_embedding.py <proteinID.hhm> <pdb.dssp.txt> <matrix_file.csv> <output txt file for embedding vector>
"""
import sys
import math
import numpy as np
import pandas as pd

hhm_file = sys.argv[1]
dssp_file = sys.argv[2]
output_file = sys.argv[3]

# add amino acid properties matrix to the final vector
## read the matrix table and turn it into a dictionary
### convert sequence to number vector
AA_table = pd.read_csv("dataset/HDMDvector.csv")
AA_array = AA_table.set_index('ID').T.to_dict('list')
### convert the list in dictionary to array
for key, value in AA_array.items():
    value = list(map(str,value))
    AA_array[key] = value

with open(hhm_file) as f, open(dssp_file) as df, open(output_file ,"w") as out:
    # read the dssp file line by line and add the HMM values in
    sequence_matrix = []
    for line in df:
        hhm_vector = []
        content = line.split()
        hhm_vector = hhm_vector + content[:2] + [content[3]] # extract sequence AA and accessible surface area from dssp.txt
        property_m = AA_array[content[1]] # assign each AA with its own property matrix, adopted from HDMD
        hhm_vector = hhm_vector + property_m # AA property features
        sequence_matrix.append(hhm_vector)
    
    for i in f:
        if i.startswith("#"):
            break
    # start from the lines with HMM values
    for i in range(4):
        f.readline()
    # skip the NULL and HMM headers
    lines = f.read().split("\n")
    # print(len(lines)) ## The file consist of three lines for each AA, first line is the HMM number against each AA,
    ## second line is the 10 conversion values, and the last line is empty. Group the three lines into one AA representative.

    for idx in range(0,int((len(lines)-2)/3)):
        # print(lines[idx*3])
        first_line = lines[idx*3].replace("*","99999") # The * symbol is like NA, so here we assigned a big number to it
        next_line = lines[idx*3+1].replace("*","99999")
        content1 = first_line.strip().split()
        content2 = next_line.strip().split()
        for val1 in content1[2:-1]:
            hhm_val1 = 10/(1 + math.exp(-1 * int(val1)/2000))
            # hhm_vector.append(str(hhm_val1))
            sequence_matrix[idx].append(str(hhm_val1))
        for val2 in content2:
            hhm_val2 = 10/(1 + math.exp(-1 * int(val2)/2000))
            # hhm_vector.append(str(hhm_val2))
            sequence_matrix[idx].append(str(hhm_val2))
        # output the vector for each AA in the protein
    for item in sequence_matrix:
        out.write("\t".join(item)+"\n")