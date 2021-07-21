"""
Created on July 21, 2021
Author: Jiali
Usage: Take the HHBlits output hhm format, and encode the protein sequence
"""
import sys
import math
import numpy as np

hhm_file = sys.argv[1]

with open(hhm_file) as f:
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
    seq_array = np.empty((0,30))
    for idx in range(0,int((len(lines)-2)/3)):
        # print(lines[idx*3])
        first_line = lines[idx*3].replace("*","99999") # The * symbol is like NA, so here we assigned a big number to it
        next_line = lines[idx*3+1].replace("*","99999")
        content1 = first_line.strip().split()
        content2 = next_line.strip().split()
        hhm_vector = []
        # print(content1)
        # print(content2)
        for val1 in content1[2:-1]:
            hhm_val1 = 10/(1 + math.exp(-1 * int(val1)/2000))
            hhm_vector.append(hhm_val1)
        for val2 in content2:
            hhm_val2 = 10/(1 + math.exp(-1 * int(val2)/2000))
            hhm_vector.append(hhm_val2)
        ## combine the numbers from first and second lines, to get a vector with 30 values, we will use this matrix to represent each AA for the protein.
        seq_array = np.vstack((seq_array, np.array(hhm_vector)))
    print(seq_array.shape)        