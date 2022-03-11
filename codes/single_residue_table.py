"""
Created on March 09, 2022
Author: Jiali
This program convert a protein sequence into single amino acid table format matching the HDX result table.
Usage: single_residue_table.py <a protein fasta file> > output_file_name.csv
"""
import sys
from Bio import SeqIO
fasta_file = sys.argv[1]

# read input sequence
fasta_seq = SeqIO.parse(open(fasta_file),"fasta")
for fasta in fasta_seq:
    seq = fasta.seq

# print(seq)


# cut sequence into a single amino acid and write into a table
for i in range(len(seq)):
    start = i+1
    end = i+2
    AA = seq[i]
    print(start, end, AA, 0, sep = ",")