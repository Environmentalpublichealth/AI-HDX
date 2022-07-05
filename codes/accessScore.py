#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 15:38:40 2021

@author: jiali
Calculate surface accessibility from PDB
Usage: python accessScore.py <ID> <PDB_file.pdb>
"""
import sys
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
PDB_ID = sys.argv[1]
PDB_file = sys.argv[2]

p = PDBParser()
structure = p.get_structure(PDB_ID,PDB_file)
model = structure[0]

dssp = DSSP(model, PDB_file,dssp='/mkdssp')

# calculate the number of AA in the protein 
aa_list = []
for chain in structure.get_chains():
    AA_len = len([_ for _ in chain.get_residues() if PDB.is_aa(_)])
    
    for id in chain.get_residues():
        if PDB.is_aa(id):
            aa_list.append(id)
print(len(aa_list))

# write the dssp result into a table
with open(PDB_file+".dssp.txt","w") as out:
    for i in range(0,len(aa_list)):
        a_key = list(dssp.keys())[i]
        for item in list(dssp[a_key]):
            out.write(str(item)+"\t")
        out.write("\n")