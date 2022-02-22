"""
Created on Feb 1, 2022
Author: Jiali
"""
from Bio import PDB
from Bio.PDB.ResidueDepth import ResidueDepth
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import residue_depth
from Bio.PDB.ResidueDepth import get_surface
import sys

PDBfile = sys.argv[1]

parser = PDBParser()
structure = parser.get_structure('xyn1',PDBfile)
model = structure[0]
surface = get_surface(model)
chain = model["A"]
# calculate the number of AA in the protein
for chain in structure.get_chains():
    AA_len = len([_ for _ in chain.get_residues() if PDB.is_aa(_)])

for i in range(1,AA_len+1):
    rd = residue_depth(chain[i], surface)
    print("residue "+ str(i)+ ": "+ str(rd))