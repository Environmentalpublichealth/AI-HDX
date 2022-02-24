"""
Created on Feb 24, 2022
Author: Jiali
Usage, replace the beta value in pdb file to HDX rates
python recode_beta. py <pdf file> <HDX.csv> <output.pdb>
"""
import sys
from Bio import PDB
from Bio.PDB import PDBParser
structurePDB = sys.argv[1]
HDX = sys.argv[2]
HDXPDB = sys.argv[3]

# calculate the AA number in protein
p = PDBParser()
PDB_ID = "1xyn"
structure = p.get_structure(PDB_ID,structurePDB)
for chain in structure.get_chains():
    AA_len = len([_ for _ in chain.get_residues() if PDB.is_aa(_)])
# AA_len is the number
# generate a list of AA HDX rates
HDXlist = ["0"]*AA_len
with open(HDX) as data:
    for HDXline in data:
        HDXcontent = HDXline.strip("\ufeff").rstrip().split(",")
        peptide = range(int(HDXcontent[0])-1,int(HDXcontent[1]))
        for i in peptide:
            HDXlist[i] = HDXcontent[3]
# replace beta value with HDX data

with open(structurePDB) as f, open(HDXPDB,"w") as out:
    for PDBline in f:
        PDBcontents = PDBline.strip("\n").split()
        if PDBcontents[0] == "ATOM":
            resID = int(PDBcontents[5])-1
            PDBcontents[10] = HDXlist[resID]
            # print(PDBcontents)
            # format the table into pdb format: https://cupnet.net/pdb-format/
            PDBtext = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(PDBcontents[0],int(PDBcontents[1]),PDBcontents[2],"",PDBcontents[3],PDBcontents[4],int(PDBcontents[5]),"",float(PDBcontents[6]),float(PDBcontents[7]),float(PDBcontents[8]),float(PDBcontents[9]),float(PDBcontents[10]),PDBcontents[11],"")
            out.write(PDBtext+"\n")
        else:
            out.write(PDBline)

# for chain in structure:
#     for residue in chain:
#         for atom in residue