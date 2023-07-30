from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import *
from Bio.PDB.Chain import Chain
from Bio.PDB.internal_coords import *
from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
from Bio.PDB.SCADIO import write_SCAD
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBIO import PDBIO
import numpy as np

ppb = PPBuilder()
parser = PDBParser()

structure_id = "thing"
filename = "ERb.preHDX.pdb"
structure = parser.get_structure(structure_id, filename)
myChain = structure[0]["A"]

resultDict = structure_rebuild_test(myChain)
assert resultDict["pass"] == True

myChain.atom_to_internal_coordinates(verbose=True)
myChain.internal_to_atom_coordinates()
io = PDBIO()
io.set_structure(myChain)
io.save("myChain.pdb")

IC_Residue.accept_atoms = IC_Residue.accept_backbone
IC_Chain.MaxPeptideBond = 4.0
myChain.internal_coord = None
write_SCAD(myChain, "myChain.scad", scale=10.0)