## Extract amino acid surface area from PDB file
Search for the right protein and get the correct PDB files from either PDB or AlphaFold database.
### Compute surface accessible area
```bash
# in the conda environment, get biopython and dssp installed
conda activate mypython3
conda install -c salilab dssp # only need to run it once
conda install -c conda-forge biopython # only need to run it once
# Calculate accessible area using the pdb file and output into a new txt file
python codes/accessScore.py PDB_ID PDB/PDB_ID_file.pdb
```
A `PDB_ID_file.pdb.dssp.txt` will be created in the same directory that has the pdb file.
