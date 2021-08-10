# Seq2HDX

## MSA encoding protein sequence
### Install HHBlits via conda
```bash
module load Anaconda3/2020.07
conda create -n hhblits
source activate hhblits
conda install -c conda-forge -c bioconda hhsuite
```
Download database
```bash
mkdir databases
cd databases
nohup wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz &
tar -xvfz UniRef30_2020_06_hhsuite.tar.gz
```
The download takes ~6 hours, because it is connected to a German server.

### Run HHBlits
```bash
module load Anaconda3/2020.07

source activate hhblits

cd /scratch/user/jialiyu/BlitSearch
for file in ./*fasta
do
BASE=$(basename $file | sed 's/.fasta//g')
hhblits -i $file -ohhm $BASE.hhm -d databases/UniRef30_2020_06/UniRef30_2020_06
mv $BASE.hhm hhm_data
done
```
Each protein is stored in a single fasta file, so I used a for loop to run all files.

Transfer `.hhm` files to local.
```bash
rsync -a grace.tamu:/scratch/user/jialiyu/BlitSearch/hhm_data .
```

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

## Combine MSA encoding and surface accessibility and AA properties into a matrix
```bash
conda activate mypython3
python codes/MSA_embedding.py <hhm file> <pdb.dssp.txt> <property_matrix.csv> <output.vector.txt>
```
Each AA of a protein will be represented by the 36 numbers in the `output.vector.txt`.           
Use the `output.vector.txt` to encode the peptide sequences in the HDX result table.

## Protein structure prediction by RoseTTAFold
Proteins with mutations and from non-model organisms rarely have protein crystal structures, not on the AlphaFold database either. We use RosseTTAFold to predict it instead. The program has been installed by the HPRC staff as a conda environment on Grace. 
```bash
git clone https://github.com/RosettaCommons/RoseTTAFold.git
# login to Grace cluster
module load Anaconda3/2021.05
source activate /sw/hprc/sw/Anaconda3/2021.05/envs/RoseTTAFold
# test if it worked
cd RoseTTAFold/example
../run_e2e_ver.sh input.fa .
```
No error message!        
**Note!** When running our own input sequence, need to adjust the MSA preposcessing bash script and `run_e2e_ver.sh` to correct the path.
