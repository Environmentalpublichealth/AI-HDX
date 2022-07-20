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
No error message print on the screen!        
**Note!** But there are errors in the log folder, need to adjust the MSA preposcessing bash script and `run_e2e_ver.sh` to correct the path of database and conda environment source.

1. edit `run_e2e_ver.sh`:
```txt
line 16: export PIPEDIR="/scratch/user/jialiyu/RoseTTAFold" # change it for the pathway that you git clone RoseTTAFold
line 30: source activate /sw/hprc/sw/Anaconda3/2021.05/envs/RoseTTAFold # use the conda env installed by the HPRC staff
line 54: DB="/scratch/data/bio/rosetta/pdb100_2021Mar03/pdb100_2021Mar03" # put the correct path for pdb database (provided by HPRC staff)
```
2. edit `input_prep/make_msa.sh`
```txt
line 12: DB="/scratch/data/bio/rosetta/UniRef30_2020_06/UniRef30_2020_06"
line 13: MYDB="/scratch/data/bio/rosetta/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
```
3. edit `input_prep/make_ss.sh`
```txt
line 11: /sw/eb/sw/csblast/2.2.3-Linux64/bin/csbuild # the correct path for csblast program
line 11: /sw/eb/sw/csblast/2.2.3-Linux64/data/K4000.crf
```
