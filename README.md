# AI-HDX

## Run AI-HDX on google colab [AI-HDX colab](https://colab.research.google.com/github/Environmentalpublichealth/AI-HDX/blob/main/AI_HDX.ipynb)

## Preprocessing
* [MSA encoding protein sequence](https://github.com/Environmentalpublichealth/AI-HDX/blob/main/Documentations/MSA_embedding.md)
* [Compute DSSP](https://github.com/Environmentalpublichealth/AI-HDX/blob/main/Documentations/dssp.md)
* Combine MSA encoding and surface accessibility and AA properties into a matrix
```bash
conda activate mypython3
python codes/MSA_embedding.py <hhm file> <pdb.dssp.txt> <output.vector.txt>
```
Each AA of a protein will be represented by the 36 numbers in the `output.vector.txt`.           
Use the `output.vector.txt` to encode the peptide sequences in the HDX result table.

If your protein is missing a structure PDB file, AlphaFold or ReseTTAFold are good sources for protein structure predictions!`

## Create peptide fragments
Using the Web-based tool https://web.expasy.org/peptide_mass/.      
* Input - protein sequence
* Output - save the output into a csv file with start, stop, seq, HDX rate like this:
```csv
8,21,DSASSPPYSVNQNL,8.93
8,21,DSASSPPYSVNQNL,13.84
34,45,YVDKLSSSGASW,6.55
34,45,YVDKLSSSGASW,8.23
46,60,HTEWTWSGGEGTVKS,9.09
46,68,HTEWTWSGGEGTVKSYSNSGVTF,5.16
```

# AI-HDX prediction. 
We will need two inputs, the sequence embedding vector `vector.output.txt` and the peptide fragment table `protein.csv`.  
Run AI-HDX on google colab:
[AI-HDX colab](https://colab.research.google.com/github/Environmentalpublichealth/AI-HDX/blob/main/AI_HDX.ipynb)

### Publication
Yu, J., Uzuner, U., Long, B., Wang, Z., Yuan, J. S., & Dai, S. Y. (2023). Artificial intelligence-based HDX (AI-HDX) prediction reveals fundamental characteristics to protein dynamics: Mechanisms on SARS-CoV-2 immune escape. Iscience, 26(4).
https://doi.org/10.1016/j.isci.2023.106282

