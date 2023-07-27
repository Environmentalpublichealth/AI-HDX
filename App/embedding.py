"""
created on July 26, 2023
Author: Jiali
The embedding function for AI-HDX shiny app
"""
# embedding protein
import pandas as pd
import numpy as np
from tensorflow.keras.preprocessing.sequence import pad_sequences

def seq_embedding(HDX_file, vector_file):
  # read HDX data file
  datafile = pd.read_csv(HDX_file, header=None)
  # filter peptides > 30 AA
  max_len = 30
  df = datafile[datafile.loc[:,2].str.len() < max_len]
  start_pos = df.loc[:,0].tolist()
  end_pos = df.loc[:,1].tolist()
  # read embedding file
  embedding = pd.read_table(vector_file, header=None)
  embed_array = embedding.loc[:,2:].to_numpy()

  row_size = len(start_pos)
  nfactor = 36

  input_array = np.empty((row_size,nfactor,max_len)) # create an empty array
  x=0
  for i, numb in enumerate(start_pos):
    seq_array = embed_array[numb-1:end_pos[i]]
    seq_arrayT = np.transpose(seq_array)
    padded_seq = pad_sequences(seq_arrayT, maxlen=max_len, padding="post",dtype="float64")
    input_array[x,:,:] = padded_seq
    x += 1
  

  return input_array, df

# import pre-processed sequence embedding vector and the peptide fragment tables. Save them into two variables, 'prot' is the peptide embedding array, 'df' is the table to store outputs.
# prot, df = seq_embedding("../example/ESR2.csv","../example/Q92731.vector.txt")
# print(df)