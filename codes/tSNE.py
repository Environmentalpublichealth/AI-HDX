"""
Created on Nov 11, 2021
Author: Jiali
Clustering on peptides embeddings 
python tSNE.py
"""

# First need to connect this to google drive
# embedding protein 1
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
    start = int(numb)
    end = int(end_pos[i])
    seq_array = embed_array[start-1:end]
    seq_arrayT = np.transpose(seq_array)
    padded_seq = pad_sequences(seq_arrayT, maxlen=max_len, padding="post",dtype="float64")
    input_array[x,:,:] = padded_seq
    x += 1
  

  label_array = df.loc[:,3]/100

  return input_array, label_array

# Import and combine training datasets
X = np.empty([0,36,30])
y = np.empty([0,1])
prot_ID = []
file_drive = "/Users/jiali/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/"

with open("/Users/jiali/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/all_data.txt") as f:
  for line in f:
    content = line.strip().split()
    seq_file = file_drive+"embedding/"+content[0]+".vector.txt"
    data_file = file_drive+"Train_August/"+content[1]+".csv"
    # print(seq_file, data_file, sep=",")
    prot, lab = seq_embedding(data_file, seq_file)
    X = np.vstack((X, prot))
    y = np.append(y, lab)
    prot_ID = prot_ID + [content[2]]*lab.shape[0]   

X_data = X.reshape(X.shape[0],36*30)
y_label = np.round(y, 1)

# from sklearn.decomposition import PCA
# pca = PCA(n_components = 50)
# pca.fit(X_data)
# pca_data = pca.components_

# run t-SNE
from sklearn.manifold import TSNE
# perplexity parameter can be changed based on the input datatset
# dataset with larger number of variables requires larger perplexity
# set this value between 5 and 50 (sklearn documentation)
# verbose=1 displays run time messages
# set n_iter sufficiently high to resolve the well stabilized cluster
# get embeddings
tsne_em = TSNE(n_components=3, perplexity=30.0, n_iter=1000, verbose=1).fit_transform(X_data)
# plot t-SNE clusters
from bioinfokit.visuz import cluster # need to get conda install bioinfokit
# cluster.tsneplot(score=tsne_em)
# plot will be saved in same directory (tsne_2d.png) 

# label the samples 
color_class = y_label
cluster.tsneplot(score=tsne_em, colorlist=prot_ID, legendpos='upper right', 
colordot=('#713e5a', '#63a375', '#edc79b', '#d57a66', '#ca6680', '#395B50', '#92AFD7', '#b0413e', '#4381c1', '#736ced', '#631a86'),
legendanchor=(1.35, 1) )

print(color_class)
# additional color options: '#de541e', '#022b3a', '#000000'