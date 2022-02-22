"""
Created on Sep 22, 2021
Author: Jiali
Guassian Network Model calculation from protein structure. 
python ANM.py <input.pdb> <HDX result table.csv>
"""


from prody import *
from pylab import *
ion()
import sys
import math
import numpy as np
from scipy import stats
import pandas as pd
from sklearn.decomposition import PCA
pca = PCA(n_components = 20)

file_in = sys.argv[1]
HDX_results = sys.argv[2]

# parse PDB structure
def GNM_fluc(filename):
    protein = parsePDB(filename)    

    calphas = protein.select("protein and name CA")
    N = len(calphas)
    # initial GNM analysis
    gnm = GNM('GNM analysis')
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(3*N-6)
    flucs = calcSqFlucts(gnm)
    eigvec = gnm.getEigvecs().round(3)
    eigval = gnm.getEigvals().round(3)
    return eigval, N

def ANM_fluc(filename):
    protein = parsePDB(filename)    

    calphas = protein.select("protein and name CA")
    N = len(calphas)
    # initial GNM analysis
    anm = ANM('ANM analysis')
    anm.buildHessian(calphas)
    anm.calcModes((100))
    # anm.calcModes((3*N-6))
    flucs = calcSqFlucts(anm)
    eigvec = anm.getEigvecs().round(3)
    eigval = anm.getEigvals().round(3)
    return eigval, N

def HDX_mean(HDX_file):
    # read HDX data file
    datafile = pd.read_csv(HDX_file, header=None)
    rate_array = datafile.loc[:,3]
    meanHDX = np.mean(rate_array)
    return meanHDX

proteins = []
with open(file_in) as f:
    for line in f:
        PID = line.strip()
        proteins.append(PID)

#proteins = ["1xyn.pdb","1ks4.pdb","filter/P36218.e2e.pdb","filter/O74705.e2e.pdb","1anf.pdb","filter/AF-P0AEX9-F1-model_v1.pdb","1gjh.pdb","filter/AF-P10415-F1-model_v1.pdb"]#,"C4T9I0.e2e.pdb", "AF-Q04207-F1-model_v1.pdb", "AF-P0A9D4-F1-model_v1.pdb", "O66490.e2e.pdb","AF-P69441-F1-model_v1.pdb"]
data = np.empty([0,100])
## Calculate flexibility from ANM for all the proteins in the list.
flex_score = []
for filename in proteins:
    print(filename)
    # residue_fluc = GNM_fluc(filename)
    # pca.fit(residue_fluc)
    # pca_fluc = pca.singular_values_
    # print(pca_fluc.shape)
    eigenvalue, N = ANM_fluc(filename)
    # print(eigenvalue)
    #scores = np.sum((eigenvalue)**2) 
    #flex_score.append(math.sqrt(scores/N))
    scores = np.mean(1/eigenvalue**2)
    flex_score.append(math.log(math.sqrt(scores)))
# print(flex_score)

    # data = np.vstack((data, res_fluc))
    data = np.vstack((data, eigenvalue))
## Calculate geometric mean for all HDX datasets in the list.
# HDXvalue = []
# with open(HDX_results) as dataset:
#     for line in dataset:
#         fileID = line.strip()
#         result = HDX_mean("/Users/jiali/Desktop/Jiali/TAMU/Dynamics/Seq2HDX/filtered/Train_August/"+fileID+".csv")
#         HDXvalue.append(result)

# ## generate scatter plot of HDX and flex score
# import matplotlib.pyplot as plt

# plt.figure(figsize=(4,4))
# for i, ID in enumerate(proteins):
#     plt.plot(flex_score[i],HDXvalue[i], "o", color="blue")
#     plt.text(flex_score[i]+0.03, HDXvalue[i], ID, fontsize=8)

# plt.xlabel("NMA flexibility")
# plt.ylabel("mean HDX rate")
# plt.savefig("cor_plot.pdf", dpi=300)

# from sklearn.cluster import KMeans
# kmeans = KMeans(init = "random", n_clusters=3, n_init = 10, max_iter=10, random_state= 0)
# kmeans.fit(data)
# print(kmeans.labels_)

# run t-SNE
from sklearn.manifold import TSNE
# perplexity parameter can be changed based on the input datatset
# dataset with larger number of variables requires larger perplexity
# set this value between 5 and 50 (sklearn documentation)
# verbose=1 displays run time messages
# set n_iter sufficiently high to resolve the well stabilized cluster
# get embeddings
tsne_em = TSNE(n_components=2, perplexity=30.0, n_iter=1000, verbose=1).fit_transform(data)
# plot t-SNE clusters
from bioinfokit.visuz import cluster # need to get conda install bioinfokit
cluster.tsneplot(score=tsne_em)
# plot will be saved in same directory (tsne_2d.png) 

# label the samples 
# color_class = df['class'].to_numpy()
# cluster.tsneplot(score=tsne_em, colorlist=color_class, legendpos='upper right', legendanchor=(1.15, 1) )