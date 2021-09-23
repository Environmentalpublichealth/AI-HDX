"""
Created on Sep 22, 2021
Author: Jiali
Guassian Network Model calculation from protein structure. 
python ANM.py <input.pdb> <output.matrix>
"""


from prody import *
from pylab import *
ion()
import sys
import numpy as np
from sklearn.decomposition import PCA
pca = PCA(n_components = 20)

# file_in = sys.argv[1]
# file_out = sys.argv[2]

# parse PDB structure
def GNM_fluc(filename):
    protein = parsePDB(filename)    

    calphas = protein.select("protein and name CA")

    # initial GNM analysis
    gnm = ANM('GNM analysis')
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(40)
    flucs = calcSqFlucts(gnm)
    eigvec = gnm.getEigvecs().round(3)
    eigval = gnm.getEigvals().round(3)
    return eigval

def ANM_fluc(filename):
    protein = parsePDB(filename)    

    calphas = protein.select("protein and name CA")

    # initial GNM analysis
    anm = ANM('GNM analysis')
    anm.buildHessian(calphas)
    anm.calcModes(40)
    flucs = calcSqFlucts(anm)
    eigvec = anm.getEigvecs().round(3)
    eigval = anm.getEigvals().round(3)
    return eigval

proteins = ["1xyn.pdb","4xi2.pdb", "6sbw.pdb", "5iv9", "1ks4.pdb","O66490.e2e.pdb","AF-P69441-F1-model_v1.pdb"]
data = np.empty([0,40])

for filename in proteins:
    print(filename)
    # residue_fluc = GNM_fluc(filename)
    # pca.fit(residue_fluc)
    # pca_fluc = pca.singular_values_
    # print(pca_fluc.shape)
    eigenvalue = ANM_fluc(filename)
    # data = np.vstack((data, res_fluc))
    data = np.vstack((data, eigenvalue))

from sklearn.cluster import KMeans
kmeans = KMeans(init = "random", n_clusters=3, n_init = 10, max_iter=10, random_state= 0)
kmeans.fit(data)
print(kmeans.labels_)

