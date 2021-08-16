from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pandas as pd # for data manipulation
import numpy as np # for data manipulation
from sklearn.linear_model import LinearRegression # for building a linear regression model
from sklearn.svm import SVR # for building support vector regression model

def vector(protein):
    aminoacids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    aminovec = []
    for i in aminoacids:
        aminovec.append(protein.count(i))
    return aminovec

def drate(uptake, maxuptake):
    return 100 * uptake/(maxuptake*0.7)

def removeOverlap(start_list, end_list):
  uniq_start = []
  uniq_end = []
  uniq_idx = []
  for idx1, val1 in enumerate(start_list[:-1]):
    overlap = False
    uniq_idx.append(idx1)
    uniq_start.append(val1)
    uniq_end.append(end_list[idx1])
    for val2 in start_list[idx1+1:]:
      if val1 >= val2 and end_list[idx1] <= end_list[idx1+1]:
        overlap = True
        break
    if overlap == True:
        uniq_start.pop()
        uniq_end.pop()
        uniq_idx.pop()
  return uniq_start, uniq_end, uniq_idx


def training(X, y):
    X, y = np.array(X), np.array(y)
    #print(); print(y.shape)

    # Creating training and test sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.30,
                                                        random_state=42)
    #print(); print(X_train.shape)
    #print(); print(X_test.shape)
    #print(); print(y_train.shape)
    #print(); print(y_test.shape)
    model2 = SVR(kernel='rbf', C=1, epsilon=10) # set kernel and hyperparameters
    svr = model2.fit(X_train, y_train)
    #print(y_test, y_train)
    #print(len(X_train))
    #print(len(X_test))
    #print(len(y_train))
    #print(len(y_test))
    #y_train[:20] = y_train[:20].reshape((20, 2))
    #X_train[0] = X_train[0].reshape((20, 1))
    #print(r2_score(X_train[0], y_train[:20]))
