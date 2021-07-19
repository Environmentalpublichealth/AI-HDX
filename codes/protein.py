import csv
import pandas as pd
import numpy as np

from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression



manager = open("manager.txt", "r")

for row in manager:
    files = row.split(": ")
    data = open(str(files[0]) + ".csv", "r")
    res = open(str(files[0]) + "res.csv", "w")

    state, sequence, maxuptake, uptake, exposure = -1, -1, -1, -1, -1
    for rows in data: 
        rows = rows.split(",")
        for i in range(len(rows)):
            if rows[i] == "State":
                state = i
            if rows[i] == "Sequence":
                sequence = i
            if rows[i] == "MaxUptake":
                maxuptake = i
            if rows[i] == "Uptake":
                uptake = i
            if rows[i] == "Exposure":
                exposure = i
        break


    if maxuptake == -1: 
        x = []
        for rows in data:
            rows = rows.split(",")
            x.append( (len(rows[sequence])) - 1 - rows[sequence].count("P")  )
        writer = csv.writer(res)
        writer.writerow(["Sequence", "State", "Exposure", "Uptake", "MaxUptake"])

    data.close()
    data = open(str(files[0]) + ".csv", "r")
    i = 0
    for rows in data:
        i+=1
        rows = rows.split(",")
        writer = csv.writer(res)
        if (files[1].strip() == rows[state]) and (int(float((rows[exposure]))) == 60 or int(float((rows[exposure]))) == 1 or int(float((rows[exposure]))) == 3600): #skipping the ones that are supposed to be skipped
            writer.writerow([rows[sequence], rows[state], rows[exposure], rows[uptake], x[i]])

    data.close()
    data = open(str(files[0]) + ".csv", "r")
    
    
    proteinNames = [] #put the states into list called proteinNames, then changed into set (which deleted all duplicates), and then changed it back to a list
    for rows in data:
        rows = rows.split(",")
        proteinNames.append(rows[state])
    proteinNames = list(set(proteinNames))
    
    datadict = {} #dictionary
    # {a:b, c:d, }
    data.close()
    data = open(str(files[0]) + ".csv", "r")
    i = 0
    for rows in data:
        i+=1
        rows = rows.split(",")
        for j in range(len(proteinNames)):
            if rows[state] == proteinNames[j]:
                datadict[str(j) + " " + str(i) + " " + str(rows[exposure])] = rows[uptake]

def polyReg(X, y):
    for i in range(len(X)):
        X[i] = float(X[i])
    for i in range(len(y)):
        y[i] = float(y[i])
    y = np.array(y)
    y = y.reshape(-1,1)
    X = np.array(X)
    X = X.reshape(-1,1)

    #lin_reg=LinearRegression()
    #lin_reg.fit(X,y)

    poly_reg=PolynomialFeatures()
    X_poly=poly_reg.fit_transform(X)
    poly_reg.fit(X_poly,y)
    lin_reg2=LinearRegression()
    lin_reg2.fit(X_poly,y)

    j = lin_reg2.predict( poly_reg.fit_transform( [X[1]] ))
    return j

X, y = [], []
for i in datadict:
    a = i.split(" ")
    
    if a[2].isalpha() == True:
        continue
    X.append( [a[2], a[0]] )
    y.append(datadict[i])

h, g = [], []
for i in range(len(X)-1):
    if x[i] == x[i+1]:
        h.append(x[0])
        g.append(y[i])
    else:
        print(polyReg(h, g))
        h, g = [], []

manager.close()
data.close()
res.close()

