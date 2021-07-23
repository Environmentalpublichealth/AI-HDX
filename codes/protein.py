import csv
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression



manager = open("manager.txt", "r") #opening the file that has the names of all the files to run

def important_variables(data):
    '''
    This iterates through the data and extracts the state, sequence, maxuptake, uptake, and exposure by its index.
    Since only the first row of values is needed. The iteration is halted after the first row is read.
    It then returns those values as a list.
    '''
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
    return [state, sequence, maxuptake, uptake, exposure]


def no_maxuptake(data, sequence):
    '''
    If no maxuptake exists, this function will be ran and a maxuptake will be generated.
    '''
    maxuptake = []
    for rows in data:
        rows = rows.split(",")
        maxuptake.append( (len(rows[sequence])) - 1 - rows[sequence].count("P")  )
    #writer = csv.writer(res)
    #writer.writerow(["Sequence", "State", "Exposure", "Uptake", "MaxUptake"])
    return maxuptake

def bad_exposures(data, state, sequence, maxuptake, uptake, exposure):
    '''
    This function will filter out the unnecessary exposure times.
    '''
    overall_data = []
    data.close()
    data = open(str(files[0]) + ".csv", "r")
    i = 0
    for rows in data:
        i+=1
        rows = rows.split(",")
        writer = csv.writer(res)
        if (files[1].strip() == rows[state]) and (int(float((rows[exposure]))) == 60 or int(float((rows[exposure]))) == 1 or int(float((rows[exposure]))) == 3600): #skipping the ones that are supposed to be skipped
            overall_data.append([rows[sequence], rows[state], rows[exposure], rows[uptake], maxuptake[i]])
    return overall_data
    
def protein_dict(data, state, exposure, uptake):
    '''
    Generates a dictionary with all the necessary information for polynomial and linear regression.
    '''
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
                datadict[str(j) + " " + str(i) + " " + str(rows[exposure])] = rows[uptake] #puts state, exposure, and uptake data into dictionary
    return datadict

def dictionary_transfer(datadict, maxuptake):
    '''
    Transfers the information from the dictionary into a type of data form that can be used for polynomial/linear regression.
    '''
    x, y = [], [] #x is exposure; y is uptake
    for i in datadict:
        a = i.split(" ")
    
        if a[2].isalpha(): #we want to skip the first row because it's title names. Sometimes the first row isn't a title name, though, so we only skip it if it contains words
            continue
        x.append( [a[2], a[0]] ) #exposure, protein number
        y.append(datadict[i]) #uptake

    h, g = [], []
    all_x_and_y_data = []

    j = 0 #***
    for i in range(len(x)-1): #delete all the lines marked with *** if you want to run the [0,0] coordinate
        j+=1 #***
        if j-1 == 0: #***
            continue #***
        if x[i][1] == x[i+1][1]: #if it is the same protein, then keep adding it to the list
            h.append(x[i][0]) #exposure
            g.append(y[i])
        else: #otherwise, plug the data into the poly_reg function
            h.append(x[i][0])
            g.append(y[i])
            all_x_and_y_data.append([h,g])
            #print(poly_reg(h, g))
            h, g = [], [] #then clear the two lists and begin the process over again
            j = 0 #***
    return all_x_and_y_data


def poly_reg(X, y):
    '''
    Polynomial regression occurs here.
    I have no idea what's going on here so no detailed docstring for this one.
    '''
    for i in range(len(X)):
        X[i] = float(X[i])
    for i in range(len(y)):
        y[i] = float(y[i])
    y = np.array(y)
    y = y.reshape(-1,1)
    X = np.array(X)
    X = X.reshape(-1,1)

    lin_reg=LinearRegression()
    lin_reg.fit(X,y)

    poly_reg=PolynomialFeatures()
    X_poly=poly_reg.fit_transform(X)
    poly_reg.fit(X_poly,y)
    lin_reg2=LinearRegression()
    lin_reg2.fit(X_poly,y)

    j = lin_reg2.predict( poly_reg.fit_transform( [X[1]] ))
    
    plt.scatter(X,y,color='red')
    plt.plot(X,lin_reg.predict(X),color='blue')
    plt.title('Truth or bluff(Linear Regression)')
    plt.xlabel('Exposure')
    plt.ylabel('Uptake')
    plt.show()

    return j

for row in manager: 
    files = row.split(": ")
    data = open(str(files[0]) + ".csv", "r")
    res = open(str(files[0]) + "res.csv", "w")
    
    a = important_variables(data) #state, sequence, maxuptake, uptake, exposure
    #if a[2] == -1:
    #    print(no_maxuptake(data, a[1]))

    #b = bad_exposures(data, a[0], a[1], a[2], a[3], a[4]) #overall_data
    
    c = protein_dict(data, a[0], a[4], a[3]) #datadict
    
    d = dictionary_transfer(c, a[2]) #all_x_and_y_data, in format of [ [x1,y1], [x2,y2] ]
    #print(d)
    for i in d:
        print(poly_reg(i[0], i[1]))

manager.close()
data.close()
res.close()
