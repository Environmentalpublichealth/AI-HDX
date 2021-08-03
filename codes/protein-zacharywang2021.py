import csv
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression

class variable_capture(object):
    def __init__(self, data):
        self.data = data
    def important_variables(self):
        '''
        This iterates through the data and extracts the state, sequence, maxuptake, uptake, and exposure by its index.
        Since only the first row of values is needed. The iteration is halted after the first row is read.
        It then returns those values as a list.
        '''
        state, sequence, maxuptake, uptake, exposure = -1, -1, -1, -1, -1
        for rows in self.data: 
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

class data_manipulation(variable_capture): # [state, sequence, maxuptake, uptake, exposure]
    def __init__(self, data, labels):
        variable_capture.__init__(self, data)

        self.state = labels[0]
        self.sequence = labels[1]
        self.maxuptake = labels[2]
        self.uptake = labels[3]
        self.exposure = labels[4]

    def no_maxuptake(self):
        '''
        If no maxuptake exists, this function will be ran and a maxuptake will be generated.
        '''
        maxuptake = []
        for rows in self.data:
            rows = rows.split(",")
            maxuptake.append( (len(rows[self.sequence])) - 1 - rows[self.sequence].count("P")  )
        return maxuptake
    def bad_exposures(self): #data, state, sequence, maxuptake, uptake, exposure
        '''
        This function will filter out the unnecessary exposure times.
        '''
        overall_data = []
        i = 0
        for rows in self.data:
            i+=1
            if (files[1].strip() == rows[self.state]) and (int(float((rows[self.exposure]))) == 60 or int(float((rows[self.exposure]))) == 1 or int(float((rows[self.exposure]))) == 3600): #skipping the ones that are supposed to be skipped
                overall_data.append([rows[self.sequence], rows[self.state], rows[self.exposure], rows[self.uptake], self.maxuptake[i]])
        return overall_data
    def protein_dict(self, data):
        '''
        Generates a dictionary with all the necessary information for polynomial and linear regression.
        '''
        
        proteinNames = [] #put the states into list called proteinNames, then changed into set (which deleted all duplicates), and then changed it back to a list
        for rows in data:
            rows = rows.split(",")
            proteinNames.append(rows[self.state])

        proteinNames = list(set(proteinNames))
        
        data.close()
        data = open(str(files[0]) + ".csv", "r") 
        datadict = {} #dictionary
        # {a:b, c:d, }
        i = 0
        for rows in data:
            i+=1
            rows = rows.split(",")
            for j in range(len(proteinNames)):
                if rows[self.state] == proteinNames[j]:
                    datadict[str(j) + " " + str(i) + " " + str(rows[self.exposure])] = rows[self.uptake]
        return datadict
    
    def dictionary_transfer(self, datadict):
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

class data_results(variable_capture):
    def __init__(self, data):
        variable_capture.__init__(self, data)

    def poly_reg(self, X, y):
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
        plt.xlabel('Position Level')
        plt.ylabel('Salary')
        plt.show()

        return j

    def protein_merge(self, start, end, sequence):
        '''
        assembles the full length protein sequence 
        '''
        
        def protein_finder(merged_protein, start, end, sequence, position): #recursive function
            #if end[position] > start[len(start)-1]: #this deals with the end of the sequence
            #    return merged_protein

            for i in range(position, len(start)):
                if start[i] - end[position] > 0 and start[i] - end[position] <= 5:
                    merged_protein += sequence[i]
                    return protein_finder(merged_protein, start, end, sequence, i)
                elif start[i] - end[position] > 5:
                    break
            return merged_protein
        
        position = 0 #position is the current position of the list
        for i in range(len(start)): #it iterates down the start list, and breaks the for loop once a positive start value is reached
            if start[i] > 0:
                position = i
                break
            
        merged_protein = sequence[position]#then it takes the first protein at the position
        return protein_finder(merged_protein, start, end, sequence, position)
manager = open("manager.txt", "r")

for row in manager: 
    files = row.split(": ")
    data = open(str(files[0]) + ".csv", "r")   
    
    a = variable_capture(data)
    
    b = a.important_variables() # [state, sequence, maxuptake, uptake, exposure]
    
    c = data_manipulation(data, b)

    #d = c.no_maxuptake()
    
    #e = c.bad_exposures()
    
    f = c.protein_dict(data)
    
    g = c.dictionary_transfer(f)
    
    h = data_results(data)

    #for i in g: #MATLAB, LINEAR REGRESSION VALUES
    #    print(h.poly_reg(i[0], i[1]))
    df = pd.read_csv(str(files[0]) + ".csv")
    index_list = df.index.tolist()
    start_name_list = df["Start"].tolist()
    end_name_list = df["End"].tolist()
    sequence_list = df["Sequence"].tolist()

    i = h.protein_merge(start_name_list, end_name_list, sequence_list)
    print(i)

    
  
