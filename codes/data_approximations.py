import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression


def linreg(xlist, ylist, predict): #xlist are exposure values, ylist are uptake values
    predictedvals = []
    for i in range(len(xlist)-1):
        if xlist[i] <= predict and xlist[i+1] > predict:
            if xlist[i] == predict:
                predictedvals.append(ylist[i])            
            else:
                try:
                    predictedvals.append( ( (ylist[i+1] - ylist[i])/(xlist[i+1] - xlist[i]) ) * (predict - xlist[i]) + ylist[i])
                except ZeroDivisionError:
                    predictedvals.append("ZeroDivisionError")
    return predictedvals

    
