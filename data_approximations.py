import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression


def linreg(xlist, ylist, predict): #xlist are exposure values, ylist are uptake values
    predictedvals = []
    for i in range(len(xlist)):
        if xlist[i] == predict:
            predictedvals.append(xlist[i])
        else:
            predictedvals.append( ( (ylist[i] - ylist[i-1])/(xlist[i] - xlist[i-1]) ) * (predict - xlist[i-1]) + ylist[i-1])
    return predictedvals
    
