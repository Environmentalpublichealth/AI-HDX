"""
Created on August 18, 2021
Author: Jiali
A linear regression prediction function, used in the main processing script.
"""
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
# check if the input value is exist, if yes, then use the original xlist, if not, add the input to form a new list.

def new_list(xlist, x_pred):
  if x_pred in xlist:
    return xlist
  else:
    xlist[-1] = x_pred
    return xlist

def linreg(xlist, ylist, predict): #xlist are exposure values, ylist are uptake values
    predictedvals = []
    X= np.array(xlist)
    y = np.array(ylist)
    y = y.reshape(-1,1)
    X = X.reshape(-1,1)

    poly_reg=PolynomialFeatures(degree=3)
    X_poly=poly_reg.fit_transform(X)
    poly_reg.fit(X_poly,y)
    lin_reg2=LinearRegression()
    lin_reg2.fit(X_poly,y)  
    
    new_xlist = new_list(xlist, predict)
    for val in new_xlist:
      X_pred = poly_reg.fit_transform( np.array([[val]]) )
      y_pred = lin_reg2.predict(X_pred)
      predictedvals.append(y_pred[0][0])
    
    return new_xlist, predictedvals