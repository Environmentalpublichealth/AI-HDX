"""
Created on July 27, 2023
Author: Jiali
Function to run AI-HDX prediction for shiny app
"""
# import all the necessary libraries
import tensorflow as tf
import pandas as pd
import numpy as np
from embedding import seq_embedding

# load the pre-trained models from
model1 = tf.keras.models.load_model("../models/NNmodel_m1")
model2 = tf.keras.models.load_model("../models/NNmodel_m2")
model3 = tf.keras.models.load_model("../models/NNmodel_m3")
model4 = tf.keras.models.load_model("../models/NNmodel_m4")
model5 = tf.keras.models.load_model("../models/NNmodel_m5")
model6 = tf.keras.models.load_model("../models/NNmodel_m6")
model7 = tf.keras.models.load_model("../models/NNmodel_m7")
model8 = tf.keras.models.load_model("../models/NNmodel_m8")
model9 = tf.keras.models.load_model("../models/NNmodel_m9")
model10 = tf.keras.models.load_model("../models/NNmodel_m10")

## Load the pre-trained confidence index score based on the error rate between predictions and ground true in validation sets
confidx = [0.435003467810392,0.584915875217703,0.691762575882027,0.707836210244849,0.626804716837474,0.591011203879881,0.532705643129479,0.426545116775042,0.251348832172361,0.36508649008649]

def pred_range(a):
  for i in range(0,10):
    if a >i*0.1 and a <= (i+1)*0.1:
      CI = confidx[i]
      return CI

def prediction(prot1, df1):
    x_test = prot1.reshape(prot1.shape[0], 30, 36)
    out_df = df1
    models = [model1, model2, model3, model4, model5, model6, model7, model8, model9, model10]
    colnames = ["model1","model2","model3","model4", "model5","model6","model7","model8","model9", "model10"]
    for idx, model in enumerate(models):
        y_pred = model.predict(x_test)
        out_df[colnames[idx]] = y_pred

# calculate average prediction values(mean) and standard deviation(SD)
    out_df["average"] = out_df[colnames].mean(axis=1)
    out_df["SD"] = out_df[colnames].std(axis=1)

# add CI score to each prediction
    scores=[]
    for idx, pred in enumerate(out_df["average"]):
        conf_score = pred_range(pred)
        scores.append(conf_score)
  
    out_df["CI"]=np.array(scores) 
    return out_df

## run embedding and prediction
# prot, df = seq_embedding("../example/ESR2.csv","../example/Q92731.vector.txt")
# output_df = prediction(prot, df) 
# print(output_df)