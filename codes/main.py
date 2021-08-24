from os import remove
import pandas as pd
from data_approximations import linreg
from protein_merge import protein_merge
from simplemodel import vector, drate, removeOverlap, training

manager = open("manager.txt", "r")

for row in manager: 
    
    files = row.split(": ")
    try:
        uptakepredict = input(files[0] +"Exposure prediction: ")
    except ValueError:
        print("Make sure your cursor is on the input line")
    df = pd.read_csv(files[0] +".csv")

    #maxuptake
    mut = []
    if "MaxUptake" not in df[0:0]:
        for index, row in df.iterrows():
            mut.append(len(str((row["Sequence"]))) - 1 - str(row["Sequence"]).count("P"))
    else:
        for index, row in df.iterrows():
            mut.append(row["MaxUptake"])

    #uptake prediction
    xlist, ylist, predict, predictedvals = [], [], float(uptakepredict), []
    for index, row in df.iterrows():
        if row["State"] == files[1].strip():
            xlist.append(row["Exposure"])
            ylist.append(row["Uptake"])
    predictedvals = linreg(xlist, ylist, predict)
    
    #protein merge
    start, end, sequence = [], [], []
    for index, row in df.iterrows():
        start.append(row["Start"])
        end.append(row["End"])
        sequence.append(row["Sequence"])
    protein_merge(start, end, sequence)

    #writing all the data to a results file
    l = []
    for index, row in df.iterrows():
        l.append(float(row["Exposure"]))
    
    i = 0; j = []; k =0
    for index, row in df.iterrows():
        if row["State"] == files[1].strip(): 
            #j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[i], predictedvals[i], row["Uptake"], row["Exposure"]])
            j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[k], row["Uptake"], row["Exposure"], None])
            if l[k] < float(uptakepredict) and l[k+1] >= float(uptakepredict):
                j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[k], row["Uptake"], uptakepredict, predictedvals[i]])
                i+=1
            
            pd.DataFrame(j, index = None, columns = ["Start", "End", "Sequence","State", "MaxUptake", "Uptake", "Exposure", "Predicted Uptake" +f' ({uptakepredict})' ]).to_csv(files[0] + ".results.csv")
            k+=1


    dfres = pd.read_csv(files[0] +".results.csv")
    #vector
    modeltesting = []
    for index, row in dfres.iterrows():
        try: 
            len(row["Sequence"])
            #print(vector(row["Sequence"]))
        except TypeError:
            break
        modeltesting.append(vector(row["Sequence"]))
    
    #calculating drate
    hdx = []
    for index, row in dfres.iterrows():
        hdx.append(drate(row["Uptake"], row["MaxUptake"]))
    #print(hdx)
    
    #remove overlap
    start_list, end_list = [], []
    for index, row in dfres.iterrows():
        start_list.append(dfres["Start"])
        end_list.append(dfres["End"])
    removeOverlap(start_list[0], end_list[0])

    #print(modeltesting)
    #training(modeltesting, hdx)
    
    
 
