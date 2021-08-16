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
    if "Maxuptake" not in df[0:0]:
        for index, row in df.iterrows():
            mut.append(len(str((row["Sequence"]))) - 1 - str(row["Sequence"]).count("P"))
    else:
        for index, row in df.iterrows():
            mut.append(row["Maxuptake"])

    #uptake prediction
    xlist, ylist, predict, predictedvals = [], [], float(uptakepredict), []
    for index, row in df.iterrows():
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
    i = 0; j = []
    '''for index, row in df.iterrows():
        j = float(row["Exposure"])
        break
    while jc !=2:
        for index, row in df.iterrows(): #0, 1, 2, 3, 0, 1, 2, 3
            if float(row["Exposure"]) == j:
                jc+=1
            loop+=1
    loop = loop - 1
    flag = True; cc = 0
    while flag != False:
        for index, row in df.iterrows():
            cc+=1
            if float(row["Exposure"]) < uptakepredict:
                flag = True
            if float(row["Exposure"]) >= uptakepredict:
                flag = False '''
    for index, row in df.iterrows():
        if row["State"] == files[1].strip():
            
            #j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[i], predictedvals[i], row["Uptake"], row["Exposure"]])
            if row["Exposure"] == uptakepredict:
                j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[i], row["Uptake"], row["Exposure"], predictedvals[i]])
            else:
                j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[i], row["Uptake"], row["Exposure"], row["Uptake"]])
            pd.DataFrame(j, index = None, columns = ["Start", "End", "Sequence","State", "Maxuptake", "Uptake", "Exposure", "Predicted Maxuptake"]).to_csv(files[0] + "results.csv")
        i+=1


    dfres = pd.read_csv(files[0] +"results.csv")
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
        hdx.append(drate(row["Uptake"], row["Maxuptake"]))
    #print(hdx)
    
    #remove overlap
    start_list, end_list = [], []
    for index, row in dfres.iterrows():
        start_list.append(dfres["Start"])
        end_list.append(dfres["End"])
    removeOverlap(start_list[0], end_list[0])

    #print(modeltesting)
    training(modeltesting, hdx)
    
    
 
