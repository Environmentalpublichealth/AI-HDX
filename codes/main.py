import pandas as pd
from data_approximations import linreg
from protein_merge import protein_merge

manager = open("manager.txt", "r")

for row in manager: 
    
    files = row.split(": ")
    uptakepredict = input(files[0] +"Exposure prediction: ")
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
    for index, row in df.iterrows():
        
        if row["State"] == files[1].strip():
            j.append([row["Start"], row["End"], row["Sequence"], row["State"], mut[i], row["Uptake"], row["Exposure"]])
        
            pd.DataFrame(j, index = None, columns = ["Start", "End", "Sequence","State", "Maxuptake", "Uptake", "Exposure"]).to_csv(files[0] + "results.csv")
        i+=1
        
        
 
