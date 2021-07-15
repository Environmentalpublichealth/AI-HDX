import csv
import pandas as pd

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


manager.close()
data.close()
res.close()