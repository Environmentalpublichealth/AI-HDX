"""
Created on July 21, 2021
Author: Jiali
Usage: calculate the D ratio using Maxuptake and uptake based on the sequence
Run:
python calculate_ptg.py <input table> <output table>
"""
import pandas as pd
import sys
in_file = sys.argv[1]
file_out = sys.argv[2]
# read file
data = pd.read_csv(in_file)
data.columns = ["start","end","sequence","uptake"]
## compute max uptake
mylist = []
for i in data["sequence"]:
  size = len(i)
  P = i.count("P")
  max_uptake = size - 1 - P
  mylist.append(max_uptake)
data["MaxUptake"] = mylist
# calculate exchange rate, considering 10% back-exchange
data["HDXrate"] = 100*(data["uptake"] / (data["MaxUptake"]*0.9))
data_sub = data[["start","end","sequence","HDXrate"]]
data_sub.to_csv(file_out, index=False, header=False)