"""
Created on August 25th, 2021
Author: Jiali
This program is to reformat the result tables after interpolate, extracting only one time and calculate the HDX rate for each peptide
Usage: python format_results.py <input data file> <output csv file>
"""

import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

with open(file_in) as f, open(file_out, "w") as out:
    header = f.readline()
    for line in f:
        content = line.strip().split(",")
        if content[8] != '':
            # print(content)
            MaxUptake = len(content[3]) - 1 - content[3].count("P")
            uptake = float(content[8])
            HDXrate = uptake/(MaxUptake * 0.8) * 100
            if HDXrate > 100:
                new_content = content[1:4] + ["100"]
            elif HDXrate < 0:
                continue
            else:
                new_content = content[1:4] + [str(HDXrate)]
            out.write(",".join(new_content)+"\n")