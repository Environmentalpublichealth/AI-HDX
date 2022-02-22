"""
Created on Oct 6th, 2021
Author: Jiali Yu
Usage: remove low quality prediction based on the confidence score
"""

import sys
file_in = sys.argv[1]
file_out = sys.argv[2]
cutoff = 70 # percentage for alphafold and ratio for rosettafold
score_pos = -2 # the second last for alphafold and the last for rosettafold

with open(file_in) as f, open(file_out, "w") as out:
    for line in f:
        if "ATOM" in line:
            content = line.strip().split()
            if float(content[-2]) > cutoff:
                out.write(line)