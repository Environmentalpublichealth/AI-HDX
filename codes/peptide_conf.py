"""
Created on Feb 16, 2022
Author: Jiali Yu
Usage: calculate peptide confidence score for structure prediction
python peptide_conf.py <structure.pdb> <peptide fragment.csv>
"""

import sys
file_in = sys.argv[1]
file_tbl = sys.argv[2]
score_pos = -2 # the second last for alphafold and the last for rosettafold

def Average(lst):
    return sum(lst) / len(lst)

with open(file_in) as f1, open(file_tbl) as f2:
    scores = []
    for line in f1:
        if "ATOM" in line:
            content = line.strip().split()
            if content[2] == "CA":
                scores.append(float(content[score_pos]))

    for pep in f2:
        fragment = pep.strip().split(",")
        start = int(fragment[2])
        end = int(fragment[3])
        pep_score = Average(scores[start-1: end])
        print(pep_score)
