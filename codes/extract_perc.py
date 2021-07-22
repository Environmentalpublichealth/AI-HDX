"""
Create on July 12th, 2021
Author: Jiali
Usage: This takes a HDX-MS result table with %D and takes 60 min time point and output into
a table with only sequence and %D
"""

import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

exposure = "3000"
with open(file_in) as f, open (file_out,"w") as out:
    header = f.readline().split(",")
    seq_idx = header.index("Sequence")
    time_idx = header.index(exposure)
    print(header)
    for line in f:
        content = line.split(",")
        out.write(content[seq_idx]+","+content[time_idx]+"\n")
        

