"""
Created on Feb 28th, 2022
Author: Jiali
For a given sequence and positions, subset it step by step into single residue
"""

# read input
seq = input('Sequence need to subset: ')
start = int(input('start position: '))
end = int(input('end position: '))

# from left to right
for i in range(1,len(seq)):
    res = seq[:i]
    new_start = start
    new_end = start + i - 1 
    print(new_start, new_end, res, sep = "\t")

# from right to left
for i in range(0,len(seq)):
    res = seq[i:]
    new_end = end
    new_start = start + i
    print(new_start, new_end, res, sep = "\t")