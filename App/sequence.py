# >sp|P00709|LALBA_HUMAN
# MRFFVPLFLVGILFPAILAKQFTKCELSQLLKDIDGYGGIALPELICTMFHTSGYDTQAIVENNESTEYGLFQISNKLWCKSSQVPQSRNICDISCDKFLDDDITDDIMCAKKILDIKGIDYWLAHKALCTEKLEQWLCEKL
'''
The function should 
1. determine sequence length, 0 - 143
2. randomly cut the sequence into fragments with 5 - 20 letters 
 0,3, MRFF
 4,10, VPLFL
 11,20, VGILFPAILA
 
 '''
import random
sequence = list("MRFFVPLFLVGILFPAILAKQFTKCELSQLLKDIDGYGGIALPELICTMFHTSGYDTQAIVENNESTEYGLFQISNKLWCKSSQVPQSRNICDISCDKFLDDDITDDIMCAKKILDIKGIDYWLAHKALCTEKLEQWLCEKL")
sequence_len = len(sequence)

# function that will split the sequence into fragments
def split_list(lst):
    # makes a list to hold the fragments
    chunks = []
    i=0
    # cuts the sequence into fragments of a random size from 5-21 
    # and appends those fragments to the chunks list
    while i < len(lst):
        chunk_size = random.randint(5, 21)
        chunk = lst[i:i + chunk_size]
        chunks.append(chunk)
        i = i + chunk_size


    return chunks

chunks = split_list(sequence)

print("Chunks:", chunks)