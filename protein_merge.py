def protein_merge(start, end, sequence):
    
    def protein_finder(merged_protein, start, end, sequence, position): #recursive function
        #if end[position] > start[len(start)-1]: #this deals with the end of the sequence
        #    return merged_protein

        for i in range(position, len(start)):
            if start[i] - end[position] > 0 and start[i] - end[position] <= 5:
                merged_protein += sequence[i]
                return protein_finder(merged_protein, start, end, sequence, i)
            elif start[i] - end[position] > 5:
                break
        return merged_protein
    
    position = 0 #position is the current position of the list
    for i in range(len(start)): #it iterates down the start list, and breaks the for loop once a positive start value is reached
        if start[i] > 0:
            position = i
            break
        
    merged_protein = sequence[position]#then it takes the first protein at the position
    return protein_finder(merged_protein, start, end, sequence, position)

