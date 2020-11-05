# Program Assignment2 of Course: Combinatorial Methods in Bioinformatics
# Author：生信1702 邵燕涵 2017317220205
# Software：PyCharm 2020.1.1 (Community Edition)
# Enviroment：python3.7

import numpy as np

seq1 = ['t', 'g', 'a', 'c', 'a', 'a', 't', 'c', 'c', 'c']
seq2 = ['t', 'g', 'a', 'g', 'c', 'a', 't', 'g', 'g', 't']
alphabet = ['a', 'c', 'g', 't']
similarity_matrix = [[2, -1, -1, -1], [-1, 2, -1, -1], [-1, -1, 2, -1], [-1, -1, -1, 2]]
pergap_score = -1


def Needleman_Wunsch_algorithm_with_linear_gap_penalty(seq1, seq2, alphabet, similarity_matrix, pergap_score):
    v_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # initialize the V table with zero
    for i in range(1, len(v_table)):
        v_table[i, 0] = pergap_score + v_table[i - 1, 0]
    for i in range(1, len(v_table[0])):
        v_table[0, i] = pergap_score + v_table[0, i - 1]
    for i in range(1, len(v_table)):
        for j in range(1, len(v_table[0])):  # For standard Smith-Waterman dynamic programing, we need to fill in
            # the whole table
            for m in range(len(alphabet)):  # First of all, we need to know the character of the entry now, so we find
                # it and give it to m and n
                if seq1[i - 1] == alphabet[m]:
                    break
            for n in range(len(alphabet)):
                if seq2[j - 1] == alphabet[n]:
                    break
            match = v_table[i - 1][j - 1] + similarity_matrix[n][m]  # Calculate the match case
            deletion = v_table[i - 1][j] + pergap_score  # Calculate the insertion case
            insertion = v_table[i][j - 1] + pergap_score  # Calculate the deletion case
            v_table[i][j] = max(match, deletion, insertion)  # Fill in the V[i][j] basing on four cases
    return v_table


# Run Needleman-Wunsch algorithm and get the V table
v_table = Needleman_Wunsch_algorithm_with_linear_gap_penalty(seq1, seq2, alphabet, similarity_matrix, pergap_score)

# i,j is the position of last entry
i = len(v_table) - 1
j = len(v_table[0]) - 1
total_ways = [[i, j]]  # This list is used to dynamically save the new road's position and refresh the former road's
# position when running back-tracking
count = 0  # Count the number of all optimal alignments
while count != len(total_ways):  # If the count is not equals to the list's length, it means there are still some
    # alignments not finish the back-tracking
    for item in range(len(total_ways)):  # For each of the position in list, we need to do refreshing or appending until
        # the back-tracking finish
        ways = []  # To save the entry's source, there have three cases: just one source, two sources and three sources
        i = total_ways[item][0]  # Starting from the last entry
        j = total_ways[item][1]
        if i == j == 0:  # If there is a [0,0] appeared in list, it means there is one alignment finish the
            # back-tracking
            count += 1
            continue
        for m in range(len(alphabet)):  # First of all, we need to know the character of the entry now, so we find it
            # and give it to m and n
            if seq1[i - 1] == alphabet[m]:
                break
        for n in range(len(alphabet)):
            if seq2[j - 1] == alphabet[n]:
                break
        if v_table[i][j] == v_table[i - 1][j - 1] + similarity_matrix[m][n]:  # S[i] aligns with T[j]
            ways.append([i - 1, j - 1])
        if v_table[i][j] == v_table[i - 1][j] + pergap_score:  # Deletion
            ways.append([i - 1, j])
        if v_table[i][j] == v_table[i][j - 1] + pergap_score:  # Instertion
            ways.append([i, j - 1])
        if len(ways) == 1:  # If the length of ways(list) is 1, it means case1: This entry just one source
            total_ways[item] = ways[0]  # We just need to refesh the position to this source's position
        elif len(ways) > 1:  # Else there is case2 or case3, it means this entry have two or three sources
            total_ways[item] = ways[0]  # Refesh the position to the first source's position
            for s in ways[1:]:  # And then append new sources' position into the list(total_ways)
                total_ways.append(s)
print('The number of all optimal global alignments is:', count)
