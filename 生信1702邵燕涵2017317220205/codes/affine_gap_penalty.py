# Program Assignment of Course: Combinatorial Methods in Bioinformatics
# Author：生信1702 邵燕涵 2017317220205
# Software：PyCharm 2020.1.1 (Community Edition)
# Enviroment：python3.7
"""
This file is the algorithm to handle affine gap penalty, including two functions for input sequences and parameter file
reading, one function for back-tracking, three functions for DP banded-DP and X-drop and last the one is main function.
To testing this algorithm all you need to do is to fill in the path of your input file parameter file and three output
file in the main function area. After that, just run this file and you will get the result! For more details please
check the annotations below.
Notice: Because the input file: input2.txt and parameter file: parameter2.txt are protein seqence, and the score for
initiating a gap are not 0, so I using them as the testing data in main function. But if you want, you can also using
the nucleotide  seqence input1.txt, input3.txt, parameter1.txt and parameter3.txt as a testing data, but remember their
score for initiating a gap are 0, maybe you need to change it first, and before that, you should also remember to
change the file path and change the protein sequence reading function to nucleotide  input sequence reading function,
because there is a little difference between them.
"""
import numpy as np  # numpy package for matrix build
import time  # time package for running time calculate


# Read in input sequences file
def read_input_file(file_path):
    with open(file_path, 'r') as f:
        list1 = f.readlines()  # Read in file row by row
    seq1 = []  # Define two list to store each character of sequence
    seq2 = []
    for i in range(len(list1)):
        if list1[i].startswith('>seq1'):  # Read in each character of sequence1
            while list1[i + 1] != '\n':  # Just read the lines with sequence
                for j in list1[i + 1]:
                    if j != '\n' and j != ' ':  # For each line, ignoring Spaces and line breaks
                        seq1.append(j)  # Add each character to list
                i = i + 1
        if list1[i].startswith('>seq2'):
            for j in range(i + 1, len(list1)):
                for m in list1[j]:
                    if m != '\n' and m != ' ':
                        seq2.append(m)
    print('\nThe input sequences file %s was read successfully!\nseq1=%s\nseq2=%s\n' % (
        file_path.split('/')[-1], ''.join(seq1), ''.join(seq2)))  # If reading process is successful, then report
    return seq1, seq2  # Return two list of sequence


# read in nucleotide sequence's parameter file
def read_nucleotide_parameter_file(file_path):
    with open(file_path, 'r') as f:
        list1 = f.readlines()  # Read in file row by row
    threshold_X = 0  # Initialize seven variables to store the threshold X for X-drop, bandwidth B,score for initiating
    # a gap, score for each base insert/delete, the set of nucleotides and the similarity matrix
    bandwidth_B = 0
    gap_initiating_score = 0
    pergap_score = 0
    alphabet = []
    similarity_matrix = []
    for i in range(len(list1)):
        if list1[i].endswith('the threshold X for X-drop\n'):  # Read in the value of threshold X for X-drop
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    threshold_X = int(list1[i][0:j])
        if list1[i].endswith('bandwidth B\n'):  # Read in the value of bandwidth B
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    bandwidth_B = int(list1[i][0:j])
        if list1[i].endswith('score for initiating a gap\n'):  # Read in the score of initiating a gap
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    gap_initiating_score = int(list1[i][0:j])
        if list1[i].endswith('score for each base insert/delete\n'):  # Read in the score for each base insert/delete
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    pergap_score = int(list1[i][0:j])
        if list1[i].endswith('Below is the alphabet\n'):  # Read in the set of neclotides as a list
            alphabet = list1[i + 1].rsplit()
        if list1[i].endswith('Below is the similarity matrix for the alphabet\n'):  # Read in the similarity matrix as a list
            for j in range(i + 1, len(list1)):
                if list1[j] != '\n':
                    similarity_matrix.append(list1[j].rsplit())
    for i in range(len(similarity_matrix)):
        for j in range(len(similarity_matrix[0])):  # Transfrom the substitution matrix list to matrix
            similarity_matrix[i][j] = int(similarity_matrix[i][j])
    print(
        'The parameter file %s was read successfully!\nthe threshold X for X-drop is:%d\nbandwidth B is:%d\nscore for '
        'initiating a gap is:%d\nscore for each base '
        'insert/delete is:%d\nalphabet:%s\nsimilarity matrix of alphabet:\n%s\n' % (
            file_path.split('/')[-1], threshold_X, bandwidth_B, gap_initiating_score, pergap_score, alphabet,
            np.array(similarity_matrix)))  # If reading process is successful, then report
    return threshold_X, bandwidth_B, gap_initiating_score, pergap_score, alphabet, similarity_matrix  # return seven
    # Variables we needed


# Read in protein sequence's parameter file
def read_protein_parameter_file(file_path):
    with open(file_path, 'r') as f:
        list1 = f.readlines()  # Read in file row by row
    threshold_X = 0  # Initialize seven variables to store the threshold X for X-drop, bandwidth B,score for initiating
    # a gap, score for each base insert/delete, the set of amino acids and the substitution matrix
    bandwidth_B = 0
    gap_initiating_score = 0
    pergap_score = 0
    alphabet = []
    matrix_title = ()
    substitution_matrix = []
    for i in range(len(list1)):  # Read in the value of threshold X for X-drop
        if list1[i].endswith('the threshold X for X-drop\n'):
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    threshold_X = int(list1[i][0:j])
        if list1[i].endswith('bandwidth B\n'):  # Read in the value of bandwidth B
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    bandwidth_B = int(list1[i][0:j])
        if list1[i].endswith('score for initiating a gap\n'):  # Read in the score of initiating a gap
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    gap_initiating_score = int(list1[i][0:j])
        if list1[i].endswith('score for each base insert/delete\n'):  # Read in the score for each base insert/delete
            for j in range(len(list1[i])):
                if list1[i][j] == ';':
                    pergap_score = int(list1[i][0:j])
        if list1[i].endswith('Below is the set of amino acids\n'):  # Read in the set of amino acids as a list
            alphabet = list1[i + 1].rsplit()
        if list1[i].endswith('substitution matrix\n'):  # Read in the substitution matrix as a list
            matrix_title = ' '.join(list1[i].rsplit()[1:])
            for j in range(i + 1, len(list1)):
                if list1[j] != '\n':
                    substitution_matrix.append(list1[j].rsplit())
    for i in range(len(substitution_matrix)):  # Transfrom the substitution matrix list to matrix
        for j in range(len(substitution_matrix[0])):
            substitution_matrix[i][j] = int(substitution_matrix[i][j])
    print(
        'The parameter file %s was read successfully!\nthe threshold X for X-drop is:%d\nbandwidth B is:%d\nscore for '
        'initiating a gap is:%d\nscore for each base '
        'insert/delete is:%d\nthe set of amino acids:%s\n%s:\n%s\n' % (
            file_path.split('/')[-1], threshold_X, bandwidth_B, gap_initiating_score, pergap_score, alphabet,
            matrix_title, np.array(substitution_matrix)))  # If reading process is successful, then report
    return threshold_X, bandwidth_B, gap_initiating_score, pergap_score, alphabet, substitution_matrix  # return seven
    # Variables we needed


# back-tracking of Smith-Waterman algorithm with affine gap penalty
def back_tracking_of_affine_gap_penalty(seq1, seq2, v_table, e_table, f_table, alphabet, substitution_matrix,
                                        pergap_score):
    output_seq1 = []  # Define two list to store the aligned sequence from back-tracking
    output_seq2 = []
    for m in range(1, len(v_table)):
        for n in range(1, len(v_table[0])):  # Because this algorithm is based on Smith-Waterman algorithm, so we need
            # to do back-tracking from the max number of the V table
            if v_table[m][n] == v_table.max():  # Find the i and j of max number's entry
                i = m
                j = n
    while v_table[i][j] != 0:  # If the entry is not the first entry then do circulation
        for m in range(len(alphabet)):  # First of all, we need to know the character of the entry now, so we find it
            # and give it to m and n
            if seq1[i - 1] == alphabet[m]:
                break
        for n in range(len(alphabet)):
            if seq2[j - 1] == alphabet[n]:
                break
        if v_table[i][j] == v_table[i - 1][j - 1] + substitution_matrix[m][n]:  # Case1: S[i] aligns with T[j]
            output_seq1.append(alphabet[m])
            output_seq2.append(alphabet[n])
            i = i - 1
            j = j - 1
        elif v_table[i][j] == e_table[i][j]:  # Case2: Instertion
            while e_table[i][j] == e_table[i][j - 1] + pergap_score:  # If there is a gap extention case then continue
                output_seq1.append('-')
                output_seq2.append(alphabet[n])
                j = j - 1
                for m in range(len(alphabet)):  # Because the j was changed above, so we need to refind the character
                    # of the entry now
                    if seq1[i - 1] == alphabet[m]:
                        break
                for n in range(len(alphabet)):
                    if seq2[j - 1] == alphabet[n]:
                        break
            output_seq1.append('-')  # The gap extention is finish
            output_seq2.append(alphabet[n])
            j = j - 1
        elif v_table[i][j] == f_table[i][j]:  # Case3: Deletion
            while f_table[i][j] == f_table[i - 1][j] + pergap_score:  # If there is a gap extention case then continue
                output_seq1.append(alphabet[m])
                output_seq2.append('-')
                i = i - 1
                for m in range(len(alphabet)):  # Because the j was changed above, so we need to refind the character
                    # of the entry now
                    if seq1[i - 1] == alphabet[m]:
                        break
                for n in range(len(alphabet)):
                    if seq2[j - 1] == alphabet[n]:
                        break
            output_seq1.append(alphabet[m])  # The gap extention is finish
            output_seq2.append('-')
            i = i - 1
    output_seq1.reverse()  # Reverse the back-tracking list so we can see it in a normal way
    output_seq2.reverse()
    output_seq1 = ''.join(output_seq1)  # Paste the list to a character string
    output_seq2 = ''.join(output_seq2)
    return output_seq1, output_seq2  # Return two strings


# Smith-Waterman algorithm with affine gap penalty
def Smith_Waterman_algorithm_with_affine_gap_penalty(seq1, seq2, gap_initiating_score, alphabet, substitution_matrix,
                                                     pergap_score, file_path):
    start_time = time.time()  # Record the start time of the program
    v_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # Initialize the V table with zero
    e_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # Initialize the E table with zero
    for i in range(len(e_table)):  # For E table, initialize the first column with negative infinity(-999) and the
        # first row by gap penalty model
        e_table[i][0] = -999
    for i in range(1, len(e_table[0])):
        e_table[0][i] = gap_initiating_score + i * pergap_score
    f_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # Initialize the F table with zero
    for i in range(len(f_table[0])):  # For F table, initialize the first row with negative infinity(-999) and the
        # first column by gap penalty model
        f_table[0][i] = -999
    for i in range(1, len(f_table)):
        f_table[i][0] = gap_initiating_score + i * pergap_score
    count = 0  # Record the counts of fill in entries
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
            ematch = e_table[i][j - 1] + pergap_score  # Calculate the two cases of formula:
            # E(i, j) = max { E(i, j-1)–s, V(i, j-1)–h–s } and give them to ematch and einsertion
            einsertion = v_table[i][j - 1] + gap_initiating_score + pergap_score
            fmatch = f_table[i - 1][j] + pergap_score  # Calculate the two cases of formula:
            # F(i, j) = max { F(i-1, j)–s, V(i-1, j)–h–s } and give them to fmatch and fdeletion
            fdeletion = v_table[i - 1][j] + gap_initiating_score + pergap_score
            e_table[i][j] = max(ematch, einsertion)  # The V[i][j] value is based on the E[i][j] and F[i][j], so we need
            # to calculate them first
            f_table[i][j] = max(fmatch, fdeletion)
            match = v_table[i - 1][j - 1] + substitution_matrix[n][m]  # Match case without gap
            v_table[i][j] = max(match, e_table[i][j], f_table[i][j], 0)  # Fill in the V[i][j] basing on four cases
            count = count + 1  # Finish a circle and count plus one
    count = count + len(v_table) + len(v_table[0]) - 1  # Finish all circulation, and the sum of the entries been filled
    # in equals to the counts in circulation plus the counts of initialize
    output_seq1, output_seq2 = back_tracking_of_affine_gap_penalty(seq1, seq2, v_table, e_table, f_table, alphabet,
                                                                   substitution_matrix, pergap_score)  # Calling
    # back-tracking function
    end_time = time.time()  # Record the end time of the program
    with open(file_path, 'w') as f:  # Write the score, entries, running time and back-tracking sequence to file
        f.write('score={0}\nentries={1}\nrunning time={2}/s\n\n>seq1\n{3}\n\n>seq2\n{4}'.format(int(np.max(v_table)),
                                             count, round(end_time - start_time, 4), output_seq1, output_seq2))
    print('The result of Smith-Waterman algorithm with affine gap penalty:\nscore=%d\nentries=%d\nrunning '
          'time=%f/s\noutput_1:%s\noutput_2:%s\nThe result was already written into the file under path:%s !\n' % (
            np.max(v_table), count, end_time - start_time, output_seq1, output_seq2, file_path))  # report the program's
    # running result


# banded Smith-Waterman algorithm with affine gap penalty
def banded_Smith_Waterman_algorithm_with_affine_gap_penalty(seq1, seq2, bandwidth_B, gap_initiating_score, alphabet,
                                                            substitution_matrix, pergap_score, file_path):
    start_time = time.time()  # Record the start time of the program
    v_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # Initialize the V table with zero
    x = []
    for i in range(bandwidth_B + 1):  # Get the set to direct the algorithm fill in the table within banded space
        x.append(i)
        x.append(i * -1)
    x.remove(0)  # If the bandwidth is 3, the x list must be [0,1,2,3,-1,-2,-3], which shows the diagonal section
    e_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # initialize the E table with zero
    for i in range(bandwidth_B + 1):  # For E table, initialize the first column of banded section with negative
        # infinity(-999) and the first row of banded section by gap penalty model
        e_table[i][0] = -999
    for i in range(1, bandwidth_B + 1):
        e_table[0][i] = gap_initiating_score + i * pergap_score
    f_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # initialize the F table with zero
    for i in range(bandwidth_B + 1):  # For F table, initialize the first row of banded section with negative
        # infinity(-999) and the first column of banded section by gap penalty model
        f_table[0][i] = -999
    for i in range(1, bandwidth_B + 1):
        f_table[i][0] = gap_initiating_score + i * pergap_score
    count = 0  # Record the counts of fill in entries
    for i in range(1, len(v_table)):
        for j in range(1, len(v_table[0])):
            if i - j in x:  # For banded Smith-Waterman dynamic programing, we just need to fill in the entries of the
                # diagonal section based on bandwidth, which follow the formula: i - j in x
                for m in range(len(alphabet)):  # First of all, we need to know the character of the entry now, so we
                    # find it and give it to m and n
                    if seq1[i - 1] == alphabet[m]:
                        break
                for n in range(len(alphabet)):
                    if seq2[j - 1] == alphabet[n]:
                        break
                ematch = e_table[i][j - 1] + pergap_score  # Calculate the two cases of formula:
                # E(i, j) = max { E(i, j-1)–s, V(i, j-1)–h–s } and give them to ematch and einsertion
                einsertion = v_table[i][j - 1] + gap_initiating_score + pergap_score
                fmatch = f_table[i - 1][j] + pergap_score  # Calculate the two cases of formula:
                # F(i, j) = max { F(i-1, j)–s, V(i-1, j)–h–s } and give them to fmatch and fdeletion
                fdeletion = v_table[i - 1][j] + gap_initiating_score + pergap_score
                e_table[i][j] = max(ematch, einsertion)  # The V[i][j] value is based on the E[i][j] and F[i][j], so we
                # need to calculate them first
                f_table[i][j] = max(fmatch, fdeletion)
                match = v_table[i - 1][j - 1] + substitution_matrix[n][m]  # Match case without gap
                v_table[i][j] = max(match, e_table[i][j], f_table[i][j], 0)  # Fill in the V[i][j] base on four cases
                count = count + 1  # Finish one circle and count plus one
    count = count + bandwidth_B * 2 + 1  # Finish all circulation, and the sum of the entries been filled in equals to
    # the counts in circulation plus the counts of initialize
    output_seq1, output_seq2 = back_tracking_of_affine_gap_penalty(seq1, seq2, v_table, e_table, f_table, alphabet,
                                                                   substitution_matrix, pergap_score)  # Calling
    # back-tracking function
    end_time = time.time()  # Record the end time of the program
    with open(file_path, 'w') as f:  # Write the score, entries, running time and back-tracking sequence to file
        f.write('score={0}\nentries={1}\nrunning time={2}/s\n\n>seq1\n{3}\n\n>seq2\n{4}'.format(int(np.max(v_table)),
                                                count, round(end_time - start_time, 4), output_seq1, output_seq2))
    print(
        'The result of banded Smith-Waterman algorithm with affine gap penalty:\nscore=%d\nentries=%d\nrunning '
        'time=%f/s\noutput_1:%s\noutput_2:%s\nThe result was already written into the file under path:%s !\n' % (
            np.max(v_table), count, end_time - start_time, output_seq1, output_seq2, file_path))  # report the program's
    # running result


# X-drop algorithm with affine gap penalty
def x_drop_algorithm_with_affine_gap_penalty(seq1, seq2, gap_initiating_score, threshold_X, pergap_score, alphabet,
                                             substitution_matrix, file_path):
    start_time = time.time()  # Record the start time of the program
    v_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # initialize the V table with zero
    e_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # initialize the E table with zero
    for i in range(len(e_table)):  # For E table, initialize the first column with negative infinity(-999) and the
        # first row by gap penalty model
        e_table[i][0] = -999
    for i in range(1, len(e_table[0])):
        e_table[0][i] = gap_initiating_score + i * pergap_score
    f_table = np.zeros((len(seq1) + 1, len(seq2) + 1))  # initialize the F table with zero
    for i in range(len(f_table[0])):  # For F table, initialize the first row with negative infinity(-999) and the
        # first column by gap penalty model
        f_table[0][i] = -999
    for i in range(1, len(f_table)):
        f_table[i][0] = gap_initiating_score + i * pergap_score
    diagonal_d = len(v_table) + len(v_table[0]) - 1  # The sum of the antidiagonal of this matrix
    count = 0  # Record the counts of fill in entries
    r = 0  # To record the max number of the matrix dynamically
    r_change = [0]  # Record each r after one circle, so we can use it to judge the consistance of the entry in the
    # former antidiagonal
    for d in range(1, diagonal_d):  # Starting from the antidiagonal one
        consistent = 0  # Recording the entries we fill in of each antidiagonal
        for i in range(len(v_table)):
            for j in range(len(v_table[0])):
                if i + j == d:  # The condition to find the entries of each antidiagonal
                    if i == 0:  # If this entry is from the first row
                        if v_table[i][j - 1] >= r_change[d - 1] - threshold_X:  # If this entry's former one is
                        # consistent, then fill in this entry basing on E table beacause we can ensure it's a insertion
                            v_table[i][j] = e_table[i][j]
                            count = count + 1  # total fill in number plus one
                            consistent = consistent + 1  # total fill in number of present antidiagonal plus one
                        else:
                            v_table[i][j] = -999  # If this entry's former one is not consistent, then we can say this
                            # is a blank space, so we fill in a number -999
                    if j == 0:  # If this entry is from the first row
                        if v_table[i - 1][j] >= r_change[d - 1] - threshold_X:  # If this entry's former one is
                        # consistent, then fill in this entry basing on F table beacause we can ensure it's a insertion
                            v_table[i][j] = f_table[i][j]
                            count = count + 1  # total fill in number plus one
                            consistent = consistent + 1  # total fill in number of present antidiagonal plus one
                        else:
                            v_table[i][j] = -999  # If this entry's former one is not consistent, then we can say this
                            # is a blank space, so we fill in a number -999
                    if i != 0 and j != 0:  # If this entry is not from the first row and column
                        if True in (v_table[i - 1][j - 1] >= r_change[d - 1] - threshold_X,
                                    v_table[i - 1][j] >= r_change[d - 1] - threshold_X,
                                    v_table[i][j - 1] >= r_change[d - 1] - threshold_X):  # If one of this entry's former
                            # one including V[i-1][j-1], V[i-1][j] and V[i][j-1] is consistent, then fill in this entry
                            # with the same ways like DP and banded-DP algorithm
                            for m in range(len(alphabet)):  # First of all, we need to know the character of the entry
                                # now, so we find it and give it to m and n
                                if seq1[i - 1] == alphabet[m]:
                                    break
                            for n in range(len(alphabet)):
                                if seq2[j - 1] == alphabet[n]:
                                    break
                            ematch = e_table[i][j - 1] + pergap_score  # Calculate the two cases of formula:
                            # E(i, j) = max { E(i, j-1)–s, V(i, j-1)–h–s } and give them to ematch and einsertion
                            einsertion = v_table[i][j - 1] + gap_initiating_score + pergap_score
                            fmatch = f_table[i - 1][j] + pergap_score  # Calculate the two cases of formula:
                            # F(i, j) = max { F(i-1, j)–s, V(i-1, j)–h–s } and give them to fmatch and fdeletion
                            fdeletion = v_table[i - 1][j] + gap_initiating_score + pergap_score
                            e_table[i][j] = max(ematch, einsertion)  # The V[i][j] value is based on the E[i][j] and
                            # F[i][j], so we need to calculate them first
                            f_table[i][j] = max(fmatch, fdeletion)
                            match = v_table[i - 1][j - 1] + substitution_matrix[n][m] # Match case without gap
                            v_table[i][j] = max(match, e_table[i][j], f_table[i][j])  # Fill in the V[i][j] base on
                            # three cases
                            count = count + 1  # total fill in number plus one
                            consistent = consistent + 1  # total fill in number of present antidiagonal plus one
                        else:
                            v_table[i][j] = -999  # If this entry's former one including V[i-1][j-1], V[i-1][j] and
                            # V[i][j-1] are not consistent, then we can say this is a blank space, so we fill in a
                            # number -999
        r_change.append(r)  # After one circle, we need to record the r value so that we can use it to judge the
        # consistence of the former one antidiagonal
        if np.max(v_table) > r:  # After we recorded the r value, we need to check all of the entries we've just filled
            # in and if we find a bigger one, then we need to update it
            r = np.max(v_table)
        if consistent == 0:  # If the total fill in number of present antidiagonal is zero, then we can say there is no
            # consistent entries of the former one antidiagonal, so we just break the whole circulation
            break
    count = count + 1  # Finish all circulation, and the sum of the entries been filled in equals to the counts in
    # circulation plus the first entry cause we initialize it as a consistent entry
    output_seq1, output_seq2 = back_tracking_of_affine_gap_penalty(seq1, seq2, v_table, e_table, f_table, alphabet,
                                                                   substitution_matrix, pergap_score)  # Calling
    # back-tracking function
    end_time = time.time()  # Record the end time of the program
    with open(file_path, 'w') as f:  # Write the score, entries, running time and back-tracking sequence to file
        f.write('score={0}\nentries={1}\nrunning time={2}/s\n\n>seq1\n{3}\n\n>seq2\n{4}'.format(int(np.max(v_table)),
                                                    count, round(end_time - start_time, 4), output_seq1, output_seq2))
    print(
        'The result of X-drop algorithm with affine gap penalty:\nscore=%d\nentries=%d\nrunning '
        'time=%f/s\noutput_1:%s\noutput_2:%s\nThe result was already written into the file under path:%s !\n' % (
            np.max(v_table), count, end_time - start_time, output_seq1, output_seq2, file_path))  # report the program's
    # running result


def main():  # Main function, using function packaging, to ensure the running memory savings, and ensure the normal use
    # of variables
    # When testing this Algorithm you need to fill in yourselves file path below, here we using input2 as an example
    input_file = '/Users/macbookair/Downloads/ass1/input2.txt'
    parameter_file = '/Users/macbookair/Downloads/ass1/parameter2.txt'
    DP_outfile = '/Users/macbookair/Downloads/out2_DP.txt'
    band_outfile = '/Users/macbookair/Downloads/out2_band.txt'
    drop_outfile = '/Users/macbookair/Downloads/out2_drop.txt'

    # Functions calling
    seq1, seq2 = read_input_file(input_file)
    threshold_X, bandwidth_B, gap_initiating_score, pergap_score, alphabet, substitution_matrix = read_protein_parameter_file(
        parameter_file)
    Smith_Waterman_algorithm_with_affine_gap_penalty(seq1, seq2, gap_initiating_score, alphabet, substitution_matrix,
                                                     pergap_score, DP_outfile)
    banded_Smith_Waterman_algorithm_with_affine_gap_penalty(seq1, seq2, bandwidth_B, gap_initiating_score, alphabet,
                                                            substitution_matrix, pergap_score, band_outfile)
    x_drop_algorithm_with_affine_gap_penalty(seq1, seq2, gap_initiating_score, threshold_X, pergap_score, alphabet,
                                             substitution_matrix, drop_outfile)


if __name__ == '__main__':
    main()
