# Program Assignment2 of Course: Combinatorial Methods in Bioinformatics
# Author：生信1702 邵燕涵 2017317220205
# Software：PyCharm 2020.1.1 (Community Edition)
# Enviroment：python3.7

D = [[1, 2], [3, 6], [4, 2], [4, 5], [5, 6], [5, 7], [6, 6], [6, 4], [7, 10], [8, 11], [8, 6]]
w = 2
A = 6
two_hit = {}  # Define a dictionary to save (x,y):(x',y') that satisfy: x'-x = y'-y > w and x'- x< A

for item in range(len(D) - 1):  # A loop from the first item in set D to the penultimate item in set D, item = (x,y)
    # First of all, between the range of item~length(D) we need to find the (x',y') that from this (x',y') on, the
    # later are satisfy: x'-x > w by dichotomizing search algorithm
    low = item
    high = len(D) - 1
    site = 0  # To save the index of the (x,y) after first dichotomizing search
    while low <= high:
        mid = (low + high) // 2
        if D[mid][0] - D[item][0] > w:
            high = mid - 1
            site = mid
        else:
            low = mid + 1

    if site == 0:  # If there is no such (x',y') satisfy: x'-x > w, because the x is ordered so it means there is no
        # more (x',y') satisfy: x'-x > w, then we finish the whole algorithm and report
        break

    # Secondly, between the range of site~length(D) we need to find the (x',y') that from this (x',y') on, the later are
    # satisfy: x'-x < A by dichotomizing search algorithm
    low = site
    high = len(D) - 1
    site2 = 0  # To save the index of the (x',y') after second dichotomizing search
    while low <= high:
        mid = (low + high) // 2
        if D[mid][0] - D[item][0] < A:
            low = mid + 1
            site2 = mid
        else:
            high = mid - 1

    twohit = []  # Define a list to save all of the (x',y') that satisfied all of the conditions
    for i in range(site, site2 + 1):  # Finally, we've already find the range makes (x',y') satisfy: w < x'-x < A, so
        # now we just need to find the (x',y') satisfy the last condition: x'-x = y'-y
        if D[i][0] - D[item][0] == D[i][1] - D[item][1]:
            twohit.append(D[i])
    if twohit:  # Append the result from this turn into dictionary two_hit
        two_hit[str(D[item])] = twohit

for m, n in two_hit.items():  # Show all of the result pairs
    print(m, ':', n)
