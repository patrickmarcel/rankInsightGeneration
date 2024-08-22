from itertools import chain, combinations



def jaccard_similarity(set1, set2):
    """
    Compute the Jaccard similarity between two sets.

    Parameters:
    - set1: first set.
    - set2: second set.

    Returns:
    - Jaccard similarity coefficient.
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union



def powerset(s):
    """
    Generates the powerset of a given set s.

    :param s: The input set
    :return: A list of subsets representing the powerset
    """
    s = list(s)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))

def countNonZeros(pairwisecomp):
    nb=0
    for c in pairwisecomp:
        if c[2] != 0:
            nb=nb+1
    return nb

def checkConsistency(pairwisecomp):
    nb=0
    for c in pairwisecomp:
        for d in pairwisecomp:
            if c[0] == d[1] and c[1] == d[0]:
                if c[2] == d[2] and c[2] != 0:
                    nb=nb+1
    return nb

def listComp(l1, l2):
    nb=0
    tabDiff=[]
    for i in range(len(l1)):
        tabDiff.append(abs(l1[i]-l2[i]))
        if l1[i] != l2[i]:
            nb=nb+1
    return nb, tabDiff

def nbTrueInList(list):
    nb=0
    for l in list:
        if l == True:
            nb=nb+1
    return nb
