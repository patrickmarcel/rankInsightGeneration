from itertools import chain, combinations


def f_measure_first_k_keys(dict1, dict2, k):
    """
    Compute the F-measure score for the first k keys of two dictionaries.

    Args:
        dict1 (dict): The first dictionary.
        dict2 (dict): The second dictionary.
        k (int): Number of keys to consider from each dictionary.

    Returns:
        float: The F-measure score, a value between 0 and 1.
    """

    if k!=0:
        # Get the first k keys from both dictionaries as sets
        keys_set1 = set(list(dict1.keys())[:k])
        keys_set2 = set(list(dict2.keys())[:k])
    else:
        keys_set1 = set(list(dict1.keys()))
        keys_set2 = set(list(dict2.keys()))

    # Calculate precision, recall, and F-measure
    intersection = keys_set1.intersection(keys_set2)
    precision = len(intersection) / len(keys_set1) if keys_set1 else 0.0
    recall = len(intersection) / len(keys_set2) if keys_set2 else 0.0
    if precision + recall == 0:
        f_measure = 0.0
    else:
        f_measure = 2 * (precision * recall) / (precision + recall)

    return precision, recall, f_measure

def jaccard_score_first_k_keys(dict1, dict2, k):
    """
    Compute the Jaccard score for the first k keys of two dictionaries.

    Args:
        dict1 (dict): The first dictionary.
        dict2 (dict): The second dictionary.
        k (int): Number of keys to consider from each dictionary.

    Returns:
        float: The Jaccard score, a value between 0 and 1.
    """
    if k!=0:
        # Get the first k keys from both dictionaries as sets
        keys_set1 = set(list(dict1.keys())[:k])
        keys_set2 = set(list(dict2.keys())[:k])
    else:
        keys_set1 = set(list(dict1.keys()))
        keys_set2 = set(list(dict2.keys()))

    # Calculate the intersection and union of the sets
    intersection = keys_set1.intersection(keys_set2)
    union = keys_set1.union(keys_set2)

    # Compute the Jaccard score
    jaccard = len(intersection) / len(union) if union else 0.0
    return jaccard


def sort_dict_by_second_entry_desc(input_dict):
    """
    Sort a dictionary by the second entry of the list contained at each key in descending order.

    Args:
        input_dict (dict): A dictionary where each value is a list.

    Returns:
        dict: A new dictionary sorted by the second entry of the lists in descending order.
    """
    # Ensure the dictionary values are lists with at least two elements
    filtered_dict = {k: v for k, v in input_dict.items() if isinstance(v, list) and len(v) > 1}

    # Sort by the second entry of the lists in descending order
    sorted_items = sorted(filtered_dict.items(), key=lambda item: item[1][1], reverse=True)

    # Convert back to a dictionary
    sorted_dict = dict(sorted_items)
    return sorted_dict



def accumulate_numbers(input_list):
    """
    Accumulates the numbers in the input list such that each position i
    contains the sum of all previous numbers (including itself).

    :param input_list: List of numbers to be accumulated
    :return: A new list with accumulated sums
    """
    accumulated_list = []
    running_sum = 0
    for num in input_list:
        running_sum += num
        accumulated_list.append(running_sum)
    return accumulated_list

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
