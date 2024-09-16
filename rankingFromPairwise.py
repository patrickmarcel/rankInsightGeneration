import numpy as np
from statStuff import welch_ttest, permutation_test

# global variable
pairwiseComparison = []


def computeRanksForAll(pairwiseComparison, Sels):
    """
    For each item in adom will count how many times it has 'won' a pairwise comparison
    :param pairwiseComparison:
    :param Sels:
    :return:
    """
    # this is borda count style
    ranks = {}
    for s in Sels:
        ranks[s] = 0
    for left, right, result, _, _ in pairwiseComparison:
        if result == 1:
            ranks.update({left: ranks[left] + 1})
        if result == -1:
            ranks.update({right: ranks[right] + 1})

    #print("ranks:", ranks)
    return ranks


def balanced_rank_estimation(pairwise_comparisons, max_iterations=1000, tol=1e-6):
    """
    Perform Balanced Rank Estimation based on pairwise comparisons.

    Parameters:
    - pairwise_comparisons: list of tuples (item1, item2, result)
      where result is 1 if item1 > item2, -1 if item1 < item2, and 0 if item1 == item2.
    - max_iterations: maximum number of iterations for the algorithm (default is 1000).
    - tol: tolerance for convergence (default is 1e-6).

    Returns:
    - ranks: dictionary with items as keys and their estimated ranks as values.
    """
    # Extract unique items from the pairwise comparisons
    items = set()
    for item1, item2, result in pairwise_comparisons:
        items.add(item1)
        items.add(item2)

    items = list(items)
    n = len(items)
    item_index = {item: i for i, item in enumerate(items)}

    # Initialize ranks to zero
    ranks = np.zeros(n)

    # Create a matrix to keep track of comparisons and results
    comparison_matrix = np.zeros((n, n))
    for item1, item2, result in pairwise_comparisons:
        i, j = item_index[item1], item_index[item2]
        comparison_matrix[i, j] = result
        comparison_matrix[j, i] = -result

    for iteration in range(max_iterations):
        old_ranks = ranks.copy()

        for i in range(n):
            sum_comparisons = 0
            sum_ranks = 0
            for j in range(n):
                if comparison_matrix[i, j] != 0:
                    sum_comparisons += 1
                    sum_ranks += old_ranks[j] + comparison_matrix[i, j]

            if sum_comparisons > 0:
                ranks[i] = sum_ranks / sum_comparisons

        if np.linalg.norm(ranks - old_ranks, ord=1) < tol:
            break

    return {items[i]: ranks[i] for i in range(n)}




def findTuple(a, b, claireTab):
    for c in claireTab:
        if (c[0] == a and c[1] == b) or (c[1] == a and c[0] == b):
            return c


def compare(a, b, claireTab):
    cl = findTuple(a, b, claireTab)

    if tuple[2]:
        # print("Welch test can be used")
        #nbWelch = nbWelch + 1
        t_stat, p_value, conclusion = welch_ttest(cl[3], cl[4])

        comp = 0  # not significant
        if p_value < 0.05 and t_stat < 0:
            comp = -1
        if p_value < 0.05 and t_stat > 0:
            comp = 1
        pairwiseComparison.append((cl[0], cl[1], comp, t_stat, float(p_value)))
    else:
        # print("Permutation test is used")
        #nbPermut = nbPermut + 1
        observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(cl[3], cl[4])

        comp = 0  # not significant
        if p_value < 0.05 and observed_t_stat < 0:
            comp = -1
        if p_value < 0.05 and observed_t_stat > 0:
            comp = 1
        pairwiseComparison.append((cl[0], cl[1], comp, observed_t_stat, float(p_value)))
    result = False
    if comp == -1 or comp == 0:
        result = True
    return result, pairwiseComparison


def compare_cost(a, b, claireTab):
    """
    This function should be implemented to return the cost of comparing two objects a and b.
    """
    # Placeholder implementation; replace with the actual cost function.
    result = 0
    for c in claireTab:
        if (c[0] == a and c[1] == b) or (c[1] == a and c[0] == b):
            if c[2]:
                result = 0
            else:
                result = 1
    return result


def split_with_min_cost(arr, claireTab):
    """
    Splits the array into two parts with the minimum comparison cost.
    """
    n = len(arr)
    min_cost = float('inf')
    best_split = None

    for i in range(1, n):
        left = arr[:i]
        right = arr[i:]
        cost = sum(compare_cost(x, y, claireTab) for x in left for y in right)

        if cost < min_cost:
            min_cost = cost
            best_split = (left, right)

    return best_split


def merge(left, right, claireTab):
    """
    Merge two sorted lists into one sorted list.
    """
    sorted_list = []
    i = j = 0

    while i < len(left) and j < len(right):
        bool = compare(left[i], right[j], claireTab)  # True if left[i] <= right[j]
        if bool:  #left[i] <= right[j]
            sorted_list.append(left[i])
            i += 1
        else:
            sorted_list.append(right[j])
            j += 1

    sorted_list.extend(left[i:])
    sorted_list.extend(right[j:])

    return sorted_list


def merge_sort(arr, claireTab):
    """
    Modified merge sort that splits the array to minimize the cost of comparisons.
    """
    if len(arr) <= 1:
        return arr

    left, right = split_with_min_cost(arr, claireTab)
    sorted_left = merge_sort(left, claireTab)
    sorted_right = merge_sort(right, claireTab)

    return merge(sorted_left, sorted_right, claireTab)

# Example usage
#arr = [10, 20, 15, 5, 25, 30]
#sorted_arr = merge_sort(arr)
#print("Sorted array:", sorted_arr)
