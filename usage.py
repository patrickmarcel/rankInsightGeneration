import numpy as np
import main
import statStuff
import plotStuff

def ExampleUsage():

# Example usage for skewness

    # Sample data
    data = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5])
    # Compute skewness
    skewness = main.compute_skewness(data)
    # Print result
    print(f"Skewness of the data: {skewness}")

# Example usage for Welch

    # Sample data
    group1 = np.array([23, 45, 56, 67, 78, 89, 90, 56, 54, 67])
    group2 = np.array([35, 56, 76, 78, 45, 34, 23, 45, 67, 88])

    # Perform Welch's t-test
    t_stat, p_value, conclusion = main.welch_ttest(group1, group2)

    # Print results
    print(f"t-statistic: {t_stat}")
    print(f"p-value: {p_value}")
    print(conclusion)

# Example usage for Brown Forsythe
    group1 = [23, 20, 25, 22, 21, 20, 23, 24]
    group2 = [22, 21, 24, 21, 20, 19, 21, 23]
    group3 = [26, 24, 27, 25, 26, 27, 28, 29]

    W, p_value = main.brown_forsythe_test(group1, group2, group3)
    print(f"Brown-Forsythe test statistic: {W}")
    print(f"P-value: {p_value}")

# Example usage for permutation test:
    np.random.seed(0)
    x1 = np.random.normal(0, 1, 30)
    x2 = np.random.normal(0.5, 1.5, 35)

    observed_t_stat, p_value, permuted_t_stats = main.permutation_test(x1, x2)
    print(f"Observed Welch's t-statistic: {observed_t_stat}")
    print(f"P-value: {p_value}")

    # Example usage of benjamini_hochberg:
    p_values = [0.01, 0.04, 0.03, 0.002, 0.05, 0.001, 0.03, 0.04]
    alpha = 0.05
    rejected, corrected_p_values = main.benjamini_hochberg(p_values, alpha)

    print("Rejected hypotheses:", rejected)
    print("Corrected p-values:", corrected_p_values)

    # Example usage BRE:
    pairwise_comparisons = [
        ('A', 'B', 1),
        ('A', 'C', -1),
        ('B', 'C', 0),
        ('A', 'D', 1),
        ('B', 'D', -1),
        ('C', 'D', 1)
    ]

    ranks = main.balanced_rank_estimation(pairwise_comparisons)
    print("Balanced Rank Estimation:", ranks)

    # Example usage
    data = [
        (1, 2, 3),
        (2, 3, 4),
        (3, 4, 5),
        (4, 5, 6),
        (5, 6, 7)
    ]

    curve_names = ["Curve 1", "Curve 2"]

    plotStuff.plot_curves(data, curve_names)

# below to trash
""" 
    # so far all tests
    # TODO only n*log(n)
    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i-1][2], S[j][2], S[i-1][1], S[j][1])
            #print("for " + S[i-1][0] + " and " + S[j][0] + " Claire test says: " + str(b))
            if b:
                #print("Welch test can be used")
                t_stat, p_value, conclusion = welch_ttest(S[i-1][3], S[j][3])
                #print(t_stat, p_value, conclusion)
                tabStat.append(t_stat)
                tabPValues.append(float(p_value))
                comp=0 # not significant
                if p_value < 0.05 and t_stat < 0:
                    comp = -1
                if p_value < 0.05 and t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[i-1][0],S[j][0],comp))
            else:
                #print("Permutation test is used")
                observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[i-1][3], S[j][3])
                #print(f"Observed Welch's t-statistic: {observed_t_stat}")
                #print(f"P-value: {p_value}")
                #print(f"conclusion: {conclusion}")
                #print(observed_t_stat, p_value, conclusion)
                tabStat.append(observed_t_stat)
                tabPValues.append(float(p_value))
                comp = 0  # not significant
                if p_value < 0.05 and observed_t_stat < 0:
                    comp = -1
                if p_value < 0.05 and observed_t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[i - 1][0], S[j][0], comp))
"""

""" 
def computeStats(queryValues, vals, conn):
    S = []

    #  statistics (skew and size) for all members
    resultValues = execute_query(conn, queryValues)
    data = []
    for row in resultValues:
        data.append(float(row[0]))
    nvalues = len(data)
    data = np.array(data)
    skewness = compute_skewness(data)
    S.append(('ALL', nvalues, skewness))

    #print("Size of the data: " + str(nvalues))
    #print(f"Skewness of the data: {skewness}")

    # compute stats for every member of the hypothesis
    for v in vals:
        queryVal = queryValues.replace(str(vals), "('" + str(v) + "')")
        # print(queryVal)
        resultValues = execute_query(conn, queryVal)
        data = []
        for row in resultValues:
            data.append(float(row[0]))
        nvalues = len(data)
        data = np.array(data)
        # Compute skewness
        skewness = compute_skewness(data)
        # Print result
        # print("Size of the data: " + str(nvalues))
        # print(f"Skewness of the data: {skewness}")
        S.append((v, nvalues, skewness))

    return S
"""