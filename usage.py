import numpy as np
import main
import statStuff

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
