import numpy as np
import scipy.stats as stats
from scipy import stats
from scipy.stats import skew
from scipy.stats import f
from scipy.stats import false_discovery_control

def benjamini_hochberg(p_values, alpha=0.05):
    return false_discovery_control(p_values)

def benjamini_hochberg_gpt(p_values, alpha=0.05):
    """
    Perform Benjamini-Hochberg correction for multiple hypothesis testing.

    Parameters:
    - p_values: list or array of p-values from individual hypothesis tests.
    - alpha: desired FDR control level (default is 0.05).

    Returns:
    - rejected: boolean array indicating which hypotheses are rejected.
    - corrected_p_values: array of adjusted p-values.
    """
    p_values = np.asarray(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]

    # Compute the Benjamini-Hochberg critical values
    bh_critical_values = np.arange(1, n + 1) * alpha / n

    # Determine the largest k where p(k) <= (k/n) * alpha
    rejected = sorted_p_values <= bh_critical_values
    max_k = np.max(np.where(rejected)) if np.any(rejected) else -1

    # Adjust p-values
    corrected_p_values = np.minimum.accumulate((n / np.arange(n, 0, -1)) * sorted_p_values)
    corrected_p_values = np.minimum(corrected_p_values, 1.0)

    # Reorder corrected p-values to match original order
    corrected_p_values = corrected_p_values[np.argsort(sorted_indices)]

    # Determine which hypotheses are rejected
    rejected = np.zeros(n, dtype=bool)
    if max_k >= 0:
        rejected[sorted_indices[:max_k + 1]] = True

    return rejected, corrected_p_values



def welch_t_statistic(x1, x2):
    """Calculate Welch's t-statistic for two samples."""
    mean1, mean2 = np.mean(x1), np.mean(x2)
    var1, var2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    n1, n2 = len(x1), len(x2)
    t_stat = (mean1 - mean2) / np.sqrt(var1 / n1 + var2 / n2)
    return t_stat


def permutation_test(x1, x2, num_permutations=10000, alpha=0.05):
    """Perform a permutation test using Welch's t-statistic as the test statistic."""
    # Calculate observed Welch's t-statistic
    observed_t_stat = welch_t_statistic(x1, x2)

    # Combine the samples
    combined = np.concatenate([x1, x2])

    # Generate the null distribution by permutation
    permuted_t_stats = []
    for _ in range(num_permutations):
        np.random.shuffle(combined)
        perm_x1 = combined[:len(x1)]
        perm_x2 = combined[len(x1):]
        permuted_t_stats.append(welch_t_statistic(perm_x1, perm_x2))

    permuted_t_stats = np.array(permuted_t_stats)

    # Calculate the p-value
    p_value = np.mean(np.abs(permuted_t_stats) >= np.abs(observed_t_stat))

    if p_value < alpha:
        conclusion = "Reject the null hypothesis: There is a significant difference between the means of the two groups."
    else:
        conclusion = "Fail to reject the null hypothesis: There is no significant difference between the means of the two groups."

    return observed_t_stat, p_value, permuted_t_stats, conclusion



def welch_ttest(sample1, sample2, alpha=0.05):
    """
    Perform Welch's t-test to determine if there is a significant difference
    between the means of two groups.

    Parameters:
    - sample1: list or numpy array of sample data for group 1
    - sample2: list or numpy array of sample data for group 2
    - alpha: significance level (default is 0.05)

    Returns:
    - t_stat: the calculated t-statistic
    - p_value: the two-tailed p-value
    - conclusion: string stating if we reject or fail to reject the null hypothesis
    """

    # Perform Welch's t-test
    t_stat, p_value = stats.ttest_ind(sample1, sample2, equal_var=False)

    # Conclusion
    if p_value < alpha:
        conclusion = "Reject the null hypothesis: There is a significant difference between the means of the two groups."
    else:
        conclusion = "Fail to reject the null hypothesis: There is no significant difference between the means of the two groups."

    return t_stat, p_value, conclusion



def compute_skewness(data):
    """
    Compute the skewness of a population.

    Parameters:
    - data: list or numpy array of sample data

    Returns:
    - skewness: the skewness of the data
    """

    # Convert data to a numpy array if it's not already
    #data = np.array(data)

    # Compute skewness
    skewness = skew(data, bias=False)

    return skewness




def brown_forsythe_test(*groups):
    """
    Perform the Brown-Forsythe test for equal variances.

    Parameters:
    *groups : array_like
        Arrays of sample data. There must be at least two arrays.

    Returns:
    W : float
        The Brown-Forsythe test statistic.
    p_value : float
        The p-value for the test.
    """
    # Number of groups
    k = len(groups)

    if k < 2:
        raise ValueError("At least two groups are required")

    # Calculate the group sizes
    n = np.array([len(group) for group in groups])
    N = np.sum(n)

    # Calculate the group medians
    medians = np.array([np.median(group) for group in groups])

    # Calculate the absolute deviations from the medians
    Z = [np.abs(group - median) for group, median in zip(groups, medians)]

    # Calculate the overall mean of the absolute deviations
    Z_flat = np.concatenate(Z)
    Z_mean = np.mean(Z_flat)

    # Calculate the between-group sum of squares
    SSB = np.sum(n * (np.array([np.mean(z) for z in Z]) - Z_mean) ** 2)

    # Calculate the within-group sum of squares
    SSW = np.sum([np.sum((z - np.mean(z)) ** 2) for z in Z])

    # Calculate the Brown-Forsythe test statistic
    dfb = k - 1
    dfw = N - k
    W = (SSB / dfb) / (SSW / dfw)

    # Calculate the p-value
    p_value = f.sf(W, dfb, dfw)

    return W, p_value



def compute_kendall_tau(rankings1, rankings2):
    """
    Compute Kendall Tau-b correlation coefficient for rankings with ties.

    Parameters:
    - rankings1: list, the first ranking.
    - rankings2: list, the second ranking.

    Returns:
    - tau: float, the Kendall Tau-b correlation coefficient.
    - p_value: float, the p-value for testing non-correlation.
    """
    tau, p_value = stats.kendalltau(rankings1, rankings2)
    return tau, p_value