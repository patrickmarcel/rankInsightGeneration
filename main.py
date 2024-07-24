import psycopg2
from psycopg2 import sql
import random
from itertools import chain, combinations
import math
import scipy.stats as stats
from scipy import stats
from scipy.stats import skew
import numpy as np
from scipy.stats import f

import pandas as pd


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



def benjamini_hochberg(p_values, alpha=0.05):
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

def connect_to_db(dbname, user, password, host='localhost', port='5432'):
    """
    Establishes a connection to the PostgreSQL database.

    :param dbname: Name of the database
    :param user: Database user
    :param password: User's password
    :param host: Database host address (default is 'localhost')
    :param port: Connection port number (default is '5432')
    :return: Connection object
    """
    try:
        conn = psycopg2.connect(
            dbname=dbname,
            user=user,
            password=password,
            host=host,
            port=port
        )
        print("Connection to database established successfully.")
        return conn
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return None


def execute_query(conn, query):
    """
    Executes a given SQL query using the established connection.

    :param conn: Connection object
    :param query: SQL query to be executed
    :return: Query result
    """
    try:
        cursor = conn.cursor()
        cursor.execute(query)
        conn.commit()
        #print("Query executed successfully.")

        try:
            result = cursor.fetchall()
            return result
        except psycopg2.ProgrammingError:
            # If the query does not return any data (like an INSERT or UPDATE)
            return None
    except Exception as e:
        print(f"Error executing query: {e}")
        return None
    finally:
        cursor.close()

def printResultSet(result):
    if result is not None:
        for row in result:
            print(row)

def close_connection(conn):
    """
    Closes the database connection.

    :param conn: Connection object
    """
    try:
        conn.close()
        print("Database connection closed.")
    except Exception as e:
        print(f"Error closing connection: {e}")

def powerset(s):
    """
    Generates the powerset of a given set s.

    :param s: The input set
    :return: A list of subsets representing the powerset
    """
    s = list(s)
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s) + 1)))


def generateGB(groupAtts):
    """
    Generates group bys and sel from the list of all categorical attributes - unused so far

    :param groupAtts: all the categorical attributes
    :return:
    """
    for g in groupAtts:
        gb = groupAtts.remove(g)
        groupbyAtt = gb
        sel = g
        for m in measures:
            meas = "sum(" + m + ")"

            queryVals = ("select distinct " + sel + " from " + table + ";")


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

def claireStat(skew1, skew2, count1, count2, threshold=0.049):
    """
    Compute Claire statistics for testing if Welch test can be used

    Parameters:
    - skew1: skewness of first samples
    - skew2: skewness of second sample
    - count1: size of first sample
    - count2: size of second sample
    - threshold: why would we want to change?

    Returns:
    - skewness: boolean indicating if Welch test can be used
    """

    stat = abs((skew1/count1) - (skew2/count2))
    if stat<threshold:
        # test can be used
        return True
    else:
        return False





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





def ExampleUsage():

# Example usage for skewness

    # Sample data
    data = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5])
    # Compute skewness
    skewness = compute_skewness(data)
    # Print result
    print(f"Skewness of the data: {skewness}")

# Example usage for Welch

    # Sample data
    group1 = np.array([23, 45, 56, 67, 78, 89, 90, 56, 54, 67])
    group2 = np.array([35, 56, 76, 78, 45, 34, 23, 45, 67, 88])

    # Perform Welch's t-test
    t_stat, p_value, conclusion = welch_ttest(group1, group2)

    # Print results
    print(f"t-statistic: {t_stat}")
    print(f"p-value: {p_value}")
    print(conclusion)

# Example usage for Brown Forsythe
    group1 = [23, 20, 25, 22, 21, 20, 23, 24]
    group2 = [22, 21, 24, 21, 20, 19, 21, 23]
    group3 = [26, 24, 27, 25, 26, 27, 28, 29]

    W, p_value = brown_forsythe_test(group1, group2, group3)
    print(f"Brown-Forsythe test statistic: {W}")
    print(f"P-value: {p_value}")

# Example usage for permutation test:
    np.random.seed(0)
    x1 = np.random.normal(0, 1, 30)
    x2 = np.random.normal(0.5, 1.5, 35)

    observed_t_stat, p_value, permuted_t_stats = permutation_test(x1, x2)
    print(f"Observed Welch's t-statistic: {observed_t_stat}")
    print(f"P-value: {p_value}")

    # Example usage of benjamini_hochberg:
    p_values = [0.01, 0.04, 0.03, 0.002, 0.05, 0.001, 0.03, 0.04]
    alpha = 0.05
    rejected, corrected_p_values = benjamini_hochberg(p_values, alpha)

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

    ranks = balanced_rank_estimation(pairwise_comparisons)
    print("Balanced Rank Estimation:", ranks)


def generateHypothesis(conn):
    # compute top sizeOfVals members - could look for best combinations?
    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " group by " + sel + " limit " + str(sizeOfVals) + ";")

    resultEmptyGb = execute_query(conn, queryEmptyGb)

    """
    if resultEmptyGb is not None:
        print("Overall the ranking is as follows:")
        for row in resultEmptyGb:
            print(row)
    """

    return resultEmptyGb



def generateRandomQuery(pwsert,hypothesis):
    nb = random.randint(0, len(pwsert) - 1)
    gb = pwsert[nb]
    strgb = ""
    for i in range(len(gb)):
        strgb = strgb + str(gb[i])
        if i != len(gb) - 1:
            strgb = strgb + ","

    print("group by is: " + strgb)

    # for debugging
    # strgb = "departure_airport"

    hyp = ""
    for i in range(len(hypothesis)):
        hyp = hyp + str(hypothesis[i])
        if i != len(hypothesis) - 1:
            hyp = hyp + ","
    queryHyp = (
            "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")

    query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
             + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
             " FROM " + table + " WHERE " + sel + " in " + str(vals) + " group by " + strgb + "," + sel + " ")

    queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "
                   + table + " WHERE " + sel + " in " + str(vals) + " group by " + strgb + "," + sel + " ) x;")


    queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + " ")

    queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
    queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")

    #return query, queryHyp, queryValues, queryExcept, strgb
    return queryValues,queryCountGb,queryCountExcept


def getValues(queryValues, vals, v, conn):
    queryVal = queryValues.replace(str(vals), "('" + v + "')")
    resultValues = execute_query(conn, queryVal)
    data = []
    for row in resultValues:
        data.append(float(row[0]))

    np.array(data)
    return data


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


def generateHypothesisTest(conn, meas, measBase, table, sel):

    queryVals = ("SELECT  " + sel + ", " + meas + " FROM " + table + " group by " + sel + " order by " + meas + "  desc limit " + str(sizeOfVals) + ";")
    #print(queryVals)
    #queryVals = ("SELECT DISTINCT " + sel +  " FROM " + table + " limit " + str(sizeOfVals) + ";")
    resultVals = execute_query(conn, queryVals)
    Vals = tuple([x[0] for x in resultVals])

    S = []
    for v in Vals:
        querySample = ("SELECT " + measBase + " FROM " + table + " TABLESAMPLE SYSTEM (10) WHERE " + sel + "='" + v + "';")

        resultSample = execute_query(conn, querySample)
        #print(resultSample)

        data = []
        for row in resultSample:
            data.append(float(row[0]))

        nvalues = len(data)
        data = np.array(data)
        skewness = compute_skewness(data)
        S.append((v, nvalues, skewness, data))

    #print(S)

    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i-1][2], S[j][2], S[i-1][1], S[j][1])
            print("for " + S[i-1][0] + " and " + S[j][0] + " Claire test says: " + str(b))
            if b:
                print("Welch test can be used")
                #sample1 = getValues(queryValues, vals, S[1][0], conn)
                #sample2 = getValues(queryValues, vals, S[i][0], conn)
                t_stat, p_value, conclusion = welch_ttest(S[i-1][3], S[j][3])
                print(t_stat, p_value, conclusion)
            else:
                print("Permutation test is used")
                observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[i-1][3], S[j][3])
                print(f"Observed Welch's t-statistic: {observed_t_stat}")
                print(f"P-value: {p_value}")
                print(f"conclusion: {conclusion}")


#def computeStatsDB():
    # for all categorical attributes
    # for all values
    # send query to collect average, size, etc.

#queryPattern="SELECT " + gb + "," + meas + " FROM " + table + " WHERE " + sel + " in " + vals + "group by " + gb +";"
#hypothesis=[('AA', 1),('UA', 2),('US', 3)]


table="fact_table"
measures=["nb_flights","departure_delay","late_aircraft"]
groupAtts=["departure_airport","date","departure_hour","flight","airline"]

measure="nb_flights"
selectionAtt="airline"
groupbyAtt=["departure_airport","date","departure_hour","flight"]
gb=""
sel="airline"
vals="('AA','UA','US')"
#meas="sum(nb_flights)"
meas="avg(departure_delay)"
measBase="departure_delay"

sizeOfVals=7 #number of members for the hypothesis
q0=""
epsilon=0.01
alpha=0.01
p=0
H=[]
nbTests=0
threshold=0.1 # 10% of tuples violating the order
n=math.log(2/alpha,10) / pow(2,epsilon*epsilon)
print("n>= " + str(n))
n=math.ceil(n)

if __name__ == "__main__":

    # Database connection parameters
    dbname = "flight_dw"
    user = ""
    password = ""
    host = "localhost"
    port = "5432"

    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    if conn:

        #test genertion hypothesis: only if statistically significant on samples
        generateHypothesisTest(conn, meas, measBase, table, sel)


        #compute powerset of categorical attributes
        pwset = powerset(groupbyAtt)

        # generating the hypothesis
        # hypothesis is a ranking of members

        # empty group by set removed from powerset
        # since it is used to generate the hypothesis
        # TODO hypothesis could also be user given
        pwset.remove(())

        hypothesis=generateHypothesis(conn);
        vals=tuple([x[0] for x in hypothesis])

        # n queries enough according to Hoeffding
        print("n: " + str(n))
        for i in range(n):

            # generate the random query
            #query, queryHyp, queryValues, queryExcept, strgb = generateRandomQuery(pwsert,hypothesis)
            queryValues,queryCountGb,queryCountExcept=generateRandomQuery(pwset, hypothesis)

            # compute skew and size for the members of vals in the result of the random query
            S=computeStats(queryValues, vals, conn)

            for i in range(2,len(S)):
                # example of statstical test for first 2 members
                b=claireStat(S[1][2], S[i][2], S[1][1], S[i][1])
                print("for " + S[1][0] + " and " + S[i][0] + " Claire test says: " + str(b))
                if b:
                    print("Welch test can be used")
                    sample1 = getValues(queryValues,vals, S[1][0], conn)
                    sample2 = getValues(queryValues,vals, S[i][0], conn)
                    t_stat, p_value, conclusion=welch_ttest(sample1, sample2)
                    print(conclusion)


            #strategy: use the db engine to check whether the hypothesis holds
            #cross join hypothesis to chosen group by set
            # compute the actual ranks of vals in the hypothesis for each group
            # except all the actual ranks with the hypothesis
            # remaining only the groups where hypothesis does not hold


            resultCountGb = execute_query(conn, queryCountGb)
            resultCountExcept = execute_query(conn, queryCountExcept)

            print("number of tuples checked: " + str( resultCountGb[0][0] ))
            print("number of exceptions: " + str(resultCountExcept[0][0]))
            print("ratio is: " + str(resultCountExcept[0][0]  / resultCountGb[0][0]))

            ratio=resultCountExcept[0][0]  / resultCountGb[0][0]


            # keep actual ratio
            # could use Bernouilli random variables instead: 1 if less than 10% errors
            if ratio > threshold:
                H.append(ratio)
            #    H.append(0)
            else:
                H.append(ratio)
            #    H.append(1)

            nbTests=nbTests+1

            print("Expected value is: " + str(sum(H)/len(H)))
            print("Size of confidence interval around p: " + str(epsilon))
            print("Probability is of making a mistake: " + str(alpha))

            # Close the connection
        close_connection(conn)


