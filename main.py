import random
from itertools import chain, combinations
import math
import numpy as np
from dbStuff import execute_query, connect_to_db, close_connection
from statStuff import benjamini_hochberg_gpt, welch_ttest, permutation_test, compute_skewness, compute_kendall_tau, benjamini_hochberg


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






def emptyGB(conn):
    # compute top sizeOfVals members - could look for best combinations?
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + " limit " + str(sizeOfVals) + ";")

    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " group by " + sel +  ";")

    resultEmptyGb = execute_query(conn, queryEmptyGb)

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

def getSample(conn, meas, measBase, table, sel, sampleSize, method):
    # sampling using postgresql: https://www.postgresql.org/docs/current/sql-select.html#SQL-FROM
    # system (faster) is block based, bernouili (slower) is row based
    # TODO try bernouilli method instead
    querySample = (
            "SELECT " + sel + ", " + measBase + " FROM " + table + " TABLESAMPLE " + method + " (" + str(sampleSize) + ");")
    # print(querySample)
    resultVals = execute_query(conn, querySample)
    print("Sample size in tuples: ", len(resultVals))
    return resultVals

def generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, method):

    resultVals = getSample(conn, meas, measBase, table, sel, sampleSize, method)

    Sels = tuple([x[0] for x in resultVals])
    Sels = list(dict.fromkeys(Sels))
    #print(Sels)
    Vals = tuple([x[1] for x in resultVals])
    #print(Vals)

    S = []
    for v in Sels:
        #querySample = ("SELECT " + measBase + " FROM " + table + " TABLESAMPLE SYSTEM (" + sampleSize + ") WHERE " + sel + "='" + v + "';")

        #resultSample = execute_query(conn, querySample)
        #print(resultSample)

        data = []
        for row in resultVals:
            if row[0]==v:
                data.append(float(row[1]))

        nvalues = len(data)
        data = np.array(data)
        skewness = compute_skewness(data)
        S.append((v, nvalues, skewness, data))

    #print(S)

    #according to  Wauthier &al JMLR 2013, nlog(n) comparisons enough for recovering the true ranking
    nbOfComparisons=len(Sels)*math.log(len(Sels),2)
    print("Number of comparisons to make: " + str(nbOfComparisons))

    tabPValues=[]
    pairwiseComparison=[]
    tabStat=[]

    # compute Claire statistics for all pairs
    claireTab=[]
    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i-1][2], S[j][2], S[i-1][1], S[j][1])
            claireTab.append((S[i-1][0], S[j][0],b, S[i - 1][3], S[j][3]))

    #print("Claire Tab: ", claireTab)

    # compute all Welch tests
    nbWelch=0
    for cl in claireTab:
        if cl[2]:
            # print("Welch test can be used")
            nbWelch=nbWelch+1
            t_stat, p_value, conclusion = welch_ttest(cl[3], cl[4])
            # print(t_stat, p_value, conclusion)
            tabStat.append(t_stat)
            tabPValues.append(float(p_value))
            comp = 0  # not significant
            if p_value < 0.05 and t_stat < 0:
                comp = -1
            if p_value < 0.05 and t_stat > 0:
                comp = 1
            pairwiseComparison.append((cl[0], cl[1], comp))
    nbremaining = nbOfComparisons - nbWelch
    nbPermut=0
    for cl in claireTab:
        if not cl[2]:
            # print("Permutation test is used")
            nbPermut=nbPermut+1
            observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(cl[3], cl[4])
            # print(f"Observed Welch's t-statistic: {observed_t_stat}")
            # print(f"P-value: {p_value}")
            # print(f"conclusion: {conclusion}")
            # print(observed_t_stat, p_value, conclusion)
            tabStat.append(observed_t_stat)
            tabPValues.append(float(p_value))
            comp = 0  # not significant
            if p_value < 0.05 and observed_t_stat < 0:
                comp = -1
            if p_value < 0.05 and observed_t_stat > 0:
                comp = 1
            pairwiseComparison.append((cl[0], cl[1], comp))
        if nbremaining - nbPermut <= 0:
            break

    print("raw pairwise comparisons: ", pairwiseComparison)
    print("size: ", len(pairwiseComparison))

    # the ones already compared
    alreadyCompared = set()
    for item1, item2, comp in pairwiseComparison:
        alreadyCompared.add(item1)
        alreadyCompared.add(item2)

    allValues = set()
    for s in Sels:
        allValues.add(s)

    difference = allValues.difference(alreadyCompared)
    if len(difference) != 0:
        # compare the ones in difference with one already compared
        for d in difference:
            for cl in claireTab:
                if not cl[0] == d or cl[1] == d:
                    # print("Permutation test is used")
                    nbPermut = nbPermut + 1
                    observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(cl[3], cl[4])
                    # print(f"Observed Welch's t-statistic: {observed_t_stat}")
                    # print(f"P-value: {p_value}")
                    # print(f"conclusion: {conclusion}")
                    # print(observed_t_stat, p_value, conclusion)
                    tabStat.append(observed_t_stat)
                    tabPValues.append(float(p_value))
                    comp = 0  # not significant
                    if p_value < 0.05 and observed_t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and observed_t_stat > 0:
                        comp = 1
                    pairwiseComparison.append((cl[0], cl[1], comp))
                break

    # Benjamini Hochberg correction
    # to be checked

    alpha = 0.05
    # rejected, corrected_p_values = benjamini_hochberg_gpt(tabPValues, alpha)
    # print("Rejected hypotheses:", rejected)
    #print("raw p-values:", tabPValues)
    # print("Corrected p-values (gpt):", corrected_p_values)
    corrected=benjamini_hochberg(tabPValues, alpha)
    #print("Corrected p-values (scipy):", corrected)

    i=0
    nbChanges=0
    for c in pairwiseComparison:
        comp=0
        if corrected[i] < 0.05 and tabStat[i] < 0:
            comp=-1
        if corrected[i] > 0.05 and tabStat[i] > 0:
            comp=1
        if c[2] != comp:
            nbChanges=nbChanges+1
            pairwiseComparison.remove(c)
            pairwiseComparison.append((c[0], c[1], comp))

    print("Number of BH corrections: ", nbChanges, " ratio: ", nbChanges/len(tabPValues), "%")

    #print(pairwiseComparison)



    # ranking
    ranks = balanced_rank_estimation(pairwiseComparison)
    #print("Balanced Rank Estimation:", ranks)

    sorted_items = sorted(ranks.items(), key=lambda item: item[1], reverse=True)
    #print(sorted_items)
    hypothesis=[]
    rank=0
    for s in sorted_items:
        if rank == 0:
            rank=1
            hypothesis.append((s[0],rank))
            val=s[1]
        else:
            if s[1] == val:
                hypothesis.append((s[0], rank))
            else:
                rank=rank+1
                hypothesis.append((s[0], rank))
    #print(hypothesis)


    return hypothesis



def hoeffdingForRank(groupbyAtt, n, hypothesis):

    print("Size of confidence interval around p: " + str(epsilon))
    print("Probability is of making a mistake: " + str(alpha))

    # n queries enough according to Hoeffding
    print("n: " + str(n))

    # compute powerset of categorical attributes
    pwset = powerset(groupbyAtt)

    # empty group by set removed from powerset
    # since it WAS used to generate the hypothesis
    # TODO hypothesis could also be user given
    pwset.remove(())

    #print("Hypothesis is:" + str(hypothesis))


    nbTests = 0
    for i in range(n):

        # generate the random query
        # query, queryHyp, queryValues, queryExcept, strgb = generateRandomQuery(pwsert,hypothesis)
        queryValues, queryCountGb, queryCountExcept = generateRandomQuery(pwset, hypothesis)


        # strategy: use the db engine to check whether the hypothesis holds
        # cross join hypothesis to chosen group by set
        # compute the actual ranks of vals in the hypothesis for each group
        # except all the actual ranks with the hypothesis
        # remaining only the groups where hypothesis does not hold

        resultCountGb = execute_query(conn, queryCountGb)
        resultCountExcept = execute_query(conn, queryCountExcept)

        print("number of tuples checked: " + str(resultCountGb[0][0]))
        print("number of exceptions: " + str(resultCountExcept[0][0]))
        print("ratio is: " + str(resultCountExcept[0][0] / resultCountGb[0][0]))

        ratio = resultCountExcept[0][0] / resultCountGb[0][0]

        # keep actual ratio
        # could use Bernouilli random variables instead: 1 if less than 10% errors
        if ratio > threshold:
            H.append(ratio)
        #    H.append(0)
        else:
            H.append(ratio)
        #    H.append(1)

        nbTests = nbTests + 1

        print("Expected value is: " + str(sum(H) / len(H)))



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
threshold=0.1 # 10% of tuples violating the order
n=math.log(2/alpha,10) / pow(2,epsilon*epsilon)
print("n>= " + str(n))
n=math.ceil(n)

sampleSize=20
samplingMethod='BERNOULLI' # or SYSTEM

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

        #generate hypothesis: ordering of members such that mean is greater (statistically significant on sample)
        hypothesis=generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, samplingMethod)

        print("Hypothesis as predicted: ", hypothesis)


        # empty group by set removed from powerset
        # since it is used to generate the hypothesis
        # TODO hypothesis could also be user given
        #pwset.remove(())

        emptyGB=emptyGB(conn);
        print("Empty GB says:", emptyGB)

        # compute kendall tau between hypothesis and emptyGB
        hypothesis.sort(key=lambda x: x[0])
        emptyGB.sort(key=lambda x: x[0])

        #print(hypothesis)
        #print(emptyGB)

        rankings_with_ties1 = [x[1] for x in hypothesis]
        rankings_with_ties2 = [x[1] for x in emptyGB]

        #print(rankings_with_ties1)
        #print(rankings_with_ties2)

        tau, p_value = compute_kendall_tau(rankings_with_ties1, rankings_with_ties2)
        print(f"Kendall Tau-c: {tau}, p-value: {p_value}")

        vals=tuple([x[0] for x in hypothesis])

        hoeffdingForRank(groupbyAtt, n, hypothesis)

        # Close the connection
        close_connection(conn)
