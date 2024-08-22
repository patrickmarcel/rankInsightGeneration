import random
import math
import numpy as np

import rankingFromPairwise
import utilities
from utilities import powerset, jaccard_similarity
from plotStuff import plot_curves
from dbStuff import execute_query, connect_to_db, close_connection
from statStuff import welch_ttest, permutation_test, compute_skewness, compute_kendall_tau, benjamini_hochberg, \
    benjamini_hochberg_statmod
import time
from rankingFromPairwise import computeRanksForAll, merge_sort

from scipy.stats import permutation_test as sp_permutation_test

import configparser
import json

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

    stat = abs((skew1 / count1) - (skew2 / count2))
    if stat < threshold:
        # test can be used
        return True
    else:
        return False


def emptyGB(conn, nb):
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + " limit " + str(sizeOfVals) + ";")
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + ";")
    #queryEmptyGb = ("SELECT " + sel + ","
    #               + " rank () over (  order by " + meas + " desc ) as rank" +
    #              " FROM " + table +  " WHERE " + sel + " in (" + hyp + ") group by " + sel +  ";")

    #hyp = ""
    #for i in range(len(vals)):
    #    hyp = hyp + "'" + str(vals[i]) + "'"
    #    if i != len(vals) - 1:
    #        hyp = hyp + ","

    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " group by " + sel + " limit " + str(nb) + ";")

    #print(queryEmptyGb)
    resultEmptyGb = execute_query(conn, queryEmptyGb)

    queryEmptyGbAll = ("SELECT " + sel + ","
                       + " rank () over (  order by " + meas + " desc ) as rank" +
                       " FROM " + table + " group by " + sel + ";")
    resultEmptyGbAll = execute_query(conn, queryEmptyGbAll)

    return resultEmptyGb, resultEmptyGbAll


def generateRandomQuery(pwsert, valsToSelect, hypothesis):
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

    #print("vals in gen queries:", valsToSelect)
    hyp = ""
    for i in range(len(hypothesis)):
        hyp = hyp + str(hypothesis[i])
        if i != len(hypothesis) - 1:
            hyp = hyp + ","
    queryHyp = (
            "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")

    query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
             + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
             " FROM " + table + " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ")

    queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "
                   + table + " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ) x;")

    queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + " ")

    queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
    queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")

    #return query, queryHyp, queryValues, queryExcept, strgb
    return queryValues, queryCountGb, queryCountExcept


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

    querySample = (
            "SELECT " + sel + ", " + measBase + " FROM " + table + " TABLESAMPLE " + method + " (" + str(
        sampleSize) + ");")
    # print(querySample)
    resultVals = execute_query(conn, querySample)
    print("Sample size in tuples: ", len(resultVals))
    return resultVals


def computeBHcorrection(pairwiseComparison, alpha=0.05):
    # rejected, corrected_p_values = benjamini_hochberg_gpt(tabPValues, alpha)
    # print("Rejected hypotheses:", rejected)
    # print("raw p-values:", tabPValues)
    # print("Corrected p-values (gpt):", corrected_p_values)
    tabPValues = []
    for p in pairwiseComparison:
        tabPValues.append(p[4])

    corrected = benjamini_hochberg(tabPValues, alpha)
    rejected, corrected2 = benjamini_hochberg_statmod(tabPValues, alpha)

    print("nb of True in rejected: ", utilities.nbTrueInList(rejected))

    pairwiseComp2 = []
    i = 0
    nbChanges = 0
    for c in pairwiseComparison:
        comp = 0  # not significant
        if corrected[i] < 0.05 and c[3] < 0:
            comp = -1
        if corrected[i] < 0.05 and c[3] > 0:
            comp = 1
        if comp != c[2]:
            nbChanges = nbChanges + 1
        pairwiseComp2.append((c[0], c[1], comp, c[2], corrected[i]))
        i = i + 1

    print("Number of BH corrections: ", nbChanges)

    print("nb non zeros after corrections: ", utilities.countNonZeros(pairwiseComp2))

    return pairwiseComp2


def generateComparisonsWithMergeSort(Sels, S):
    # compute Claire statistics for all pairs
    claireTab = []
    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i - 1][2], S[j][2], S[i - 1][1], S[j][1])
            claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))

    #claireComp = [(x[0],x[1],x[2]) for x in claireTab]
    #print("claireComp: ", claireTab)
    print("merge: ", merge_sort(Sels, claireTab))
    print("pairwise: ", rankingFromPairwise.pairwiseComparison)

    return rankingFromPairwise.pairwiseComparison


def generateAllComparisons(Sels, S, nbOfComparisons):
    #tabPValues = []
    pairwiseComparison = []
    #tabStat = []

    # compute Claire statistics for all pairs
    claireTab = []
    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i - 1][2], S[j][2], S[i - 1][1], S[j][1])
            claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))

            if b:
                # print("Welch test can be used")
                t_stat, p_value, conclusion = welch_ttest(S[i - 1][3], S[j][3])
                # print(t_stat, p_value, conclusion)
                #tabStat.append(t_stat)
                #tabPValues.append(float(p_value))
                comp = 0  # not significant
                if p_value < 0.05 and t_stat < 0:
                    comp = -1
                if p_value < 0.05 and t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[i - 1][0], S[j][0], comp, t_stat, float(p_value)))
            else:
                # print("Permutation test is used")
                observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[i - 1][3], S[j][3])
                # print(f"Observed Welch's t-statistic: {observed_t_stat}")
                # print(f"P-value: {p_value}")
                # print(f"conclusion: {conclusion}")
                # print(observed_t_stat, p_value, conclusion)
                #tabStat.append(observed_t_stat)
                #tabPValues.append(float(p_value))
                comp = 0  # not significant
                if p_value < 0.05 and observed_t_stat < 0:
                    comp = -1
                if p_value < 0.05 and observed_t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[i - 1][0], S[j][0], comp, observed_t_stat, float(p_value)))

    pairwiseComparison = computeBHcorrection(pairwiseComparison, 0.05)
    return pairwiseComparison


def generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, method):
    resultVals = getSample(conn, meas, measBase, table, sel, sampleSize, method)

    # get adom values
    Sels = list(set([x[0] for x in resultVals]))

    #analyse sample for each adom value: value, nb of measures, skewness, and tuples
    S = []
    for v in Sels:

        data = []
        for row in resultVals:
            if row[0] == v:
                data.append(float(row[1]))

        nvalues = len(data)
        data = np.array(data)
        skewness = compute_skewness(data)
        S.append((v, nvalues, skewness, data))

    #print(S)

    # nlog(n) comparisons enough for recovering the true ranking when commarisons are certain (not noisy)
    # we should try less
    nbOfComparisons = len(Sels) * math.log(len(Sels), 2)
    #print("Number of comparisons to make: " + str(nbOfComparisons))

    #pairwiseComparison=generateAllComparisons(Sels, S, nbOfComparisons)

    pairwiseComparison = generateComparisonsWithMergeSort(Sels, S)

    # ranking
    #ranks = balanced_rank_estimation(pairwiseComparison)
    #print("Balanced Rank Estimation:", ranks)
    ranks = computeRanksForAll(pairwiseComparison, Sels)

    sorted_items = sorted(ranks.items(), key=lambda item: item[1], reverse=True)

    # Construct a rank from the number of comparison won for each adom values
    hypothesis = []
    rank = 0
    for s in sorted_items:
        if rank == 0:
            rank = 1
            hypothesis.append((s[0], rank))
            val = s[1]
        else:
            if s[1] == val:
                hypothesis.append((s[0], rank))
                val = s[1]
            else:
                rank = rank + 1
                hypothesis.append((s[0], rank))
                val = s[1]

    return hypothesis


def hoeffdingForRank(groupbyAtt, n, valsToSelect, limitedHyp):
    print("Size of confidence interval around p: " + str(epsilon))
    print("Probability is of making a mistake: " + str(alpha))

    # n queries enough according to Hoeffding
    print("n: " + str(n))

    # compute powerset of categorical attributes
    pwset = powerset(groupbyAtt)

    # empty group by set removed from powerset
    # since it WAS used to generate the hypothesis

    pwset.remove(())

    #print("Hypothesis is:" + str(hypothesis))

    nbTests = 0
    for i in range(n):

        # generate the random query
        # query, queryHyp, queryValues, queryExcept, strgb = generateRandomQuery(pwsert,hypothesis)
        queryValues, queryCountGb, queryCountExcept = generateRandomQuery(pwset, valsToSelect, limitedHyp)

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

        #nbTests = nbTests + 1

    expectedValue = sum(H) / len(H)
    print("Expected value is: " + str(sum(H) / len(H)))
    return expectedValue


#parameters

# number of values of adom to consider - top ones after hypothesis is generated
nbAdomVals = 5

# for Hoeffding
epsilon = 0.01
alpha = 0.01
p = 0
H = []
threshold = 0.1  # 10% of tuples violating the order
n = math.log(2 / alpha, 10) / pow(2, epsilon * epsilon)
n = math.ceil(n)

# for DB sampling
sampleSize = 30
samplingMethod = 'SYSTEM'  # or SYSTEM

nbruns = 2

# TODO
#  check stability of hypothesis - BRE was generated by GPT...
#   also ranking is either all nb 1 or total order? also tau is nan when all 1's
#  check expected values, too good to be true
#  loop over sample size from 5 to 95 or so
#  change sel and measures
#  use other databases
#  hypothesis could also be user given

if __name__ == "__main__":

    config = configparser.ConfigParser()
    config.read('configs/flights.ini')
    USER = "AC"
    # Database connection parameters
    dbname = config[USER]['dbname']
    user = config[USER]['user']
    password = config[USER]['password']
    host = config[USER]['host']
    port = int(config[USER]['port'])
    # Cube info
    table = config["Common"]['table']
    measures = json.loads(config.get("Common", "measures"))
    groupbyAtt = json.loads(config.get("Common", "groupbyAtt"))
    sel = config["Common"]['sel']
    meas = config["Common"]['meas']
    measBase = config["Common"]['measBase']

    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    if conn:

        tabHypo = []
        resultRuns = []
        for i in range(nbruns):
            rankingFromPairwise.pairwiseComparison = []

            start_time = time.time()
            #generate hypothesis: ordering of members such that mean is greater (statistically significant on sample)
            hypothesis = generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, samplingMethod)

            #print("Hypothesis as predicted: ", hypothesis)

            # limit hypothesis to top nbAdomVals
            limitedHyp = []
            valsToSelect = []
            j = 0
            for h in hypothesis:
                if (h[1] <= nbAdomVals and j < nbAdomVals):
                    limitedHyp.append(h)
                    valsToSelect.append(h[0])
                    j = j + 1
            #print("Hypothesis limited: ", limitedHyp)
            #print("vals: ",valsToSelect)

            # should only be done once
            emptyGBresult, emptyGBresultAll = emptyGB(conn, nbAdomVals)
            print("Empty GB says:", emptyGBresult)

            # compute kendall tau between hypothesis and emptyGB
            limitedHyp.sort(key=lambda x: x[0])
            emptyGBresult.sort(key=lambda x: x[0])
            hypothesis.sort(key=lambda x: x[0])
            emptyGBresultAll.sort(key=lambda x: x[0])

            # record all hypotheses
            tabHypo.append(limitedHyp)

            # should also compute tau between hypothesis of different runs

            #print(hypothesis)
            #print(emptyGB)

            rankings_with_ties1 = [x[1] for x in hypothesis]
            rankings_with_ties2 = [x[1] for x in emptyGBresultAll]

            #print(rankings_with_ties1)
            #print(rankings_with_ties2)

            tau, p_value = compute_kendall_tau(rankings_with_ties1, rankings_with_ties2)
            print("Tau-c between hypothesis and emptyGBAll: ", tau, "p-value: ", p_value)

            #todo should also compute for limitedHyp and emptyGB

            #vals=tuple([x[0] for x in limitedHyp])
            #print("********** vals nnd vals2sel ")
            #print(vals)
            #print(tuple(valsToSelect))

            expected = hoeffdingForRank(groupbyAtt, n, tuple(valsToSelect), limitedHyp)

            end_time = time.time()
            elapsed_time = end_time - start_time

            resultRuns.append((i, float(tau), expected, elapsed_time))

        print(resultRuns)
        print(tabHypo)
        # compute hamming dist in tabhypo or jaccard

        names = ['tau', 'expected', 'time']
        title = 'Sample size=' + str(sampleSize)
        plot_curves(resultRuns, names, 'time', 'expected', title)
        # Close the connection
        close_connection(conn)
