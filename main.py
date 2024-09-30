import random
import math
import numpy as np

import dbStuff
import rankingFromPairwise
import utilities
from utilities import powerset
from plotStuff import plot_curves
from dbStuff import execute_query, connect_to_db, close_connection, getSample, emptyGB
from statStuff import welch_ttest, permutation_test, compute_skewness, compute_kendall_tau, benjamini_hochberg, \
    benjamini_hochberg_statmod, claireStat
import time
from rankingFromPairwise import computeRanksForAll, merge_sort

import configparser
import json

import bernstein


# ------  Debug ?  ------------
DEBUG_FLAG = True


def get_state_sample(conn, measBase, table, sel, sampleSize, state):

    querySample = "SELECT "+sel+", "+measBase+" FROM "+table+" where "+sel+" = '"+state+"' limit "+str(sampleSize)+";"
    #print(querySample)
    resultVals = execute_query(conn, querySample)
    return resultVals

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
    resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=False)
    #resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=DEBUG_FLAG)

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

#This function estimates the number of violations in all the cube of R
#by randomly drawing tuples from the materialized cuboids (R included)
#It uses Hoeffding concentration inequality for bounding the number of draws according to a confidence interval
def estimateViolations(conn, meas, measBase, table, sel, cuboids, ranking, epsilon = 0.1, alpha = 0.1):
    #n is number of draws
    n = math.log(2 / alpha, 10) / (2 * epsilon * epsilon)
    n = math.ceil(n)

    estimates=0

    for i in range(n):
        nCuboid = random.randint(1, len(cuboids)) #+1 is for R itself
        print("nCuboid: ",nCuboid)
        #if nCuboid == len(cuboids):
        #    #draw in R
        #    #TODO draw tuples where only diff is on sel attribute!
        #    tuples=getSample(conn, measBase, table, sel, 5)
        #else:
        #draw in cuboid nCuboid
        view=cuboids[nCuboid][0]
        tuples=getSample(conn, "avg", view, sel, 5)
        if checkViolation(tuples, ranking) == True:
            estimates=estimates+1

    return n, estimates

# returns the rank of value in ranking
# returns 0 if value not found
def getRank(value, ranking):
    for r in ranking:
        if r[0] == value:
            rank=r[1]
    return rank

#checks if tuples violate ranking, return True if this is the case
def checkViolation(tuples, ranking):
    print("tuples: ", tuples)
    print("ranking: ", ranking)

    meast1 = tuples[0][1]
    meast2 = tuples[1][1]

    valuet1 = tuples[0][0]
    valuet2 = tuples[1][0]

    rankt1 = getRank(valuet1,ranking)
    rankt2 = getRank(valuet2, ranking)

    if meast1<meast2 and rankt1<rankt2:
        return True
    else:
        return False

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

# returns true if tuple violates ranking
# TODO check!
def countViolations(conn, viewName, ranking):
    viewDef=dbStuff.getDefOfMV(conn, viewName)
    strgb=viewDef.split("GROUP BY ")[1].split(";")[0]
    queryHyp = (
            "select * from (select  * from " + viewName +  ") t1  cross join (values " + ranking + ") as t2 ")

    query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
             + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
             " FROM " + table +  " group by " + strgb + "," + sel + " ")

    #queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "+ table + " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ) x;")

    queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + " ")

    #queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
    queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")
    return dbStuff.execute_query(conn, queryCountExcept)




#draws n views, return those having less than threshold violations
# TODO check!
def azuma(conn, n, threshold, ranking):
    print(n + " draws, you have "+ (100-n) +"% of chances to get " + math.sqrt(2*n*math.log(n))+ " cuboids with acceptable violations")
    res = dbStuff.getMVnames(conn)
    #print(res)
    for i in range(n):
        nb = random.randint(0, len(res) - 1)
        viewName = res[nb][0]
        print(viewName)
        res.remove(nb)
        view=dbStuff.execute_query(conn, "select * from "+viewName + ";")
        print(view)
        nbViolations=countViolations(view, ranking)
        tabView=[]
        if nbViolations < threshold:
            tabView.append(view)
        return tabView


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
    # The DB wee want
    config.read('configs/flights.ini')
    # The system this is running on
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
    function = config["Common"]['function']


    # number of values of adom to consider - top ones after hypothesis is generated
    nbAdomVals = 5

    # for Hoeffding
    epsilon = 0.1
    alpha = 0.1
    p = 0
    H = []
    threshold = 0.1  # 10% of tuples violating the order
    n = math.log(2 / alpha, 10) / (2* epsilon * epsilon)
    n = math.ceil(n)

    # for DB sampling
    sampleSize = 0.2
    samplingMethod = 'SYSTEM_ROWS'  # or SYSTEM

    if DEBUG_FLAG:
        nbruns = 1
    else:
        nbruns = 10

    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    # to always have the same order in group bys, with sel attribute last
    groupbyAtt.sort()

    # fetch the congressional sample
    adom = [x[0] for x in execute_query(conn, "select distinct  "+sel+" from "+table+";")]
    table_size = execute_query(conn, "select count(1) from " + table + ";")[0][0]

    sample_size = int(table_size * sampleSize)
    alpha = 0.25
    house_size = sample_size * alpha
    senate_size = sample_size * (1 - alpha)

    house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)

    senate = []
    state_sample_size = int(senate_size / len(adom))
    for state in adom:
        senate.extend(get_state_sample(conn, measBase, table, sel, state_sample_size, state))

    congress = house + senate
    # END - fetch the congressional sample

    buckets = {s: [] for s in adom}
    skews = dict()
    for item in congress:
        buckets[item[0]].append(item[1])
        skews[item[0]] = compute_skewness(item[1])

    # do all welch tests
    param_budget = 20
    param_budget = int(param_budget / 2)

    from scipy.stats import ttest_ind

    welch_matrix = [[1 for j in adom] for i in adom]
    w_comparisons = []
    w_comparisons_rej = []

    for i in range(len(adom)):
        for j in range(i + 1):
            left = adom[i]
            right = adom[j]
            res = ttest_ind(buckets[left], buckets[right], equal_var=False)
            welch_matrix[i][j] = res.pvalue
            welch_matrix[j][i] = res.pvalue
            stat_c = claireStat(skews[left], skews[right], len(left), len(right))
            if res.pvalue < 0.05:
                if res.statistic < 0:
                    w_comparisons.append((left, right, stat_c))
                else:
                    w_comparisons.append((right, left, stat_c))
            else:
                if res.statistic < 0:
                    w_comparisons_rej.append((left, right, stat_c))
                else:
                    w_comparisons_rej.append((right, left, stat_c))

    print("NB de comparaisons significatives (welch)", len(w_comparisons))
    #print_comp_list(sorted(w_comparisons, key=lambda x: x[0] + x[1]))
    by_prox_to_threshold = sorted(w_comparisons, key=lambda x: abs(0.05 - x[2]), reverse=True)
    # print(by_prox_to_threshold)

    final = by_prox_to_threshold[param_budget:]
    to_redo = by_prox_to_threshold[:param_budget]

    to_redo.extend(sorted(w_comparisons_rej, key=lambda x: abs(0.05 - x[2]), reverse=True)[:param_budget])

    for left, right, _ in to_redo:
        res = permutation_test(buckets[left], buckets[right])
        if res[3].startswith("Reject"):
            if res[1] > 0:
                final.append((left, right, -1))
            else:
                final.append((right, left, -1))

    print("NB de comparaisons significatives (welch + X param)", len(final))
     #print_comp_list(sorted(final, key=lambda x: x[0] + x[1]))

    # borda hypothesis
    patrick_format = [(a, b, 1, None, None) for (a, b, c) in final]
    hypothesis = computeRanksForAll(patrick_format, adom).items()

    print("Hypothesis as predicted: ", hypothesis)
    limitedHyp = []
    valsToSelect = []
    j = 0
    for h in hypothesis:
        if (h[1] <= nbAdomVals and j < nbAdomVals):
            limitedHyp.append(h)
            valsToSelect.append(h[0])
            j = j + 1
    print("Hypothesis limited: ", limitedHyp)
    print("vals: ",valsToSelect)

    emptyGBresult, emptyGBresultAll = emptyGB(conn, nbAdomVals, table, sel, meas)
    print("Empty GB says:", emptyGBresult)

    proba=0.1
    error=0.3
    # without replacement!
    sizeofsample=int(bernstein.sizeOfSampleHoeffding(proba ,error))+1
    print(sizeofsample)
    pwrset=dbStuff.getCuboidsOfAtt(groupbyAtt,sel)
    print(str(tuple(valsToSelect)))
    queryCountviolations,queryCountCuboid=bernstein.getSample(proba, error, pwrset, sel, measBase, function, table, tuple(valsToSelect), limitedHyp)
    #queryCountviolations,queryCountCuboid=bernstein.getSample(proba, error, pwrset, sel, measBase, function, table, tuple(valsToSelect), emptyGBresult)

    for i in range(len(queryCountviolations)):
        print(queryCountviolations[i])
        print(queryCountCuboid[i])
        v=dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
        c=dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]
        print(v)
        print(c)
        print(v/c)



    ''' 
    if conn:

        tabHypo = []
        resultRuns = []
        for i in range(nbruns):
            rankingFromPairwise.pairwiseComparison = []

            start_time = time.time()
            #generate hypothesis: ordering of members such that mean is greater (statistically significant on sample)
            hypothesis = generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, samplingMethod)

            # below some testing stuff
            print("Hypothesis as predicted: ", hypothesis)

            dbStuff.dropAllMVs(conn)
            dbStuff.createMV(conn, groupbyAtt, sel, meas, table, 0.5)
            tabView=dbStuff.getMVnames(conn)
            n, nbV=estimateViolations(conn, meas, measBase, table, sel, tabView, hypothesis)

            print("on " + str(n) + " draws, there are " + str(nbV) + " violations")
            print("violation rate is: ", nbV/n)
            
            
            
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
            emptyGBresult, emptyGBresultAll = emptyGB(conn, nbAdomVals, table, sel, meas)
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

        


        if not DEBUG_FLAG:
            names = ['tau', 'expected', 'time']
            title = 'Sample size=' + str(sampleSize)
            plot_curves(resultRuns, names, 'time', 'expected', title)

      
        # Close the connection
        close_connection(conn)
    '''

