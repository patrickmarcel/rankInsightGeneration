import random
import math
import numpy as np
from copy import copy

import dbStuff
import RankingFromPairwise
import utilities
from utilities import powerset
from plotStuff import plot_curves
from dbStuff import execute_query, connect_to_db, close_connection, getSample, emptyGB
from statStuff import welch_ttest, permutation_test, compute_skewness, compute_kendall_tau, benjamini_hochberg, \
    benjamini_hochberg_statmod, claireStat
import time
from RankingFromPairwise import computeRanksForAll

from main import generateComparisonsWithMergeSort

from statsmodels.stats.multitest import fdrcorrection



import configparser
import json


def get_state_sample(conn, measBase, table, sel, sampleSize, state):

    querySample = "SELECT "+sel+", "+measBase+" FROM "+table+" where "+sel+" = '"+state+"' limit "+str(sampleSize)+";"
    #print(querySample)
    resultVals = execute_query(conn, querySample)
    return resultVals

def generateHypothesisTest_from_sample(conn, meas, measBase, table, sel, sample):

    # get adom values
    Sels = list(set([x[0] for x in sample]))

    buckets = {s: [] for s in Sels}
    for item in sample:
        buckets[item[0]].append(item[1])

    #analyse sample for each adom value: value, nb of measures, skewness, and tuples
    S = []
    for v in Sels:
        data = buckets[v]
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

def print_comp_list(l):
    if len(l[0]) == 2:
        l = [str(a) + ">" + str(b) for a, b in l]
    else:
        l = [str(a) + ">" + str(b) for a, b, c in l]
    print(" | ".join(l))

def strip_comp_list(l):
    if len(l[0]) == 2:
        l = [str(a) + ">" + str(b) for a, b in l]
    else:
        l = [str(a) + ">" + str(b) for a, b, c in l]
    return l

def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(set(list1)) + len(set(list2))) - intersection
    return float(intersection) / union


if __name__ == "__main__":

    config = configparser.ConfigParser()
    # The DB wee want
    config.read('configs/flights.ini')
    # The system this is running on
    USER = "AC"
    SAMPLE_SIZE_REL = 0.2

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

    conn = connect_to_db(dbname, user, password, host, port)

    # collect some stats
    table_size = execute_query(conn, "select count(1) from " + table + ";")[0][0]
    adom = [x[0] for x in execute_query(conn, "select distinct  "+sel+" from "+table+";")]
    print("NB de groupes",len(adom))

    # fetch the congressional sample
    sample_size = int(table_size*SAMPLE_SIZE_REL)
    alpha = 0.25
    house_size = sample_size*alpha
    senate_size = sample_size * (1-alpha)

    house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)

    senate = []
    state_sample_size = int(senate_size/len(adom))
    for state in adom:
        senate.extend(get_state_sample(conn, measBase, table, sel, state_sample_size, state))

    congress = house + senate
    # END - fetch the congressional sample

    buckets = {s: [] for s in adom}
    skews = dict()
    for item in congress:
        buckets[item[0]].append(item[1])
    for k in buckets.keys():
        skews[k] = compute_skewness(buckets[k])


    # do all welch tests
    param_budget = 20
    param_budget = int(param_budget/2)

    from scipy.stats import ttest_ind

    raw_comparisons = []

    for i in range(len(adom)):
        for j in range(i + 1, len(adom)):
            left = adom[i]
            right = adom[j]
            res = ttest_ind(buckets[left], buckets[right], equal_var=False)
            stat_c = claireStat(skews[left], skews[right], len(left), len(right))
            if res.statistic < 0:
                raw_comparisons.append((left, right, stat_c, res.pvalue ))
            else:
                raw_comparisons.append((right, left, stat_c, res.pvalue ))

    w_comparisons = []
    w_comparisons_rej = []
    print(raw_comparisons)
    rejected, corrected = fdrcorrection([x[3] for x in raw_comparisons], alpha=0.05)
    for i in range(len(raw_comparisons)):
        if rejected[i]:
            w_comparisons_rej.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))
        else:
            w_comparisons.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))

    print("NB de comparaisons significatives (welch)", len(w_comparisons))
    print_comp_list(sorted(w_comparisons,key=lambda x : x[0]+x[1]))
    by_prox_to_threshold = sorted(w_comparisons, key=lambda x: abs(0.05-x[2]), reverse=True)
    #print(by_prox_to_threshold)

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
    print_comp_list(sorted(final, key=lambda x: x[0]+x[1]))
    """
    # get actual ranking for the sample
    comparisons = [] # (a , b) means a greater than b
    # Exhaustive method with non-parametric
    matrix = [[1 for j in adom] for i in adom]
    for i in range(len(adom)):
        for j in range(i  + 1):
            left = adom[i]
            right = adom[j]
            res = permutation_test(buckets[left], buckets[right])
            matrix[i][j] = res[1]
            matrix[j][i] = res[1]
            if res[3].startswith("Reject"):
                if res[1] > 0:
                    comparisons.append((left, right))
                else:
                    comparisons.append((right, left))

    print("NB de comparaisons significatives (exhaustif)", len(comparisons))
    print_comp_list(sorted(comparisons,key=lambda x : x[0]+x[1]))

    print("Only welch", jaccard_similarity(strip_comp_list(comparisons), strip_comp_list(w_comparisons)))
    print("Welch + param", jaccard_similarity(strip_comp_list(comparisons), strip_comp_list(final)))
    """
    #borda hypothesis
    patrick_format = [(a,b, 1, None, None) for (a, b, c) in final]
    other_patrick_format = computeRanksForAll(patrick_format, adom).items()