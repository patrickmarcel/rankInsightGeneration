import random
import math
import numpy as np
from copy import copy

import dbStuff
import rankingFromPairwise
import utilities
from utilities import powerset
from plotStuff import plot_curves
from dbStuff import execute_query, connect_to_db, close_connection, getSample, emptyGB
from statStuff import welch_ttest, permutation_test, compute_skewness, compute_kendall_tau, benjamini_hochberg, \
    benjamini_hochberg_statmod, claireStat
import time
from rankingFromPairwise import computeRanksForAll

from main import generateComparisonsWithMergeSort


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


if __name__ == "__main__":

    config = configparser.ConfigParser()
    # The DB wee want
    config.read('configs/flights.ini')
    # The system this is running on
    USER = "AC"
    SAMPLE_SIZE_REL = 0.1

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

    # fetch the congressional sample
    sample_size = int(table_size*SAMPLE_SIZE_REL)
    alpha = 0.5
    house_size = sample_size*alpha
    senate_size = sample_size * (1-alpha)

    house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)

    senate = []
    state_sample_size = int(senate_size/len(adom))
    for state in adom:
        senate.extend(get_state_sample(conn, measBase, table, sel, state_sample_size, state))

    congress = house + senate

    # get actual ranking for the sample
    comparisons = [] # (a , b) means a greater than b

    buckets = {s: [] for s in adom}
    for item in congress:
        buckets[item[0]].append(item[1])

    #stats mathod


    # preference method
    for i in range(len(adom)):
        for j in range(i  + 1):
            left = adom[i]
            right = adom[j]
            res = permutation_test(buckets[left], buckets[right])
            if res[3].startswith("Reject"):
                if res[1] > 0:
                    comparisons.append((left, right))
                else:
                    comparisons.append((right, left))

    print(comparisons)
    from pwlistorder import agg_preferences, eval_ordering, minconflict, pagerank

    # creating the dictionary of aggregated preferences
    dct_prefs = agg_preferences(comparisons)

    starting_order = copy(adom)
    random.shuffle(starting_order)

    # running min-conflict, requires a starting point
    ordering_minconflict = minconflict(dct_prefs, starting_order)

    print("min conflict local search", ordering_minconflict)

    h = generateHypothesisTest_from_sample(conn, meas, measBase, table, sel, congress)

