import random
import math
import statistics
import numpy as np
import dbStuff
import statStuff
import utilities
from plotStuff import plot_curves_with_error_bars
from dbStuff import execute_query, connect_to_db, close_connection, getSample
from statStuff import permutation_test, compute_skewness, claireStat
import time
from rankingFromPairwise import computeRanksForAll, generateHypothesisTest

from statsmodels.stats.multitest import fdrcorrection


import configparser
import json

import bernstein


# ------  Debug ?  ------------
DEBUG_FLAG = True


def countViolations(conn,query,hypothesis):
    #print(query)
    hyp=[a for (a,b) in hypothesis]
    #print('hyp:',hyp)
    v=0
    res=dbStuff.execute_query(conn,query)
    for r in res:
        #print(r)
        #print('this is s', r[1])
        s=r[1].split(",")
        #print('this is s', s)
        tau=statStuff.compute_kendall_tau(s,hyp)[0]
        #print('tau:',tau)
        v=v+( (1-tau) *((len(hyp)*(len(hyp)+1))/2) )
    return v

def get_state_sample(conn, measBase, table, sel, sampleSize, state):

    querySample = "SELECT "+sel+", "+measBase+" FROM "+table+" where "+sel+" = '"+str(state)+"' limit "+str(sampleSize)+";"
    #print('stat query:', querySample)
    resultVals = execute_query(conn, querySample)
    return resultVals


def fetchCongressionalSample(conn,sel,table,measBase,sampleSize, adom_restr=None):
    # fetch the congressional sample
    if adom_restr:
        adom = adom_restr
    else:
        adom = [x[0] for x in execute_query(conn, "select distinct  " + sel + " from " + table + ";")]
    table_size = execute_query(conn, "select count(1) from " + table + ";")[0][0]

    sample_size = int(table_size * sampleSize)
    alpha = 0.10
    house_size = sample_size * alpha
    senate_size = sample_size * (1 - alpha)

    house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)

    senate = []
    state_sample_size = int(senate_size / len(adom))
    for state in adom:
        senate.extend(get_state_sample(conn, measBase, table, sel, state_sample_size, state))

    if adom_restr:
        house = list(filter(lambda x: x[0] in adom_restr, house))
    congress = house + senate
    # END - fetch the congressional sample
    return adom, congress


def getHypothesisCongressionalSampling(adom,congress):

    buckets = {s: [] for s in adom}
    skews = dict()
    for item in congress:
        buckets[item[0]].append(item[1])
    for k in buckets.keys():
        skews[k] = compute_skewness(buckets[k])

    # do all welch tests
    param_budget = 20
    param_budget = int(param_budget / 2)

    from scipy.stats import ttest_ind

    raw_comparisons = []

    for i in range(len(adom)):
        for j in range(i + 1, len(adom)):
            left = adom[i]
            right = adom[j]
            res = ttest_ind(buckets[left], buckets[right], equal_var=False)
            stat_c = claireStat(skews[left], skews[right], len(buckets[left]), len(buckets[right]))
            if res.statistic < 0:
                raw_comparisons.append((left, right, stat_c, res.pvalue ))
            else:
                raw_comparisons.append((right, left, stat_c, res.pvalue ))

    w_comparisons = []
    w_comparisons_rej = []
    #print(raw_comparisons)
    rejected, corrected = fdrcorrection([x[3] for x in raw_comparisons], alpha=0.05)
    for i in range(len(raw_comparisons)):
        if rejected[i]:
            w_comparisons_rej.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))
        else:
            w_comparisons.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))

    print("NB de comparaisons significatives (welch)", len(w_comparisons))
    # print_comp_list(sorted(w_comparisons, key=lambda x: x[0] + x[1]))
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
    # print_comp_list(sorted(final, key=lambda x: x[0] + x[1]))

    # borda hypothesis
    patrick_format = [(a, b, 1, None, None) for (a, b, c) in final]
    #print('alex pairwise:',patrick_format)
    hypothesis = computeRanksForAll(patrick_format, adom).items()
    #print(hypothesis)
    hypothesis = sorted(
        hypothesis,
        key=lambda x: x[1],
        reverse=True
    )
    #print(hypothesis)
    #hypothesis = [(a, b + 1) for (a, b) in hypothesis]
    correctHyp=[]
    i=1
    prevb=-1
    for (a,b) in hypothesis:
        if prevb==-1:
            currentRank = 1
            correctHyp.append((a,currentRank))
            prevb=b
        else:
            if b==prevb:
                correctHyp.append((a, currentRank))
            else:
                currentRank=currentRank+1
                correctHyp.append((a, currentRank))
                prevb=b

    print('Alex hypothesis:',correctHyp)
    # checking all comparisons
    valsToSelect=('HA','00','NK')
    correctHyp=generateHypothesisTest(conn, meas, measBase, table, sel, 9307,'SYSTEM_ROWS',valsToSelect)
    print('all comp. hypothesis:',correctHyp)
    return correctHyp


def test(conn, nbAdomVals, ratioViolations, proba, error, percentOfLattice, groupbyAtt, sel, measBase, function,table,sampleSize,comparison=False,generateIndex=False):
    #sampling
    start_time = time.time()
    adom, congress=fetchCongressionalSample(conn,sel,table,measBase,sampleSize, adom_restr=prefs)
    end_time = time.time()
    samplingTime = end_time - start_time
    print('sampling time:',samplingTime)

    # compute hypothesis
    start_time = time.time()
    hypothesis = getHypothesisCongressionalSampling(adom,congress)
    end_time = time.time()
    hypothesisGenerationTime = end_time - start_time
    print('hypothesis generation time:', hypothesisGenerationTime)

    print("Hypothesis as predicted: ", hypothesis)
    limitedHyp = []
    valsToSelect = []
    j = 0
    for h in hypothesis:
        if (h[1] <= nbAdomVals and j < nbAdomVals):
            limitedHyp.append(h)
            valsToSelect.append(h[0])
            j = j + 1
    #print("Hypothesis limited to choosen values: ", limitedHyp)

    # print("vals: ",valsToSelect)

    # just for checking on groupBy sel
    # emptyGBresult, emptyGBresultAll = emptyGB(conn, nbAdomVals, table, sel, meas)
    # print("Empty GB says:", emptyGBresult)
    # valsEmptyGB=[a for (a, b) in emptyGBresult]
    # print(valsEmptyGB)


    # generate and get all materialized cuboids
    dbStuff.dropAllMVs(conn)
    dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice)
    mvnames = dbStuff.getMVnames(conn)

    # generate hash index on sel attribute
    if generateIndex==True:
        dbStuff.generateHashIndex(conn,table,sel)

    #validation of hypothesis
    start_time = time.time()

    sizeofsample = int(bernstein.sizeOfSampleHoeffding(proba, error)) + 1
    print('size of sample according to Hoeffding:', sizeofsample)

    # total number of cuboids
    N = len(utilities.powerset(groupbyAtt))
    print('size of sample according to Bardenet:',
          int(bernstein.sizeOfSampleHoeffdingSerflingFromBardenet(proba, error, N)) + 1)

    pwrset = dbStuff.getCuboidsOfAtt(groupbyAtt, sel)
    print(str(tuple(valsToSelect)))
    ranks, queryCountviolations, queryCountCuboid, cuboid = bernstein.getSample(proba, error, pwrset, sel, measBase, function,
                                                                         table, tuple(valsToSelect), limitedHyp,
                                                                         mvnames,False,True)
    # queryCountviolations, queryCountCuboid, cuboid=bernstein.getSample(proba, error, pwrset, sel, measBase, function, table, tuple(valsEmptyGB), emptyGBresult, mvnames)

    tabRandomVar = []
    nbViewOK = 0
    for i in range(len(queryCountviolations)):
        #print(queryCountviolations[i])
        #print(queryCountCuboid[i])
        #print(ranks[i])
        v=countViolations(conn,ranks[i],hypothesis)
        #v = dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
        c = dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]
        # print(v)
        print(c)
        if c!=0:
            print(v/c, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)
            if v / c < ratioViolations:
                tabRandomVar.append(1)
                nbViewOK = nbViewOK + 1
            else:
                tabRandomVar.append(0)
        else:
            print("inconclusive, not enough tuples in cuboid for select values")


    end_time = time.time()
    validationTime = end_time - start_time
    print('validation time:', validationTime)

    variance = np.var(tabRandomVar)
    # print('variance: ', variance)
    prediction = nbViewOK / sizeofsample
    predictionNbOk = prediction * len(pwrset)
    print('nb of views ok: ', nbViewOK, 'out of ', sizeofsample, 'views, i.e., rate of:', nbViewOK / sizeofsample)
    print('predicting number of views ok:', predictionNbOk)

    nbErrors = 2
    print('probability of making ', nbErrors, ' errors: ', bernstein.bernsteinBound(variance, nbErrors))
    print('the error (according to Bernstein) for sum and confidence interval of size', proba, ' is: ',
          bernstein.bersteinError(proba, variance))
    bennetError=bernstein.bennetErrorOnAvg(proba, variance, sizeofsample)
    print('the error (according to Bennet) for avg and confidence interval of size', proba, ' is: ',
          bernstein.bennetErrorOnAvg(proba, variance, sizeofsample))
    print('the error (empirical bennet) for avg and confidence interval of size', proba, ' is: ',
          bernstein.empiricalBennetFromMaurer(proba, variance, sizeofsample))
    print('the error (according to bardenet) for avg and confidence interval of size', proba, ' is: ',
          bernstein.empiricalBernsteinFromBardenet(proba, variance, sizeofsample, N))



    if comparison==True:
        # comparison with ground truth
        dbStuff.dropAllMVs(conn)
        nbMVs = dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, 1)

        ranks, queryCountviolations, queryCountCuboid, cuboid = bernstein.generateAllqueries(pwrset, sel, measBase, function,
                                                                                      table, tuple(valsToSelect),
                                                                                      limitedHyp, mvnames)

        tabRandomVar = []
        nbViewOK = 0
        for i in range(len(queryCountviolations)):
            #print('gt violations:',queryCountviolations[i])
            #print('gt count:',queryCountCuboid[i])
            #v = dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
            v = countViolations(conn, ranks[i], hypothesis)
            c = dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]
            # print(v)
            # print(c)
            if c != 0:
                print(v / c, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)
                if v / c < ratioViolations:
                    tabRandomVar.append(1)
                    nbViewOK = nbViewOK + 1
                else:
                    tabRandomVar.append(0)
            else:
                print("inconclusive, not enough tuples in cuboid for select values")

        variance = np.var(tabRandomVar)
        # print('variance: ', variance)
        print('*** comparison to ground truth ***')
        print('nb of views ok: ', nbViewOK, 'out of ', nbMVs, 'views, i.e., rate of:', nbViewOK / nbMVs)
        gtratio= nbViewOK / nbMVs

        realError=abs(prediction - (nbViewOK / nbMVs))
        print('Error on avg is: ', abs(prediction - (nbViewOK / nbMVs)))

        print('Error on sum is: ', abs(nbViewOK - predictionNbOk))

        print('the error (according to Bennet) for avg and confidence interval of size', proba, ' is: ',
              bernstein.bennetErrorOnAvg(proba, variance, sizeofsample))
        print('the error (according to Bernstein) for confidence interval of size', proba, ' is: ',
              bernstein.bersteinError(proba, variance))

        return prediction,bennetError,realError,gtratio
    else:
        return bennetError, samplingTime, hypothesisGenerationTime, validationTime


# TODO
#
#  change sel and measures
#  use other databases
#

if __name__ == "__main__":

    config = configparser.ConfigParser()

    # The DB wee want
    config.read('configs/flights.ini')
    #config.read('configs/ssb.ini')
    # The system this is running on
    USER = "PM"

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
    prefs = json.loads(config.get("Common", "preferred"))
    if len(prefs) == 0:
        prefs = None


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
    sampleSize = 1
    samplingMethod = 'SYSTEM_ROWS'  # or SYSTEM

    if DEBUG_FLAG:
        nbruns = 1
    else:
        nbruns = 10

    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    # to always have the same order in group bys, with sel attribute last
    groupbyAtt.sort()

    ratioViolations = 0.4
    ratioCuboidOK = 0.8

    proba = 0.2
    error = 0.4  # rate

    percentOfLattice=0.3

    nbWrongRanking=0
    resultRuns=[]

    # do we compare to ground truth?
    comparison = True

    nbOfRuns=10

    data=[]

    if comparison==True:

        listPred=[]
        devPred=[]
        listError=[]
        devError=[]
        listWR=[]
        devWR=[]
        listBennet=[]
        devBennet=[]

        paramTested = 'Percent of lattice'
        tabTest=(0.1, 0.25, 0.5, 0.75, 1)

        for percentOfLattice in tabTest:
        #for sampleSize in (0.1, 0.25, 0.5, 0.75, 1):
        #for nbAdomVals in range(2,10):

            predictionTab=[]
            realErrorTab=[]
            nbWrongRankingTab=[]
            bennetTab = []

            for i in range(nbOfRuns):

                prediction,bennetError,realError,gtratio=test(conn, nbAdomVals, ratioViolations, proba, error, percentOfLattice, groupbyAtt, sel, measBase, function,table, sampleSize, comparison)
                #resultRuns.append((percentOfLattice,prediction,bennetError,realError))

                predictionTab.append(prediction)
                bennetTab.append(bennetError)
                realErrorTab.append(realError)
                if gtratio <ratioCuboidOK:
                    nbWrongRanking=1
                else:
                    nbWrongRanking = 0
                nbWrongRankingTab.append(nbWrongRanking)

            meanPred=statistics.mean(predictionTab)
            stdevPred = statistics.stdev(predictionTab)
            meanBen= statistics.mean(bennetTab)
            stdevBen = statistics.stdev(bennetTab)
            meanError=statistics.mean(realErrorTab)
            stdevError = statistics.stdev(realErrorTab)
            meanWRTab = statistics.mean(nbWrongRankingTab)
            stdevWRTab = statistics.stdev(nbWrongRankingTab)


            listPred.append(meanPred)
            devPred.append(stdevPred)
            listBennet.append(meanBen)
            devBennet.append(stdevBen)
            listError.append(meanError)
            devError.append(stdevError)
            listWR.append(meanWRTab)
            devWR.append(stdevWRTab)

        # Example usage:
        data = [
            {'x': tabTest, 'y':listPred,  'yerr': devPred, 'label': 'prediction'},
            {'x': tabTest, 'y': listError, 'yerr': devError, 'label': 'error'},
            {'x': tabTest, 'y': listWR, 'yerr': devWR, 'label': 'wrong prediction'},
            {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Bennet error'}
        ]

        plot_curves_with_error_bars(data, x_label=paramTested, y_label='Error',
                                    title='prediction and errors')
        #print('Number of incorrect hypothesis:', nbWrongRanking)
        #names = ['prediction', 'bennet', 'error']
        #title = 'top-' + str(nbAdomVals)
        #plot_curves(resultRuns, names, 'percentoflattice', 'error', title)
    else:

        listBennet = []
        devBennet = []
        listSampling=[]
        devSampling=[]
        listHypo=[]
        devHypo=[]
        listValid=[]
        devValid=[]

        tabTest=(0.1, 0.2, 0.3, 0.4, 0.5)
        paramTested='Percent of Lattice'
        #paramTested='Sample size'

        for percentOfLattice in tabTest:
        #for sampleSize in tabTest:
        # for nbAdomVals in range(2,10):

            benTab=[]
            samplingTab=[]
            hypoTab=[]
            validTab = []

            for i in range(nbOfRuns):
                bennetError, samplingTime, hypothesisTime, validationTime = test(conn, nbAdomVals, ratioViolations, proba, error,
                                                               percentOfLattice, groupbyAtt, sel, measBase, function,
                                                               table, sampleSize, comparison)
                benTab.append(bennetError)
                samplingTab.append(samplingTime)
                hypoTab.append(hypothesisTime)
                validTab.append(validationTime)

            #resultRuns.append((percentOfLattice, bennetError, hypothesisTime, validationTime))
            meanBen = statistics.mean(benTab)
            stdevBen = statistics.stdev(benTab)
            meanSamp = statistics.mean(samplingTab)
            stdevSamp = statistics.stdev(samplingTab)
            meanHypo = statistics.mean(hypoTab)
            stdevHypo = statistics.stdev(hypoTab)
            meanValid = statistics.mean(validTab)
            stdevValid = statistics.stdev(validTab)

            listBennet.append(meanBen)
            devBennet.append(stdevBen)
            listSampling.append(meanSamp)
            devSampling.append(stdevSamp)
            listHypo.append(meanHypo)
            devHypo.append(stdevHypo)
            listValid.append(meanValid)
            devValid.append(stdevValid)

        data = [
            {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Bennet error'},
            {'x': tabTest, 'y': listSampling, 'yerr': devSampling, 'label': 'Sampling time'},
            {'x': tabTest, 'y': listHypo, 'yerr': devHypo, 'label': 'Hypothesis time'},
            {'x': tabTest, 'y': listValid, 'yerr': devValid, 'label': 'Generation time'}
        ]

        plot_curves_with_error_bars(data, x_label=paramTested, y_label='Time (s)',
                                    title='Times')
        #names = ['error', 'hypothesis', 'validation']
        #title = 'top-' + str(nbAdomVals)
        #plot_curves(resultRuns, names, 'percentoflattice', 'time', title)



    # Close the connection
    close_connection(conn)



