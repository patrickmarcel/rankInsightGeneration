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

def  compareHypToGB(hypothesis, conn, measBase,function, sel, vals):
    #query="Select " + sel +" from  " + sel + " where " + sel + " in " +  str(vals) + " order by " + measBase + " desc;"
    query="select 'all',string_agg(" + sel + ",',') from (SELECT " + sel + ", " + function + "(" + measBase + "),  rank () over (  order by " + function + "(" + measBase + ") desc ) as rank FROM \"" + sel + "\" WHERE " + sel + " in " + str(vals) +" group by " +  sel + " order by rank);"
    #print(query)
    v,ratio=countViolations(conn,query,hypothesis)
    #print(v)
    print("hypothesis compared to group by ",sel," has ",v," violations")
    return v


def countViolations(conn,query,hypothesis):
    #print(query)
    hyp=[a for (a,b) in hypothesis]
    #print('hyp:',hyp)
    v=0
    res=dbStuff.execute_query(conn,query)
    normalize=0
    for r in res:
        #print(r)
        #print('this is s', r[1])
        s=r[-1].split(",")
        #print('this is s', s)
        normalize=normalize+(len(s) * ( len(s) -1 ))/2
        if len(s)==len(hyp):
            #tau=statStuff.compute_kendall_tau(s,hyp)[0]
            #if tau!=1:
            #    v=v+1
            tau=statStuff.normalised_kendall_tau_distance(s,hyp)
            #print('tau:',tau)
            v=v+tau
        else:
            # s is smaller
            hyp2=hyp.copy()
            #print('hyp2:',hyp2,' and s:',s)
            for e in hyp:
                #print('e:',e)
                if e not in s:
                    hyp2.remove(e)
            #print('hyp2:',hyp2)
            #tau = statStuff.compute_kendall_tau(s, hyp2)[0]
            #if tau != 1:
            #    v = v + 1
            #print("s:",s)
            #print("hyp2:",hyp2)
            tau = statStuff.normalised_kendall_tau_distance(s, hyp2)
            #print('tau:',tau)
            v = v + tau
    #if len(res)!=0:
    #    ratio=v/len(res)
    #else:
    #    ratio=0
    if len(res) != 0:
        ratio=v/normalize
    else:
        ratio=0
    return v,ratio



def getHypothesisAllComparisons(conn, meas, measBase, table, sel,valsToSelect, sampleSize, method='SYSTEM_ROWS'):
    # checking all comparisons
    correctHyp = generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, method, valsToSelect)
    #print('all comp. hypothesis:', correctHyp)
    return correctHyp

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

    #print('Alex hypothesis:',correctHyp)
    return correctHyp


def test(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice, groupbyAtt, sel, measBase, function,table,sampleSize,comparison=False,generateIndex=False,allComparison=True):

    if allComparison==False:
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
    else:
        # sampling and hypothesis
        start_time = time.time()
        hypothesis=getHypothesisAllComparisons(conn, meas, measBase, table, sel, tuple(prefs), sampleSize, method='SYSTEM_ROWS')
        end_time = time.time()
        samplingTime = end_time - start_time
        hypothesisGenerationTime= samplingTime
        #print('sampling time:', samplingTime)
        print('hypothesis generation time:', hypothesisGenerationTime)

    print("Hypothesis predicted: ", hypothesis)


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

    # generate hash index on sel attribute
    if generateIndex == True:
        dbStuff.generateHashIndex(conn, table, sel)
    else:
        dbStuff.dropIndex(conn, table, sel)
    # generate and get all materialized cuboids
    dbStuff.dropAllMVs(conn)
    dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,generateIndex)
    mvnames = dbStuff.getMVnames(conn)



    #validation of hypothesis
    start_time = time.time()

    sizeofsample = int(bernstein.sizeOfSampleHoeffding(proba, error)) + 1
    print('size of sample according to Hoeffding:', sizeofsample)

    # total number of cuboids
    N = len(utilities.powerset(groupbyAtt))
    #print('size of sample according to Bardenet:',
    #      int(bernstein.sizeOfSampleHoeffdingSerflingFromBardenet(proba, error, N)) + 1)

    pwrset = dbStuff.getCuboidsOfAtt(groupbyAtt, sel)
    #print(str(tuple(valsToSelect)))

    print("Generating sample of aggregate queries")
    ranks, queryCountviolations, queryCountCuboid, cuboid = bernstein.getSample(proba, error, pwrset, sel, measBase, function,
                                                                         table, tuple(valsToSelect), limitedHyp,
                                                                         mvnames,False,True)
    # queryCountviolations, queryCountCuboid, cuboid=bernstein.getSample(proba, error, pwrset, sel, measBase, function, table, tuple(valsEmptyGB), emptyGBresult, mvnames)


    print("Computing violations")

    tabRandomVar = []
    nbViewOK = 0

    nbInconclusive=0

    for i in range(len(queryCountviolations)):
        #print(queryCountviolations[i])
        #print(queryCountCuboid[i])
        #print(ranks[i])
        v,ratio=countViolations(conn,ranks[i],hypothesis)
        #v = dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
        c = dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]
        # print(v)
        #print(c)
        if c!=0:
            #OLD print(v/c, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)

            #print(ratio, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)

            if ratio < ratioViolations:
                tabRandomVar.append(1)
                nbViewOK = nbViewOK + 1
            else:
                tabRandomVar.append(0)
        else:
            #print("inconclusive, not enough tuples in cuboid for select values")
            sizeofsample = sizeofsample - 1
            nbInconclusive=nbInconclusive+1


    end_time = time.time()
    validationTime = end_time - start_time
    print('validation time:', validationTime)

    print('number of inconclusive: ',nbInconclusive, ' ratio: ',nbInconclusive/len(queryCountviolations))

    variance = np.var(tabRandomVar)
    # print('variance: ', variance)
    prediction = nbViewOK / sizeofsample
    predictionNbOk = prediction * len(pwrset)
    print('nb of views ok: ', nbViewOK, 'out of ', sizeofsample, 'views, i.e., rate of:', nbViewOK / sizeofsample)
    print('predicting number of views ok:', predictionNbOk)

    nbErrors = 2
    #print('probability of making ', nbErrors, ' errors: ', bernstein.bernsteinBound(variance, nbErrors))
    #print('the error (according to Bernstein) for sum and confidence interval of size', proba, ' is: ',
    #      bernstein.bersteinError(proba, variance))
    bennetError=bernstein.bennetErrorOnAvg(proba, variance, sizeofsample)
    print('the error (according to Bennet) for avg and confidence interval of size', proba, ' is: ',
          bernstein.bennetErrorOnAvg(proba, variance, sizeofsample))
    print('the error (empirical bennet) for avg and confidence interval of size', proba, ' is: ',
          bernstein.empiricalBennetFromMaurer(proba, variance, sizeofsample))
    ###
    ### IF REPORTING EMPIRICAL ERROR
    ### UNCOMMENT BELOW
    #bennetError = bernstein.empiricalBennetFromMaurer(proba, variance, sizeofsample)

    #print('the error (according to bardenet) for avg and confidence interval of size', proba, ' is: ',
    #      bernstein.empiricalBernsteinFromBardenet(proba, variance, sizeofsample, N))



    if comparison==True:
        # comparison with ground truth
        print('*** comparison to ground truth ***')

        dbStuff.dropAllMVs(conn)
        nbMVs = dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, 1,generateIndex)
        #print("nb views generated:",nbMVs)

        compareHypToGB(hypothesis, conn, measBase,function, sel,tuple(valsToSelect))
        ranks, queryCountviolations, queryCountCuboid, cuboid = bernstein.generateAllqueries(pwrset, sel, measBase, function,
                                                                                      table, tuple(valsToSelect),
                                                                                      limitedHyp, mvnames)
        nbInconclusive=0
        tabRandomVar = []
        nbViewOK = 0
        for i in range(len(queryCountviolations)):
            #print('gt violations:',queryCountviolations[i])
            #print('gt count:',queryCountCuboid[i])
            #v = dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
            v,ratio = countViolations(conn, ranks[i], hypothesis)
            c = dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]
            # print(v)
            # print(c)
            if c != 0:
                #OLD print(v / c, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)

                #print(ratio, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)
                if ratio < ratioViolations:
                    tabRandomVar.append(1)
                    nbViewOK = nbViewOK + 1
                else:
                    tabRandomVar.append(0)
            else:
                #print("inconclusive, not enough tuples in cuboid for select values")
                nbMVs=nbMVs-1
                nbInconclusive=nbInconclusive+1

        variance = np.var(tabRandomVar)
        # print('variance: ', variance)

        print('number of inconclusive: ', nbInconclusive, ' ratio: ', nbInconclusive / len(queryCountviolations))

        print('nb of views ok: ', nbViewOK, 'out of ', nbMVs, 'views, i.e., rate of:', nbViewOK / nbMVs)
        gtratio= nbViewOK / nbMVs

        realError=abs(prediction - (nbViewOK / nbMVs))
        print('Error is: ', abs(prediction - (nbViewOK / nbMVs)))

        #print('Error on sum is: ', abs(nbViewOK - predictionNbOk))

        #print('the error (according to Bennet) for avg and confidence interval of size', proba, ' is: ',
              #bernstein.bennetErrorOnAvg(proba, variance, sizeofsample))
        #print('the error (according to Bernstein) for confidence interval of size', proba, ' is: ',
              #bernstein.bersteinError(proba, variance))

        return prediction,bennetError,realError,gtratio
    else:
        return bennetError, samplingTime, hypothesisGenerationTime, validationTime


# TODO
#


if __name__ == "__main__":

    config = configparser.ConfigParser()

    # The DB wee want
    #config.read('configs/flights.ini')
    config.read('configs/ssb.ini')
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
    #print(tuple(prefs))
    if len(prefs) == 0:
        prefs = None

    if DEBUG_FLAG:
        nbruns = 1
    else:
        nbruns = 10


    # for Hoeffding - OLD
    epsilon = 0.1
    alpha = 0.1
    p = 0
    H = []
    threshold = 0.1  # 10% of tuples violating the order
    n = math.log(2 / alpha, 10) / (2* epsilon * epsilon)
    n = math.ceil(n)

    ###
    ### PARAMETERS
    ###
    proba = 0.1
    error = 0.4  # rate

    # number of values of adom to consider - default = all of prefs
    nbAdomVals = len(prefs)

    # for DB sampling
    sampleSize = 1
    samplingMethod = 'SYSTEM_ROWS'  # or SYSTEM

    # ratio max of violations in a cuboid
    ratioViolations = 0.4

    # ratio min of cuboids with raio violations < ratioViolations
    ratioCuboidOK = 0.8

    # percentage of the lattice to generate
    percentOfLattice = 0.3

    # do we generate indexes?
    generateIndex = False

    # do we compare to ground truth?
    comparison = False

    #do we generate all comparisons?
    allComparisons = True

    #number of runs
    nbOfRuns = 10


    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    #get size of R
    sizeOfR=dbStuff.getSizeOf(conn,table)
    #print(sizeOfR)

    # to always have the same order in group bys, with sel attribute last
    groupbyAtt.sort()

    nbWrongRanking=0
    resultRuns=[]
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

        paramTested = 'Sample size'
        #paramTested = 'Percent of lattice'
        #tabTest=(0.1, 0.25, 0.5, 0.75, 1)
        tabTest=(0.1, 0.25, 0.5, 0.75, 1)


        #for percentOfLattice in tabTest:
        for sampleSize in tabTest:
        #for nbAdomVals in range(2,10):
            sampleSize=sampleSize*sizeOfR

            print("--- TESTING VALUE:",sampleSize)
            predictionTab=[]
            realErrorTab=[]
            nbWrongRankingTab=[]
            bennetTab = []

            for i in range(nbOfRuns):

                print("-----RUN: ",i)
                prediction,bennetError,realError,gtratio=test(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice, groupbyAtt,
                                                              sel, measBase, function,table, sampleSize, comparison,generateIndex,allComparisons)
                #resultRuns.append((percentOfLattice,prediction,bennetError,realError))

                predictionTab.append(prediction)
                bennetTab.append(bennetError)
                realErrorTab.append(realError)
                print("Desired cuboid ratio is:",ratioCuboidOK,". We predicted ratio of: ",prediction,". Real ratio is: ",gtratio)
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
            {'x': tabTest, 'y': listError, 'yerr': devError, 'label': 'real error'},
            {'x': tabTest, 'y': listWR, 'yerr': devWR, 'label': 'unvalidated prediction'},
            {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Bennet theoretical error'}
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
        #paramTested='Percent of Lattice'
        paramTested='Sample size'

        #for percentOfLattice in tabTest:
        for sampleSize in tabTest:
        # for nbAdomVals in range(2,10):

            sampleSize = sampleSize * sizeOfR

            print("--- TESTING VALUE:", percentOfLattice)

            benTab=[]
            samplingTab=[]
            hypoTab=[]
            validTab = []

            for i in range(nbOfRuns):
                print("-----RUN: ",i)
                bennetError, samplingTime, hypothesisTime, validationTime = test(conn, nbAdomVals, prefs, ratioViolations, proba, error,
                                                               percentOfLattice, groupbyAtt, sel, measBase, function,
                                                               table, sampleSize, comparison,generateIndex,allComparisons)
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
            {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Bennet theoretical error'},
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



