import math
import statistics
import numpy as np
import time
import configparser
import json
from statsmodels.stats.multitest import fdrcorrection

import dbStuff
import plotStuff
import rankingFromPairwise
import statStuff
import utilities
from plotStuff import plot_curves_with_error_bars
from dbStuff import execute_query, connect_to_db, close_connection, getSample
from statStuff import permutation_test, compute_skewness, claireStat
from rankingFromPairwise import computeRanksForAll, generateHypothesisTest
import bounders
import tests
from tqdm import tqdm
import pandas as pd

# ------  Debug ?  ------------
DEBUG_FLAG = True

def generateHypothesisTestDolap(conn, meas, measBase, table, sel, sampleSize, method, valsToSelect=None):
    # sampling
    start_time = time.time()
    resultVals = getSample(conn, measBase, table, sel, sampleSize, method, False, valsToSelect)
    end_time = time.time()
    samplingTime = end_time - start_time
    #print('sampling time:', samplingTime)

    #resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=DEBUG_FLAG)
    #print(resultVals)

    start_time = time.time()

    # get adom values
    Sels = list(set([x[0] for x in resultVals]))
    #print('Sels:',Sels)
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

    #print('S:',S)

    # nlog(n) comparisons enough for recovering the true ranking when comparisons are certain (not noisy)
    # we should try less
    #print(len(Sels))
    #nbOfComparisons = len(Sels) * math.log(len(Sels), 2)
    #print("Number of comparisons to make: " + str(nbOfComparisons))

    pairwiseComparison=rankingFromPairwise.generateAllComparisons(Sels, S)

    #separation=computeSeparationJMLR18(pairwiseComparison,len(valsToSelect))

    #print("pairwise comparisons:")
    #for p in pairwiseComparison:
    #    print("p: ", p)

    #pairwiseComparison = generateComparisonsWithMergeSort(Sels, S)

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

    end_time = time.time()
    hypothesisGenerationTime = end_time - start_time
    #print('Hypothesis generation time:', hypothesisGenerationTime)
    if pairwiseComparison!=[]:
        pvalue=float(pairwiseComparison[0][4])
    else:
        pvalue=1000 #change me
    return hypothesis, samplingTime, hypothesisGenerationTime,pvalue


def  compareHypToGB(hypothesis, conn, measBase,function, sel, vals,mvnames, table):
    #query="Select " + sel +" from  " + sel + " where " + sel + " in " +  str(vals) + " order by " + measBase + " desc;"
    materialized=bounders.findMV(mvnames,sel,table)
    ##print("MATERIALIZED: ",materialized)
    query="select 'all',string_agg(" + sel + ",',') from (SELECT " + sel + ", " + function + "(" + measBase + "),  rank () over (  order by " + function + "(" + measBase + ") desc ) as rank FROM \"" + materialized + "\" WHERE " + sel + " in " + str(vals) +" group by " +  sel + " order by rank);"
    ##print(query)
    v,ratio,qtime=countViolations(conn,query,hypothesis)
    ##print(v)
    #print("hypothesis compared to group by ",sel," has ",v," violations")
    return v



def countViolationsDOLAP(conn,query,hypothesis):
    ##print(query)
    hyp=[str(a) for (a,b) in hypothesis]
    #print('hyp:',hyp)
    v=0
    res=dbStuff.execute_query(conn,query)
    normalize=0
    for r in res:
        #print(r)
        ##print('this is s', r[1])
        s=r[-1].split(",")
        #print('this is s', s)
        normalize=normalize+(len(s) * ( len(s) -1 ))/2
        if len(s)==len(hyp):

            tau,pvalue=statStuff.compute_kendall_tau(s,hyp)

            ##print('tau:',tau)
            tau=(tau+1)/2
            v=v+tau
        else:
            # s is smaller
            hyp2=hyp.copy()
            ##print('hyp2:',hyp2,' and s:',s)
            for e in hyp:
                ##print('e:',e)
                if e not in s:
                    hyp2.remove(e)

            tau,pvalue = statStuff.compute_kendall_tau(s, hyp2)

            ##print('tau:',tau)
            tau = (tau + 1) / 2
            v = v + tau

    if len(res) >1:
        ratio=v/normalize
    else:
        ratio=0
    return v,ratio,0



def countViolations(conn,query,hypothesis):
    ##print(query)
    hyp=[str(a) for (a,b) in hypothesis]
    ##print('hyp:',hyp)
    v=0
    start_time_q = time.time()
    res=dbStuff.execute_query(conn,query)
    end_time_q = time.time()
    querytime=end_time_q-start_time_q
    normalize=0
    for r in res:
        ##print(r)
        ##print('this is s', r[1])
        s=r[-1].split(",")
        ##print('this is s', s)
        normalize=normalize+(len(s) * ( len(s) -1 ))/2
        if len(s)==len(hyp):
            #tau=statStuff.compute_kendall_tau(s,hyp)[0]
            #if tau!=1:
            #    v=v+1
            #tau=statStuff.normalised_kendall_tau_distance(s,hyp)

            tau,pvalue=statStuff.compute_kendall_tau(s,hyp)

            ##print('tau:',tau)
            tau=(tau+1)/2
            v=v+tau
        else:
            # s is smaller
            hyp2=hyp.copy()
            ##print('hyp2:',hyp2,' and s:',s)
            for e in hyp:
                ##print('e:',e)
                if e not in s:
                    hyp2.remove(e)
            ##print('hyp2:',hyp2)
            #tau = statStuff.compute_kendall_tau(s, hyp2)[0]
            #if tau != 1:
            #    v = v + 1
            ##print("s:",s)
            ##print("hyp2:",hyp2)
            #tau = statStuff.normalised_kendall_tau_distance(s, hyp2)

            tau,pvalue = statStuff.compute_kendall_tau(s, hyp2)

            ##print('tau:',tau)
            tau = (tau + 1) / 2
            v = v + tau
    #if len(res)!=0:
    #    ratio=v/len(res)
    #else:
    #    ratio=0
    if len(res) >1:
        ratio=v/normalize
    else:
        ratio=0
    return v,ratio,querytime





def get_state_sample(conn, measBase, table, sel, sampleSize, state):

    querySample = "SELECT "+sel+", "+measBase+" FROM "+table+" where "+sel+" = '"+str(state)+"' limit "+str(sampleSize)+";"
    ##print('stat query:', querySample)
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
    ##print(raw_comparisons)
    rejected, corrected = fdrcorrection([x[3] for x in raw_comparisons], alpha=0.05)
    for i in range(len(raw_comparisons)):
        if rejected[i]:
            w_comparisons_rej.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))
        else:
            w_comparisons.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))

    #print("NB de comparaisons significatives (welch)", len(w_comparisons))
    # #print_comp_list(sorted(w_comparisons, key=lambda x: x[0] + x[1]))
    by_prox_to_threshold = sorted(w_comparisons, key=lambda x: abs(0.05 - x[2]), reverse=True)
    # #print(by_prox_to_threshold)

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

    #print("NB de comparaisons significatives (welch + X param)", len(final))
    # #print_comp_list(sorted(final, key=lambda x: x[0] + x[1]))

    # borda hypothesis
    patrick_format = [(a, b, 1, None, None) for (a, b, c) in final]
    ##print('alex pairwise:',patrick_format)
    hypothesis = computeRanksForAll(patrick_format, adom).items()
    ##print(hypothesis)
    hypothesis = sorted(
        hypothesis,
        key=lambda x: x[1],
        reverse=True
    )
    ##print(hypothesis)
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

    ##print('Alex hypothesis:',correctHyp)
    return correctHyp


def materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex):
    # generate and get all materialized cuboids
    #print("Creating views")
    dbStuff.dropAllMVs(conn)
    dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex)
    mvnames = dbStuff.getMVnames(conn)

    aggQueries = dbStuff.getAggQueriesOverMV(mvnames, sel)
    #print("Materializing ", len(mvnames), " views: ", mvnames)
    # #print("queries: ",aggQueries)
    #print("Number of aggregate queries over the MVs: ", len(aggQueries))
    return mvnames,aggQueries

def materializeViewsWithoutIndex(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice):
    # generate and get all materialized cuboids
    #print("Creating views")
    dbStuff.dropAllMVs(conn)
    dbStuff.createMVWithoutIndex(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice)
    mvnames = dbStuff.getMVnames(conn)

    aggQueries = dbStuff.getAggQueriesOverMV(mvnames, sel)
    #print("Materializing ", len(mvnames), " views: ", mvnames)
    # #print("queries: ",aggQueries)
    #print("Number of aggregate queries over the MVs: ", len(aggQueries))
    return mvnames,aggQueries

def getHypothesisAllComparisons(conn, meas, measBase, table, sel, valsToSelect, sampleSize, method='SYSTEM_ROWS'):
        # checking all comparisons
        correctHyp, samplingTime, hypothesisGenerationTime, pvalue = generateHypothesisTestDolap(conn, meas, measBase,
                                                                                                 table, sel, sampleSize,
                                                                                                 method, valsToSelect)
        ##print('all comp. hypothesis:', correctHyp)
        return correctHyp, samplingTime, hypothesisGenerationTime, pvalue

def hypothesisGeneration(conn, prefs, sel, measBase, meas, table, sampleSize, allComparison):
    if allComparison == False:
        # sampling
        start_time = time.time()
        adom, congress = fetchCongressionalSample(conn, sel, table, measBase, sampleSize, adom_restr=prefs)
        end_time = time.time()
        samplingTime = end_time - start_time
        #print('sampling time:', samplingTime)

        # compute hypothesis
        start_time = time.time()
        hypothesis = getHypothesisCongressionalSampling(adom, congress)
        end_time = time.time()
        hypothesisGenerationTime = end_time - start_time
        #print('Hypothesis generation time:', hypothesisGenerationTime)
    else:
        # sampling and hypothesis
        #start_time = time.time()
        hypothesis,samplingTime, hypothesisGenerationTime,pvalue = getHypothesisAllComparisons(conn, meas, measBase, table, sel, tuple(prefs), sampleSize,
                                                 method='SYSTEM_ROWS')
        #end_time = time.time()
        #samplingTime = end_time - start_time
        #hypothesisGenerationTime = samplingTime
        # #print('sampling time:', samplingTime)
        ##print('Hypothesis generation time:', hypothesisGenerationTime)
    return hypothesis,hypothesisGenerationTime,samplingTime, pvalue



def test(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice, groupbyAtt, sel, measBase, meas, function,table,
         sampleSize,comparison,generateIndex,allComparison,ratioOfQuerySample,mvnames,aggQueries,currentSample,cumulate):
    ##print("Sample size: ",sampleSize)

    hypothesis,hypothesisGenerationTime,samplingTime,pvalue=hypothesisGeneration(conn, prefs, sel, measBase, meas, table, sampleSize, allComparison)
    ##print("Hypothesis predicted: ", hypothesis)

    # only ok if hypothesis is a<B or a>b
    if len(hypothesis) <2 or hypothesis[0][1]==hypothesis[1][1]:
        return 99,99,99,99,currentSample #ugly
    else:
        #print("Hypothesis predicted: ", hypothesis)

        limitedHyp = []
        valsToSelect = []
        j = 0
        for h in hypothesis:
            if (h[1] <= nbAdomVals and j < nbAdomVals):
                limitedHyp.append(h)
                valsToSelect.append(h[0])
                j = j + 1

        start_time = time.time()


        sizeofquerysample = int(ratioOfQuerySample * len(aggQueries))
        if sizeofquerysample==0:
            sizeofquerysample=1
        ##print("ratio: ",ratioOfQuerySample)
        ##print("len agg: ",len(aggQueries))
        ##print('Size of query sample:', sizeofquerysample)


        pwrset=aggQueries
        ##print("pwrset:",pwrset)

        ##print("Sampling aggregate queries")
        if cumulate==False:
            ranks, queryCountviolations, queryCountCuboid, cuboid, newpset = bounders.getSample(pwrset, sel, measBase, function,
                                                                                 table, tuple(valsToSelect), limitedHyp,
                                                                                 mvnames,False,False,
                                                                                       sizeofquerysample)
            # queryCountviolations, queryCountCuboid, cuboid=bernstein.getSample(proba, error, pwrset, sel, measBase, function, table, tuple(valsEmptyGB), emptyGBresult, mvnames)
        else:
            if currentSample=={}:
                ranks, queryCountviolations, queryCountCuboid, cuboid, newpset = bounders.getSample(pwrset, sel, measBase, function,
                                                                                           table, tuple(valsToSelect),
                                                                                           limitedHyp,
                                                                                           mvnames, False, False,
                                                                                           sizeofquerysample)
                currentSample["ranks"]=ranks
                currentSample["queryCountviolations"]=queryCountviolations
                currentSample["queryCountCuboid"]=queryCountCuboid
                currentSample["cuboid"]=cuboid
                currentSample["pset"]=newpset
            else:
                ranksTemp, queryCountviolationsTemp, queryCountCuboidTemp, cuboidTemp, newpset  = bounders.getMoreRandomQueries(sizeofquerysample,currentSample,
                                                                                                                                sel, measBase, function,
                                                                                           table, tuple(valsToSelect),
                                                                                           limitedHyp,
                                                                                           mvnames, False, False)
                ##print(ranksTemp)
                currentSample["ranks"].extend(ranksTemp)
                ##print(currentSample["ranks"])
                currentSample["queryCountviolations"].extend(queryCountviolationsTemp)
                currentSample["queryCountCuboid"].extend(queryCountCuboidTemp)
                currentSample["cuboid"].extend(cuboidTemp)
                currentSample["pset"] = newpset
                ranks=currentSample["ranks"]
                queryCountviolations=currentSample["queryCountviolations"]
                queryCountCuboid=currentSample["queryCountCuboid"]
                cuboid= currentSample["cuboid"]


        #print("Validating: computing violations")

        tabRandomVar = []
        nbViewOK = 0

        nbInconclusive=0
        sizeofsample=sizeofquerysample

        totalQueryTime=0

        for i in range(len(ranks)):

            v,ratio,qtime=countViolations(conn,ranks[i],hypothesis)
            totalQueryTime=totalQueryTime+qtime
            c = dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]

            if c!=0:

                if ratio < ratioViolations:
                    tabRandomVar.append(1)
                    nbViewOK = nbViewOK + 1
                else:
                    tabRandomVar.append(0)
            else:
                sizeofsample = sizeofsample - 1
                nbInconclusive=nbInconclusive+1


        end_time = time.time()
        validationTime = end_time - start_time
        #print('Validation time:', validationTime)

        #print('Number of inconclusive: ',nbInconclusive, ' ratio: ',nbInconclusive/len(queryCountviolations))

        variance = np.var(tabRandomVar)
        # #print('variance: ', variance)
        # check if sizeofsample=0!
        if sizeofsample==0:
            prediction = 0
            #print("WARNING: Nothing conclusive")
            bennetError=0
        else:
            prediction = nbViewOK / sizeofsample

            predictionNbOk = prediction * len(pwrset)
            #print('Number of views ok: ', nbViewOK, 'out of ', sizeofsample, 'views, i.e., rate of:', nbViewOK / sizeofsample)
            #print('Prediction is:', predictionNbOk)

            #nbErrors = 2
            ##print('probability of making ', nbErrors, ' errors: ', bernstein.bernsteinBound(variance, nbErrors))
            ##print('the error (according to Bernstein) for sum and confidence interval of size', proba, ' is: ',
            #      bernstein.bersteinError(proba, variance))
            bennetError=bounders.bennetErrorOnAvg(proba, variance, sizeofsample)
            #print('The error (according to Bennet) for confidence interval of size', proba, ' is: ', bounders.bennetErrorOnAvg(proba, variance, sizeofsample))



        if comparison==True:
            # comparison with ground truth
            #print('*** comparison to ground truth ***')

            #dbStuff.dropAllMVs(conn)
            #nbMVs = dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, 1,generateIndex)
            nbMVs=len(pwrset)

            #compareHypToGB(hypothesis, conn, measBase,function, sel,tuple(valsToSelect),mvnames,table)

            ranks, queryCountviolations, queryCountCuboid, cuboid = bounders.generateAllqueriesOnMVs(pwrset, sel, measBase, function,
                                                                                          table, tuple(valsToSelect),
                                                                                          limitedHyp, mvnames)


            nbInconclusive=0
            tabRandomVar = []
            nbViewOK = 0
            for i in range(len(queryCountviolations)):
                ##print('gt violations:',queryCountviolations[i])
                ##print('gt count:',queryCountCuboid[i])
                #v = dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
                v,ratio,qtime = countViolations(conn, ranks[i], hypothesis)
                c = dbStuff.execute_query(conn, queryCountCuboid[i])[0][0]
                # #print(v)
                # #print(c)
                if c != 0:
                    #OLD #print(v / c, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)

                    ##print(ratio, " violation rate in cuboid ", cuboid[i], " of size: ", c, ". Number of violations: ", v)
                    if ratio < ratioViolations:
                        tabRandomVar.append(1)
                        nbViewOK = nbViewOK + 1
                    else:
                        tabRandomVar.append(0)
                else:
                    ##print("inconclusive, not enough tuples in cuboid for select values")
                    nbMVs=nbMVs-1
                    nbInconclusive=nbInconclusive+1

            #variance = np.var(tabRandomVar)
            # #print('variance: ', variance)

            #print('number of inconclusive: ', nbInconclusive, ' ratio: ', nbInconclusive / len(queryCountviolations))

            #print('nb of views ok: ', nbViewOK, 'out of ', nbMVs, 'views, i.e., rate of:', nbViewOK / nbMVs)
            gtratio= nbViewOK / nbMVs

            realError=abs(prediction - (nbViewOK / nbMVs))
            #print('Error is: ', abs(prediction - (nbViewOK / nbMVs)))


            return prediction,bennetError,realError,gtratio,currentSample
        else:
            #return totalQueryTime, samplingTime, hypothesisGenerationTime, validationTime

            #if no comparison only outputs bennet error and prediction
            return prediction, bennetError, bennetError, prediction, currentSample




def groundTruthForQueriesOverMVs(pairs,groupbyAtt,pred=-1,error=1):
    dict = {}
    #ratioOfQuerySample=1
    #initsampleSize=1
    currentSample={}
    #sampleSize = initsampleSize * sizeOfR
    for p in pairs:

        predT, bennetError, minErrorT, gtratio,currentSample = test(conn, nbAdomVals, p,
                                                           ratioViolations, proba, error,
                                                           percentOfLattice, groupbyAtt,
                                                           sel, measBase, meas, function, table,
                                                           sizeOfR, comparison,
                                                           generateIndex, allComparisons,
                                                           1, mvnames,
                                                           aggQueries, currentSample,
                                                           cumulate=False)

        if minErrorT != 99:
            ##print(meanError)
            #e = 0
            #minErrorT = meanError
            #predT = meanPred

            if True:
            #if predT >=pred and minErrorT < error:
                dict[p] = [minErrorT, predT]

    dict = utilities.sort_dict_by_second_entry_desc(dict)
    return dict


# returns the pairs found on all the lattice
def groundTruthAllLatice(pairs,groupbyAtt,pred=-1,error=1):
    dict = {}
    #ratioOfQuerySample=1
    #initsampleSize=1
    #percentOfLattice=1
    currentSample = {}
    #sampleSize = initsampleSize * sizeOfR

    mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, 1,
                                           generateIndex)

    for p in pairs:

        predT, bennetError, minErrorT, gtratio,currentSample = test(conn, nbAdomVals, p,
                                                        ratioViolations, proba, error,
                                                        1, groupbyAtt,
                                                        sel, measBase, meas, function, table,
                                                        sizeOfR, comparison,
                                                        generateIndex, allComparisons,
                                                        1, mvnames,
                                                        aggQueries, currentSample,
                                                        cumulate=False)


        if minErrorT != 99:
            ##print(meanError)
            #e = 0
            #minErrorT = meanError
            #predT = meanPred

            if predT >=pred and minErrorT < error:
                dict[p] = [minErrorT, predT]


    dict = utilities.sort_dict_by_second_entry_desc(dict)
    dbStuff.dropAllMVs(conn)
    return dict




def plotRuns(dictRuns, tabTest, nbruns, x_label='Size of query sample', y_label='F-measure',
                                          title='F-measure by sample size'):
    dataRuns = []
    for t in tabTest:
        if t in dictRuns:
            x = dictRuns[t]
            tabmean = []
            tabstdev = []
            for y in range(len(tabTest)):
                tabtemp = []
                for r in range(nbruns):
                    tabtemp.append(x[r][y])
                # #print(tabtemp)
                tabmean.append(statistics.mean(tabtemp))
                tabstdev.append(statistics.stdev(tabtemp))
            # #print(tabmean)
            dataRuns.append({'x': tabTest, 'y': tabmean, 'yerr': tabstdev, 'label': t})
    plotStuff.plot_curves_with_error_bars(dataRuns, x_label, y_label,
                                          title)



def comparisonToGT(groupbyAtt,allLattice):
    sel = groupbyAtt[0]

    groupbyAtt = groupbyAtt[1:]

    nbpairs = 90

    pairs = dbStuff.generateAllPairs(conn, sel, table, nbpairs)


    #tabTest = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
    tabTest = (0.6, 1)

    if allLattice:
        dictGT = groundTruthAllLatice(pairs,groupbyAtt)
        mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,
                                           generateIndex)
    else:
        mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,
                                               generateIndex)
        dictGT = groundTruthForQueriesOverMVs(pairs,groupbyAtt, -1,1)


    column_names = ['Runs', 'Initial Sample', 'Query Sample', 'Pair', 'Error', 'F1', 'Recall', 'Precision', 'Recall@10']

    # Create an empty DataFrame with the specified columns
    df = pd.DataFrame(columns=column_names)

    for nr in tqdm(range(nbruns)):

        for initsampleSize in tabTest:
            #print("INIT SAMPLE SIZE: ", initsampleSize)

            sampleSize = initsampleSize * sizeOfR

            for p in pairs:

                dict={}
                currentSample = {}

                for ratioOfQuerySample in tabTest:

                    predT, bennetError, minErrorT, gtratio, currentSample = test(conn, nbAdomVals, p,
                                                                                 ratioViolations, proba, error,
                                                                                 percentOfLattice, groupbyAtt,
                                                                                 sel, measBase, meas, function, table,
                                                                                 sampleSize, comparison,
                                                                                 generateIndex, allComparisons,
                                                                                 ratioOfQuerySample, mvnames,
                                                                                 aggQueries, currentSample,
                                                                                 cumulate=True)
                    if minErrorT != 99:
                            dict[p] = [minErrorT, predT]

                dict = utilities.sort_dict_by_second_entry_desc(dict)

                p, r, f = utilities.f_measure_first_k_keys(dict, dictGT, 0)

                p10, r10, f10 = utilities.f_measure_first_k_keys(dict, dictGT, 30)

                df.loc[len(df)] = [nr, initsampleSize, ratioOfQuerySample, p, minErrorT, f, r, p, r10]

    df.to_csv(fileResults)


#
# Main

if __name__ == "__main__":

    config = configparser.ConfigParser()

    current_time=time.localtime()
    formatted_time = time.strftime("%d-%m-%y:%H:%M:%S", current_time)
    fileResults='results/res_'+formatted_time+'.csv'

    # The DB we want
    config.read('configs/flightsDolap.ini')
    #config.read('configs/flightsquarterDolap.ini')
    #config.read('configs/ssbDolap.ini')
    #config.read('configs/flights1923Dolap.ini')
    #config.read('configs/flights1923.ini')
    #config.read('configs/artificial.ini')
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
    ##print(tuple(prefs))
    if len(prefs) == 0:
        prefs = None


    ###
    ### PARAMETERS
    ###

    # for query sample size according to Hoeffding
    proba = 0.1 #probability of making an error
    error = 0.1 #error

    # number of values of adom to consider - default = all of prefs
    nbAdomVals = len(prefs)

    # for sampling fact table with Postgresql
    initsampleSize = 0.4
    samplingMethod = 'SYSTEM_ROWS'  # or SYSTEM

    # ratio max of violations in a cuboid
    ratioViolations = 0.4

    # ratio min of cuboids with ratio violations < ratioViolations
    ratioCuboidOK = 0.8

    # percentage of the lattice to generate
    percentOfLattice = 0.4

    # do we generate indexes?
    # possible values:
    # True (create index on sel attribute),
    # False (no index),
    # mc (one multicolumn index, sel first), so far only over views (not fact table)
    generateIndex = 'mc'
    #generateIndex = False

    # do we compare to ground truth? Otherwise, efficiency is tested
    comparison = False

    # do we generate all comparisons?
    allComparisons = True

    # ratio of sample of queries for validation
    ratioOfQuerySample = 0.4

    # number of runs
    nbruns=2

    ###
    ### END OF PARAMETERS
    ###


    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    # get size of fact table
    sizeOfR=dbStuff.getSizeOf(conn,table)
    ##print(sizeOfR)

    # to always have the same order in group bys, with sel attribute last
    groupbyAtt.sort()

    nbWrongRanking=0
    resultRuns=[]


    if comparison == True:
        #comparisonToGT(groupbyAtt, True)


        sel = groupbyAtt[0]
        groupbyAtt = groupbyAtt[1:]
        # #print(groupbyAtt)

        # comparison = false if we don't want empirical error
        # comparison = True if we want both empirical and Bennet error
        comparison = True
        nbpairs = 90
        paramTested = list(range(nbpairs))

        pairs = dbStuff.generateAllPairs(conn, sel, table, nbpairs)
        #dict = {}

        sampleSize = 1
        minError = 0.1  # threshold
        pred = 0.4
        #maxPred = 0


        #ratioOfQuerySample = 0.5
        tabTest = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
        #tabTest = (0.6, 1)



        # do we want GT for error/pred on queries over MVs or all lattice?

        dictGT = groundTruthAllLatice()

        mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,
                                               generateIndex)

        #dictGT = groundTruthForQueriesOverMVs( -1,1)

        dictRuns={}
        dictRunsErr={}

        column_names = ['Runs', 'Initial Sample', 'Query Sample', 'Pair','Error','F1', 'Recall', 'Precision', 'Recall@10']

        # Create an empty DataFrame with the specified columns
        df = pd.DataFrame(columns=column_names)

        for nr in tqdm(range(nbruns)):
            ##print("RUN: ",nr)

            data = []
            dataErrorsAllPairs=[]

            for initsampleSize in tabTest:

                #for percentOfLattice in tabTest:
                # comment if not percent of lattice tested
                #mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table,percentOfLattice, generateIndex)

                dataPairs = []
                dataError=[]
                dataStdevError=[]

                sampleSize = initsampleSize * sizeOfR
                #currentSample = {}

                for ratioOfQuerySample in tabTest:
                    print("INIT SAMPLE : ", initsampleSize," QUERY SAMPLE : ", ratioOfQuerySample)
                    dict = {}
                    # data = []
                    #timings = []
                    tabError=[]


                    for p in pairs:

                        currentSample = {}

                        predT, bennetError, minErrorT, gtratio,currentSample = test(conn, nbAdomVals, p,
                                                                                 ratioViolations, proba, error,
                                                                                 percentOfLattice, groupbyAtt,
                                                                                 sel, measBase, meas, function, table,
                                                                                 sampleSize, comparison,
                                                                                 generateIndex, allComparisons,
                                                                                 ratioOfQuerySample, mvnames,
                                                                                 aggQueries, currentSample,
                                                                                 cumulate=True)


                        ##print("output of Test: ", p, meanError, meanPred, meanBennet)

                        if minErrorT != 99:
                            tabError.append(minErrorT)

                            dict[p] = [minErrorT, predT]


                    dict = utilities.sort_dict_by_second_entry_desc(dict)

                    #scoreComp = utilities.jaccard_score_first_k_keys(dict, dictGT, 0)
                    precision, recall, f1 = utilities.f_measure_first_k_keys(dict, dictGT, 0)
                    scoreComp = f1

                    p10, r10, f10 = utilities.f_measure_first_k_keys(dict, dictGT, 10)
                    scoreComp = r10

                    # if we want the number of pairs
                    # dataPairs.append(len(dict))

                    # if we want the comparison with GT
                    dataPairs.append(scoreComp)

                    avgError=statistics.mean(tabError)
                    stdevError=statistics.stdev(tabError)
                    dataError.append(avgError)
                    dataStdevError.append((stdevError))

                    #['Runs', 'Initial Sample', 'Query Sample', 'Pair','Error','F1', 'Recall', 'Precision', 'Recall@10']
                    df.loc[len(df)]=[nr,initsampleSize,ratioOfQuerySample,p,minErrorT,f1,recall,precision,r10]



                # plots number of pairs with error<0.1 by size of query sample
                stdevPairs = [0] * len(tabTest)

                # change percentOfLattice by initsampleSize when changing external for loop
                data.append({'x': tabTest, 'y': dataPairs, 'yerr': stdevPairs, 'label': initsampleSize})
                dataErrorsAllPairs.append({'x': tabTest, 'y': dataError, 'yerr': dataStdevError, 'label': initsampleSize})

                if initsampleSize in dictRuns:
                    dictRuns[initsampleSize].append(dataPairs)
                    dictRunsErr[initsampleSize].append(dataError)
                else:
                    dictRuns[initsampleSize]=[dataPairs]
                    dictRunsErr[initsampleSize] = [dataError]


        plotRuns(dictRuns, tabTest, nbruns, x_label='Size of query sample', y_label='Recall@10',
                 title='Recall@10 by sample size')
        plotRuns(dictRunsErr, tabTest, nbruns, x_label='Size of query sample', y_label='Error',
                 title='Error by sample size')

        df.to_csv(fileResults)


    else:

        sel = groupbyAtt[0]
        groupbyAtt = groupbyAtt[1:]
        #sel ='departure_airport'
        #groupbyAtt = ['airline', 'date', 'departure_hour', 'flight']

        #41616 for airport, 91 for airline
        nbpairs = 90

        paramTested = list(range(nbpairs))
        pairs = dbStuff.generateAllPairs(conn, sel, table, nbpairs)

        # do we generate indexes?
        # possible values:
        # True (create index on sel attribute),
        # False (no index),
        # mc (one multicolumn index, sel first), so far only over views (not fact table)
        generateIndex = 'mc'
        # generateIndex = False

        mvnames, aggQueries = materializeViewsWithoutIndex(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice)

        data=[]
        for generateIndex in [False, True, 'mc']:
        #for generateIndex in ['mc']:
        #for percentOfLattice in [0.4,0.5,0.6,0.7]:

            #mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,generateIndex)
            dbStuff.dropAllIndexOnMVs(conn,mvnames)
            dbStuff.generateIndexesOnMVs(conn, sel, mvnames, generateIndex)

            #ratioOfQuerySample = 0.5
            #tabTest = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)


            timings=[]
            currentSample = {}
            for p in pairs:

                start_time = time.time()


                sampleSize = initsampleSize * sizeOfR

                prediction, bennetError, realError, gtratio,currentSample = test(conn, nbAdomVals, p,
                                                                   ratioViolations, proba, error,
                                                                   percentOfLattice, groupbyAtt,
                                                                   sel, measBase, meas, function, table,
                                                                   sampleSize, comparison,
                                                                   generateIndex, allComparisons,
                                                                   ratioOfQuerySample, mvnames,
                                                                   aggQueries, currentSample,
                                                                   cumulate=False)

                end_time = time.time()
                timings.append(end_time - start_time)

            # TIMINGS
            timings = utilities.accumulate_numbers(timings)
            # #print(timings)
            stdevTiming = [0] * nbpairs
            #print()
            data.append(
                {'x': paramTested, 'y': timings, 'yerr': stdevTiming, 'label': generateIndex}
            )


        plotStuff.plot_curves_with_error_bars(data, x_label='Number of pairs', y_label='Time (s)',title='Times')


    # Clean and cose the connection
    dbStuff.dropAllMVs(conn)
    close_connection(conn)



