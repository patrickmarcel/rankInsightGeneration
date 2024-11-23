import math
import statistics
import numpy as np
import time
import configparser
import json
from statsmodels.stats.multitest import fdrcorrection

import dbStuff
import plotStuff
import statStuff
import utilities
from plotStuff import plot_curves_with_error_bars
from dbStuff import execute_query, connect_to_db, close_connection, getSample
from statStuff import permutation_test, compute_skewness, claireStat
from rankingFromPairwise import computeRanksForAll, generateHypothesisTest
import bounders
import tests


# ------  Debug ?  ------------
DEBUG_FLAG = True


def  compareHypToGB(hypothesis, conn, measBase,function, sel, vals,mvnames, table):
    #query="Select " + sel +" from  " + sel + " where " + sel + " in " +  str(vals) + " order by " + measBase + " desc;"
    materialized=bounders.findMV(mvnames,sel,table)
    #print("MATERIALIZED: ",materialized)
    query="select 'all',string_agg(" + sel + ",',') from (SELECT " + sel + ", " + function + "(" + measBase + "),  rank () over (  order by " + function + "(" + measBase + ") desc ) as rank FROM \"" + materialized + "\" WHERE " + sel + " in " + str(vals) +" group by " +  sel + " order by rank);"
    #print(query)
    v,ratio,qtime=countViolations(conn,query,hypothesis)
    #print(v)
    print("hypothesis compared to group by ",sel," has ",v," violations")
    return v



def countViolations(conn,query,hypothesis):
    #print(query)
    hyp=[str(a) for (a,b) in hypothesis]
    #print('hyp:',hyp)
    v=0
    start_time_q = time.time()
    res=dbStuff.execute_query(conn,query)
    end_time_q = time.time()
    querytime=end_time_q-start_time_q
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
            #tau=statStuff.normalised_kendall_tau_distance(s,hyp)

            tau,pvalue=statStuff.compute_kendall_tau(s,hyp)

            #print('tau:',tau)
            tau=(tau+1)/2
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
            #tau = statStuff.normalised_kendall_tau_distance(s, hyp2)

            tau,pvalue = statStuff.compute_kendall_tau(s, hyp2)

            #print('tau:',tau)
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



def getHypothesisAllComparisons(conn, meas, measBase, table, sel,valsToSelect, sampleSize, method='SYSTEM_ROWS'):
    # checking all comparisons
    correctHyp,samplingTime, hypothesisGenerationTime = generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, method, valsToSelect)
    #print('all comp. hypothesis:', correctHyp)
    return correctHyp, samplingTime, hypothesisGenerationTime

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


def materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex):
    # generate and get all materialized cuboids
    print("Creating views")
    dbStuff.dropAllMVs(conn)
    dbStuff.createMV(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex)
    mvnames = dbStuff.getMVnames(conn)

    aggQueries = dbStuff.getAggQueriesOverMV(mvnames, sel)
    print("Materializing ", len(mvnames), " views: ", mvnames)
    # print("queries: ",aggQueries)
    print("Number of aggregate queries over the MVs: ", len(aggQueries))
    return mvnames,aggQueries

def hypothesisGeneration(conn, prefs, sel, measBase, meas, table, sampleSize, allComparison):
    if allComparison == False:
        # sampling
        start_time = time.time()
        adom, congress = fetchCongressionalSample(conn, sel, table, measBase, sampleSize, adom_restr=prefs)
        end_time = time.time()
        samplingTime = end_time - start_time
        print('sampling time:', samplingTime)

        # compute hypothesis
        start_time = time.time()
        hypothesis = getHypothesisCongressionalSampling(adom, congress)
        end_time = time.time()
        hypothesisGenerationTime = end_time - start_time
        print('Hypothesis generation time:', hypothesisGenerationTime)
    else:
        # sampling and hypothesis
        #start_time = time.time()
        hypothesis,samplingTime, hypothesisGenerationTime = getHypothesisAllComparisons(conn, meas, measBase, table, sel, tuple(prefs), sampleSize,
                                                 method='SYSTEM_ROWS')
        #end_time = time.time()
        #samplingTime = end_time - start_time
        #hypothesisGenerationTime = samplingTime
        # print('sampling time:', samplingTime)
        #print('Hypothesis generation time:', hypothesisGenerationTime)
    return hypothesis,hypothesisGenerationTime,samplingTime


def test(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice, groupbyAtt, sel, measBase, meas, function,table,
         sampleSize,comparison,generateIndex,allComparison,ratioOfQuerySample,mvnames,aggQueries,currentSample,cumulate):
    #print("Sample size: ",sampleSize)

    hypothesis,hypothesisGenerationTime,samplingTime=hypothesisGeneration(conn, prefs, sel, measBase, meas, table, sampleSize, allComparison)
    print("Hypothesis predicted: ", hypothesis)

    # only ok if hypothesis is a<B or a>b
    if len(hypothesis) <2 or hypothesis[0][1]==hypothesis[1][1]:
        return 99,99,99,99 #ugly
    else:

        limitedHyp = []
        valsToSelect = []
        j = 0
        for h in hypothesis:
            if (h[1] <= nbAdomVals and j < nbAdomVals):
                limitedHyp.append(h)
                valsToSelect.append(h[0])
                j = j + 1
        #print("Hypothesis limited to choosen values: ", limitedHyp)
        #print("vals: ",valsToSelect)

        # just for checking on groupBy sel
        # emptyGBresult, emptyGBresultAll = emptyGB(conn, nbAdomVals, table, sel, meas)
        # print("Empty GB says:", emptyGBresult)
        # valsEmptyGB=[a for (a, b) in emptyGBresult]
        # print(valsEmptyGB)

        """
        #generate index on sel attribute over table R
        if generateIndex == True :
            print('Creating indexes')
            dbStuff.dropAllIndex(conn,table)
            dbStuff.generateHashIndex(conn, table, sel)
        else:
            if generateIndex == 'mc':
                print('Creating indexes')
                dbStuff.dropAllIndex(conn, table)
                dbStuff.generateHashIndex(conn, table, sel)
                gat=''
                for g in groupbyAtt:
                    gat=gat+g+','
                gat=gat+sel
                dbStuff.generateMulticolIndex(conn, table, gat, sel)
            else: #false
                dbStuff.dropAllIndex(conn, table)
                #dbStuff.dropIndex(conn, table, sel)
        """

        #MVs creation was here


        #validation of hypothesis
        start_time = time.time()

        #size of query sample
        #sizeofsample = int(bounders.sizeOfSampleHoeffding(proba, error)) + 1
        #print('size of query sample according to Hoeffding:', sizeofsample)
        #sizeofsample = sizeofquerysample

        sizeofquerysample = int(ratioOfQuerySample * len(aggQueries))
        if sizeofquerysample==0:
            sizeofquerysample=1
        #print("ratio: ",ratioOfQuerySample)
        #print("len agg: ",len(aggQueries))
        print('Size of query sample:', sizeofquerysample)


        # total number of cuboids
        #N = len(utilities.powerset(groupbyAtt))
        #print('size of sample according to Bardenet:',
        #      int(bounders.sizeOfSampleHoeffdingSerflingFromBardenet(proba, error, N)) + 1)

        #pwrset = dbStuff.getCuboidsOfAtt(groupbyAtt, sel)
        pwrset=aggQueries
        #print("pwrset:",pwrset)

        print("Sampling aggregate queries")
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
                ranksTemp, queryCountviolationsTemp, queryCountCuboidTemp, cuboidTemp, newpset  = bounders.getMoreRandamQueries(sizeofquerysample,currentSample,
                                                                                                                                sel, measBase, function,
                                                                                           table, tuple(valsToSelect),
                                                                                           limitedHyp,
                                                                                           mvnames, False, False)
                #print(ranksTemp)
                currentSample["ranks"].extend(ranksTemp)
                #print(currentSample["ranks"])
                currentSample["queryCountviolations"].extend(queryCountviolationsTemp)
                currentSample["queryCountCuboid"].extend(queryCountCuboidTemp)
                currentSample["cuboid"].extend(cuboidTemp)
                currentSample["pset"] = newpset
                ranks=currentSample["ranks"]
                queryCountviolations=currentSample["queryCountviolations"]
                queryCountCuboid=currentSample["queryCountCuboid"]
                cuboid= currentSample["cuboid"]


        print("Validating: computing violations")

        tabRandomVar = []
        nbViewOK = 0

        nbInconclusive=0
        sizeofsample=sizeofquerysample

        totalQueryTime=0

        for i in range(len(ranks)):
            #print(len(ranks))
            #print(i,ranks[i])
            v,ratio,qtime=countViolations(conn,ranks[i],hypothesis)
            totalQueryTime=totalQueryTime+qtime
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
        print('Validation time:', validationTime)

        print('Number of inconclusive: ',nbInconclusive, ' ratio: ',nbInconclusive/len(queryCountviolations))

        variance = np.var(tabRandomVar)
        # print('variance: ', variance)
        # check if sizeofsample=0!
        if sizeofsample==0:
            prediction = 0
            print("WARNING: Nothing conclusive")
            bennetError=0
        else:
            prediction = nbViewOK / sizeofsample

            predictionNbOk = prediction * len(pwrset)
            print('Number of views ok: ', nbViewOK, 'out of ', sizeofsample, 'views, i.e., rate of:', nbViewOK / sizeofsample)
            print('Prediction is:', predictionNbOk)

            #nbErrors = 2
            #print('probability of making ', nbErrors, ' errors: ', bernstein.bernsteinBound(variance, nbErrors))
            #print('the error (according to Bernstein) for sum and confidence interval of size', proba, ' is: ',
            #      bernstein.bersteinError(proba, variance))
            bennetError=bounders.bennetErrorOnAvg(proba, variance, sizeofsample)
            print('The error (according to Bennet) for confidence interval of size', proba, ' is: ',
                  bounders.bennetErrorOnAvg(proba, variance, sizeofsample))



            #if sizeofsample>1:
            #    print('The error (empirical Bennet from Maurer and Pontil) for confidence interval of size', proba, ' is: ',
            #      bounders.empiricalBennetFromMaurer(proba, variance, sizeofsample))

            ###
            ### IF REPORTING EMPIRICAL ERROR
            ### UNCOMMENT BELOW
            #bennetError = bernstein.empiricalBennetFromMaurer(proba, variance, sizeofsample)

            #print('the error (according to bardenet) for avg and confidence interval of size', proba, ' is: ',
            #      bernstein.empiricalBernsteinFromBardenet(proba, variance, sizeofsample, N))



        if comparison==True:
            # comparison with ground truth
            print('*** comparison to ground truth ***')

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
                #print('gt violations:',queryCountviolations[i])
                #print('gt count:',queryCountCuboid[i])
                #v = dbStuff.execute_query(conn, queryCountviolations[i])[0][0]
                v,ratio,qtime = countViolations(conn, ranks[i], hypothesis)
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
            #return totalQueryTime, samplingTime, hypothesisGenerationTime, validationTime

            #if no comparison only outputs bennet error and prediction
            return prediction, bennetError, bennetError, prediction


def groundTruth(minError):
    dict = {}
    ratioOfQuerySample=1
    initsampleSize=1
    for p in pairs:

        meanError, meanPred, meanBennet = tests.testAccuracyQuerySampleSizeDOLAP(tabTest, mvnames, aggQueries, nbruns,
                                                                                 conn,
                                                                                 nbAdomVals, p, ratioViolations, proba,
                                                                                 error, percentOfLattice,
                                                                                 groupbyAtt, sel,
                                                                                 measBase, meas, function, table,
                                                                                 comparison, generateIndex,
                                                                                 allComparisons, initsampleSize,
                                                                                 sizeOfR, ratioCuboidOK,
                                                                                 ratioOfQuerySample, cumulate=True)


        if meanError != []:
            print(meanError)
            e = 0
            while meanError[e] >= minError and e < len(meanError) - 1:
                # print(e)
                e = e + 1
            sampleSizeT = e / 10
            sampleSizeT = ratioOfQuerySample
            minErrorT = meanError[e]
            predT = meanPred[e]
            # if minErrorT < minError and predT > pred and minErrorT < 0.1 and sampleSizeT > 0:
            if minErrorT < minError:
                dict[p] = [minErrorT, predT]


    dict = utilities.sort_dict_by_second_entry_desc(dict)
    return dict

#
# Main

if __name__ == "__main__":

    config = configparser.ConfigParser()

    # The DB wee want
    #config.read('configs/flights1923.ini')
    #config.read('configs/flightsquarterDolap.ini')
    config.read('configs/flightsDolap.ini')
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
    #print(tuple(prefs))
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
    initsampleSize = 0.5
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
    comparison = True

    # do we generate all comparisons?
    allComparisons = True

    # ratio of sample of queries for validation
    ratioOfQuerySample = 0.4

    # number of runs
    #nbOfRuns = 1 ## no more used
    nbruns=1

    ###
    ### END OF PARAMETERS
    ###


    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    # get size of fact table
    sizeOfR=dbStuff.getSizeOf(conn,table)
    #print(sizeOfR)

    # to always have the same order in group bys, with sel attribute last
    groupbyAtt.sort()

    nbWrongRanking=0
    resultRuns=[]
    timings=[]

    if comparison==True:

        # todo for on measures
        # todo for testedAtt in groupbyAtt:
        # todo compare to ground truth

        # size of query sample according to Hoeffding
        sizeHoeffding = int(bounders.sizeOfSampleHoeffding(proba, error))
        print('size of query sample according to Hoeffding for proba=',proba," and error=",error,": ",sizeHoeffding)



        sel=groupbyAtt[0]
        groupbyAtt=groupbyAtt[1:]
        #print(groupbyAtt)




        # comparison = false if we don't want empirical error
        # comparison = True if we want both empirical and Bennet error
        comparison = True
        nbpairs=90
        paramTested=list(range(nbpairs))

        pairs=dbStuff.generateAllPairs(conn, sel, table,nbpairs)
        dict={}

        sampleSize = 1
        minError = 0.2 #threshold
        pred = 0
        maxPred = 0

        mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,
                                               generateIndex)

        # total number of cuboids
        N = len(aggQueries)
        print('size of sample according to Bardenet:',
              int(bounders.sizeOfSampleHoeffdingSerflingFromBardenet(proba, error, N)))

        ratioOfQuerySample = 0.5
        tabTest = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
        #tabTest = (0.8, 0.9, 1)

        dictGT=groundTruth(minError)

        data = []
        for initsampleSize in tabTest:
            dataPairs = []
            for ratioOfQuerySample in tabTest:
                dict = {}
                #data = []
                timings=[]
                for p in pairs:
                    start_time = time.time()

                    meanError, meanPred, meanBennet=tests.testAccuracyQuerySampleSizeDOLAP(tabTest,mvnames, aggQueries,nbruns, conn,
                                                              nbAdomVals, p, ratioViolations, proba, error, percentOfLattice,
                                                              groupbyAtt, sel,
                                measBase, meas, function, table, comparison, generateIndex,
                                allComparisons, initsampleSize, sizeOfR, ratioCuboidOK, ratioOfQuerySample, cumulate=True)

                    print(p,meanError,meanPred)
                    if meanError!=[]:
                        print(meanError)
                        e=0
                        while meanError[e] >=minError and e<len(meanError)-1:
                            #print(e)
                            e=e+1
                        sampleSizeT=e/10
                        sampleSizeT=ratioOfQuerySample
                        minErrorT=meanError[e]
                        predT=meanPred[e]

                        #if minErrorT < minError and predT > pred and minErrorT < 0.1 and sampleSizeT > 0:
                        if minErrorT < minError:
                            #minError = minErrorT
                            #pred = predT
                            dict[p] = [minErrorT, predT]

                        #print("TEST",predT, maxPred, sampleSizeT, minErrorT)
                        #if predT>maxPred and sampleSizeT>0 and minErrorT<0.1 :
                        #    dict["best"]=[p,sampleSizeT,minErrorT,predT]
                        #    maxPred=predT

                    end_time = time.time()
                    timings.append(end_time - start_time)
                dict=utilities.sort_dict_by_second_entry_desc(dict)
                print("Best: ", dict)
                print("Number of pairs with error < 0.1 (size of dict):",len(dict))

                scoreComp=utilities.jaccard_score_first_k_keys(dict,dictGT,0)
                p,r,f=utilities.f_measure_first_k_keys(dict,dictGT,0)
                scoreComp = f

                # if we want the number of pairs
                #dataPairs.append(len(dict))

                # if we want the comparison with GT
                dataPairs.append(scoreComp)


                #TIMINGS
                timings=utilities.accumulate_numbers(timings)
                #print(timings)
                stdevTiming=[0]*nbpairs
                #data = [
                #    {'x': paramTested, 'y': timings, 'yerr': stdevTiming, 'label': 'Number of pairs'}
                #]

                # uncomment if plot timings
                #plotStuff.plot_curves_with_error_bars(data, x_label='Number of pairs', y_label='Time (s)',title='Times')

            #plots number of pairs with error<0.1 by size of query sample
            stdevPairs = [0] * len(tabTest)
            #print(dataPairs)
            #print(tabTest)

            #data = [
            #    {'x': tabTest, 'y': dataPairs, 'yerr': stdevPairs, 'label': 'Sample size'}
            #]
            data.append({'x': tabTest, 'y': dataPairs, 'yerr': stdevPairs, 'label': initsampleSize})

        plotStuff.plot_curves_with_error_bars(data, x_label='Size of query sample', y_label='F-measure', title='Pairs by sample')






        #tests.testAccuracyInitSampleSize(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice,
        #                           groupbyAtt, sel, measBase, meas, function, table, comparison, generateIndex,
        #                           allComparisons,
        #                           initsampleSize, sizeOfR, nbOfRuns, ratioCuboidOK,
        #                           ratioOfQuerySample, cumulate=False)
    #else:
    #    tests.testTimingsQuerySampleSize(nbruns,conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice, groupbyAtt, sel,
    #                measBase, meas, function,table, comparison, generateIndex,
    #               allComparisons, initsampleSize, sizeOfR, ratioOfQuerySample, cumulate=True)
        #tests.testTimingsLattice(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice,
        #                   groupbyAtt, sel, measBase, meas, function, table, comparison, generateIndex,
        #                   allComparisons,
        #                   initsampleSize, sizeOfR, nbOfRuns,
        #                   ratioOfQuerySample, cumulate=False)



    # Clean and cose the connection
    dbStuff.dropAllMVs(conn)
    close_connection(conn)



