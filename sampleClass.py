from Config import Config
from Lattice import Lattice

import random
import time
import pandas as pd
from tqdm import tqdm, trange

import statStuff
import utilities
import duckdb
import sqlite3


from dolap import countViolations
from dbStuff import dropAllMVs, getAggQueriesOverMV, createMV, getMVnames, connect_to_db, generateAllPairs, \
    execute_query, getSizeOf, dropAllIndexOnMVs, generateIndexesOnMVs
from bounders import generateAllqueriesOnMVs, findMV
from  Hypothesis import Hypothesis

class Sample:

    def __init__(self, conn,attInGB, selAtt, meas, measBase, function, table, generateIndex=False):
        # create all MVs and aggregates for ground truth
        self.conn=conn
        self.sel=selAtt
        self.meas=meas
        self.measBase=measBase
        self.function=function
        self.table=table
        dropAllMVs(conn)
        #create all MVs
        print("Creating materialized views")
        createMV(conn, attInGB, selAtt, measBase, function, table, 1, generateIndex)
        self.allMVs=getMVnames(conn)
        #create all aggregate queries
        self.aggOverMV = getAggQueriesOverMV(self.allMVs,sel)
        # currentSample is the set of cuboid names in the sample
        self.currentSample=[]
        self.remaining=[]
        self.aggOverMC=[]
        self.MC=[]


    def getCurrentSample(self):
        return self.currentSample

    def generateSampleOfAggQueries(self, ratio):
        sizeOfSample = int(ratio*len(self.aggOverMC))
        if sizeOfSample==0:
            sizeOfSample=1
        remaining = self.aggOverMC.copy()
        #pwrset=self.allMVs
        # if withBias:
        #     tabNb = []
        #     tabWeights = []
        #     sumLength = 0
        #     for p in pwrset:
        #         sumLength = sumLength + len(p)
        #     maxlength = len(pwrset[-1]) + 1
        #     i = 0
        #     for p in pwrset:
        #         w = (maxlength - len(p)) / sumLength
        #         tabWeights.append(w)
        #         tabNb.append(i)
        #         i = i + 1
        chosen = []
        for j in range(sizeOfSample):
            #if withBias:
            #    nb = max(random.choices(tabNb, tabWeights))
            #else:
            nb = random.randint(0, len(remaining)-1)
            if nb==0:
                nb=1
            gb = remaining[nb]
            remaining.remove(gb)
            chosen.append(gb)
        self.currentSample = chosen
        self.remaining = remaining


    def increaseSample(self,newRatio):
        currentSize=len(self.currentSample)
        desriredSize=newRatio*len(self.aggOverMC)
        diff=int(desriredSize-currentSize)
        if len(self.currentSample)==0:
            remaining = self.aggOverMC.copy()
        else:
            remaining = self.remaining.copy()
        chosen = []
        for j in range(diff):
            nb = random.randint(0, len(remaining)-1)
            gb = remaining[nb]
            remaining.remove(gb)
            chosen.append(gb)
        self.currentSample.extend(chosen)
        self.remaining = remaining


    def generateRandomMC(self,percentOfLattice):
        sizeOfSample=int(percentOfLattice*len(self.allMVs))
        remaining=self.allMVs.copy()
        chosen=[]
        for j in range(sizeOfSample):
            nb = random.randint(0, len(remaining)-1)
            gb = remaining[nb]
            remaining.remove(gb)
            chosen.append(gb)
        self.MC = chosen
        self.aggOverMC=getAggQueriesOverMV(chosen,self.sel)

    def getMC(self):
        return self.MC

    def getAggOverMC(self):
        return self.aggOverMC

    def getAggOverAllLattice(self):
        return self.aggOverMV

    def getAllLattice(self):
        return self.allMVs

    def clean(self):
        dropAllMVs(self.conn)

    def getGTallLattice(self, pairs,sizeOfR,ratioViolations):
        print("Computing ground truth on all the lattice")
        dictGT={}
        H=Hypothesis()
        for p in pairs:
            sampleSize = sizeOfR

            hypothesis, hypothesisGenerationTime, samplingTime, pvalue = H.hypothesisGeneration(self.conn, p, self.sel,
                                                                                              self.measBase, self.meas,
                                                                                              self.table, sampleSize,
                                                                                              allComparison=True)
            if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
            #if len(hypothesis) == 2:
                valsToSelect = []
                for h in hypothesis:
                    valsToSelect.append(h[0])
                nbMVs = len(self.getAggOverAllLattice())

                ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(self.getAggOverAllLattice(),
                                                                                                self.sel,
                                                                                                self.measBase,
                                                                                                self.function,
                                                                                                self.table,
                                                                                                tuple(valsToSelect),
                                                                                                hypothesis,
                                                                                                self.getAllLattice())

                realRatio, nbViewOK, sizeofsample, nbInconclusive= H.checkViolations(ratioViolations, conn, ranks, hypothesis, queryCountCuboid, queryCountviolations, nbMVs)
                dictGT[p] = [nbInconclusive, realRatio]

        dictGT = utilities.sort_dict_by_second_entry_desc(dictGT)
        print("Number of comparisons found:", len(dictGT))
        return dictGT

    def getGTQueriesOverMC(self, pairs,sizeOfR,ratioViolations):
        print("Computing ground truth on queries over materialzed views")
        dictGT={}
        H=Hypothesis()
        for p in pairs:
            sampleSize = sizeOfR
            hypothesis, hypothesisGenerationTime, samplingTime,pvalue = H.hypothesisGeneration(self.conn, p, self.sel, self.measBase, self.meas,
                                                                                      self.table, sampleSize,
                                                                                      allComparison=True)
            if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
            #if len(hypothesis) == 2:
                valsToSelect = []
                for h in hypothesis:
                    valsToSelect.append(h[0])
                nbMVs = len(self.getAggOverMC())

                ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(self.getAggOverMC(),
                                                                                                self.sel,
                                                                                                self.measBase,
                                                                                                self.function,
                                                                                                self.table,
                                                                                                tuple(valsToSelect),
                                                                                                hypothesis,
                                                                                                self.getMC())

                realRatio, nbViewOK, sizeofsample, nbInconclusive = H.checkViolations(ratioViolations, conn, ranks, hypothesis, queryCountCuboid,
                                                          queryCountviolations, nbMVs)
                dictGT[p] = [nbInconclusive, realRatio]

        dictGT = utilities.sort_dict_by_second_entry_desc(dictGT)
        print("Number of comparisons found:",len(dictGT))
        return dictGT




def runComparisons():
    dfError.to_csv(fileResultsError)
    dfF1.to_csv(fileResultsF1)
    config = Config()
    s1 = Sample(conn, groupbyAtt, sel, config.meas, config.measBase, config.function, config.table)
    #s1.generateRandomMC(0.4)

    #new tests on ratioCuboidOK and ratioViolations
    for nr in tqdm(range(nbruns), desc='Runs'):
        for ratioViolations in [0.1,0.2,0.3,0.4,0.5]:
            dictGTLattice = s1.getGTallLattice(pairs, sizeOfR, ratioViolations)
    #cut here

    #dictGTLattice = s1.getGTallLattice(pairs, sizeOfR, ratioViolations)
    #dictGTMC = s1.getGTQueriesOverMC(pairs, sizeOfR, ratioViolations)


    #for nr in tqdm(range(nbruns), desc='Runs'):
        # for nr in trange(nbruns, desc='Runs'):

        #for initsampleRatio in tqdm(tabTest, desc='Init sample'):
        #for initsampleRatio in tqdm([0.01,0.1,1]):

            s1 = Sample(conn, groupbyAtt, sel, config.meas, config.measBase, config.function, config.table, 'cl')
            s1.generateRandomMC(0.4)
            dictGTMC = s1.getGTQueriesOverMC(pairs, sizeOfR, ratioViolations)

            inc=0.4
            for ratioCuboidOK in [0.4,0.6,0.8,1]:
            #for inc in [0.4]:
            #for inc in tabTest:
                # if cumulate
                s1.increaseSample(inc)
                # else
                # s1.generateSampleOfAggQueries(inc)

                dict = {}
                H=Hypothesis()
                for p in pairs:

                    # generate candidate
                    sampleSize = initsampleRatio * sizeOfR
                    hypothesis, hypothesisGenerationTime, samplingTime,pvalue = H.hypothesisGeneration(conn, p, sel, config.measBase,
                                                                                              config.meas,
                                                                                              config.table, sampleSize,
                                                                                              allComparison=True)

                    # only ok if hypothesis is a<b or a>b
                    # turn it off to check if we miss some pairs
                    # if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
                    #if len(hypothesis) == 2 and pvalue <0.01
                    if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
                        #print(hypothesis,pvalue)
                        # if True:

                        valsToSelect = []
                        for h in hypothesis:
                            valsToSelect.append(h[0])

                        ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                            s1.getCurrentSample(), sel, config.measBase, config.function, config.table, tuple(valsToSelect), hypothesis,
                            s1.getMC())

                        # compute violations over sample
                        sizeofsample = len(s1.getCurrentSample())

                        prediction, nbViewOK, sizeofsample, nbInconclusive=H.checkViolations(ratioViolations, conn, ranks, hypothesis, queryCountCuboid,queryCountviolations,sizeofsample)


                        # compute violations over ground truth
                        # first over all queries over MC
                        nbMVs = len(s1.getAggOverMC())

                        ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                            s1.getAggOverMC(), sel,
                            config.measBase,
                            config.function,
                            config.table, tuple(valsToSelect),
                            hypothesis,
                            s1.getMC())

                        gtratio, nbViewOK, nbMVs, nbInconclusive=H.checkViolations(ratioViolations, conn, ranks, hypothesis, queryCountCuboid,queryCountviolations,nbMVs)

                        errorOnMC = abs(prediction - (nbViewOK / nbMVs))

                        # then for all lattice
                        nbMVs = len(s1.getAggOverAllLattice())

                        ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                            s1.getAggOverAllLattice(), sel,
                            config.measBase,
                            config.function,
                            config.table,
                            tuple(valsToSelect),
                            hypothesis,
                            s1.getAllLattice())

                        gtratio, nbViewOK, nbMVs, nbInconclusive=H.checkViolations(ratioViolations, conn, ranks, hypothesis, queryCountCuboid,queryCountviolations,nbMVs)

                        errorOnLattice = abs(prediction - (nbViewOK / nbMVs))

                        dfError.loc[len(dfError)] = [nr, initsampleRatio, inc, p, errorOnMC, errorOnLattice, prediction,ratioCuboidOK,ratioViolations]
                        dict[p] = [errorOnLattice, prediction]
                    # flush csv

                dict = utilities.sort_dict_by_second_entry_desc(dict)
                precisionLattice, recallLattice, f1Lattice = utilities.f_measure_first_k_keys(dict, dictGTLattice, 0)
                precisionQueries, recallQueries, f1Queries = utilities.f_measure_first_k_keys(dict, dictGTMC, 0)

                pkL, rkL, fkL = utilities.f_measure_first_k_keys(dict, dictGTLattice, k)
                pkQ, rkQ, fkQ = utilities.f_measure_first_k_keys(dict, dictGTMC, k)
                dfF1.loc[len(dfF1)] = [nr, initsampleRatio, inc, precisionLattice, recallLattice, f1Lattice, rkL,
                                       precisionQueries, recallQueries, f1Queries, rkQ, k,len(dict),H.getNbWelch(),H.getNbPerm(),ratioCuboidOK,ratioViolations]

        dfError.to_csv(fileResultsError,mode='a',header=False)
        dfF1.to_csv(fileResultsF1,mode='a',header=False)
    s1.clean()


def runTimings():
    dfTimes.to_csv(fileResultsTimes)
    for nr in tqdm(range(nbruns), desc="runs"):
        s1 = Sample(conn, groupbyAtt, sel, config.meas, config.measBase, config.function, config.table)
        percentOfLattice=0.4
        s1.generateRandomMC(percentOfLattice)
        mvnames = s1.getMC()

        #generateIndex='cl'
        #for test in ['Welch','Permutation']:
        # if we want to use Claire statistics (stat), only Welch (Welch) or only permutation (Permutation) tests for generating the hypothesis
        test='stat'
        #for ratioCuboidOK in tabTest:
        #for percentOfLattice in tqdm(tabTest, desc="lattice", leave=False):
        #    s1.generateRandomMC(percentOfLattice)
        #    mvnames = s1.getMC()
        #for generateIndex in tqdm([False, True, 'mc','cl','mc-cl'], desc="index", leave=False):
        for generateIndex in tqdm([True,'mc'], desc="index", leave=False):
        #for generateIndex in ['mc-cl']:
            dropAllIndexOnMVs(conn, mvnames)
            generateIndexesOnMVs(conn, sel, mvnames, generateIndex)

            s1.generateSampleOfAggQueries(initsampleRatio)
            timings = 0
            count=0
            H = Hypothesis()
            H.setTest(test)
            for p in pairs:
            # for p in tqdm(pairs, desc="pairs", leave=False):
                start_time = time.time()

                # generate candidate
                sampleSize = initsampleRatio * sizeOfR
                #H=Hypothesis()
                hypothesis, hypothesisGenerationTime, samplingTime,pvalue = H.hypothesisGeneration(conn, p, sel,
                                                                                                  config.measBase,
                                                                                                  config.meas,
                                                                                                  config.table, sampleSize,
                                                                                                  allComparison=True)

                # only ok if hypothesis is a<b or a>b
                # turn it off to check if we miss some pairs
                if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
                #if len(hypothesis) == 2:
                    # if True:

                    valsToSelect = []
                    for h in hypothesis:
                        valsToSelect.append(h[0])

                    ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                        s1.getCurrentSample(), sel, config.measBase, config.function, config.table, tuple(valsToSelect), hypothesis,
                        s1.getMC())

                    #print(ranks)
                    # compute violations over sample
                    sizeofsample = len(s1.getCurrentSample())

                    nbViewOK,sizeofsample,nbInconclusive=H.checkViolationsOpt(ratioViolations, ratioCuboidOK, conn, ranks, hypothesis, queryCountCuboid,sizeofsample)

                    if sizeofsample == 0:
                        prediction = 0
                        bennetError = 0
                    else:
                        prediction = nbViewOK / sizeofsample
                end_time = time.time()
                timings=timings + (end_time - start_time)
                count=count+1
                dfTimes.loc[len(dfTimes)] = [nr, generateIndex, count, timings, ratioCuboidOK,percentOfLattice,test]

        dfTimes.to_csv(fileResultsTimes,mode='a',header=False)
    s1.clean()



def runTimingsByCuboids():
    dfTimes.to_csv(fileResultsTimes)
    for nr in tqdm(range(nbruns), desc="runs"):
        s1 = Sample(conn, groupbyAtt, sel, config.meas, config.measBase, config.function, config.table)
        #percentOfLattice = 0.4
        #s1.generateRandomMC(percentOfLattice)
        #mvnames = s1.getMC()

        #generateIndex='group'
        generateIndex = 'false'
        test = 'stat'

        sampleRatio=0.4
        #for percentOfLattice in tqdm(tabTest, desc="percent of lattice", leave=False):
        #for ratioCuboidOK in [0.2,0.4,0.6,0.8]:
        for ratioViolations in [0.1, 0.2,0.3,0.4,0.5]:
        #for generateIndex in tqdm([False, 'group'], desc="index", leave=False):
            for sampleRatio in tqdm(tabTest, desc="sample ratio", leave=False):
        #for test in ['Welch','Permutation']:

                s1.generateRandomMC(0.4)
                mvnames = s1.getMC()
                dropAllIndexOnMVs(conn, mvnames)
                generateIndexesOnMVs(conn, sel, mvnames, generateIndex)

                s1.generateSampleOfAggQueries(sampleRatio)
                timings = 0
                count = 0
                H = Hypothesis()
                H.setTest(test)

                start_time = time.time()

                # generate all hypothesis
                sampleSize = initsampleRatio * sizeOfR

                # generate the query for the hypothesis
                tabHypotheses=[]
                print("Generating all hypotheses")
                for p in pairs:
                    hypothesis, hypothesisGenerationTime, samplingTime, pvalue = H.hypothesisGeneration(conn, p, sel,
                                                                                                        config.measBase,
                                                                                                        config.meas,
                                                                                                        config.table, sampleSize,
                                                                                                        allComparison=True)

                    if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
                        tabHypotheses.append(hypothesis)

                #nb of queries to send=size of sample
                dictViolations={}
                for cuboidName in s1.getCurrentSample():
                    print("Validating on ", cuboidName)


                    # compute violations over current query
                    strgb = ""
                    for i in range(len(cuboidName)):
                        strgb = strgb + str(cuboidName[i])
                        if i != len(cuboidName) - 1:
                            strgb = strgb + ","
                    materialized = findMV(s1.getMC(), strgb, config.table)

                    proj=""
                    cond=""
                    for g in cuboidName[:-1]:
                        proj=proj+"c1."+str(g)+","
                        cond=cond+" and c1."+str(g) + "=c2."+str(g)

                    queryMat = ("select " +strgb+ ", "+ config.meas + "as " +config.measBase + " from \"" + materialized + "\" group by " + strgb)
                    if materialized != strgb:
                        querySigns=("select "+proj+ " c1." + sel + " as " + sel +"_1, c2." + sel + " as " + sel + "_2, "
                                    "sign(c1." + config.measBase + "- c2." + config.measBase + ") "
                                    "from ("+ queryMat + ") c1, (" + queryMat + ") c2 where c1."+sel + " < c2."+sel + cond)
                    else:
                        querySigns = ("select " + proj + " c1." + sel + " as " + sel + "_1, c2." + sel + " as " + sel + "_2, "
                                      "sign(c1." + config.measBase + "- c2." + config.measBase + ") "
                                      "from \"" + strgb + "\" c1, \"" + strgb + "\" c2 where c1." + sel + " < c2." + sel + cond)

                    queryComputeRatio=("select " + sel +"_1," + sel + "_2, (pos-neg)/cnt::float "
                                    "from (select " + sel +"_1," + sel + "_2, count(*) as cnt, count(*) filter (where sign=1) as pos,count(*) filter(where sign=-1) as neg "
                                    "from ("+querySigns + ") r group by " + sel +"_1," + sel + "_2) x")
                    #check validation on airline_code

                    queryResult=execute_query(conn, queryComputeRatio)
                    dictViolations[cuboidName]=queryResult

                #compute the scores

                dictScore={}
                for h in tabHypotheses:
                    nbViewOK = 0
                    for cuboidName in s1.getCurrentSample():
                        scoreInC=0
                        queryResult=dictViolations[cuboidName]
                        for r in queryResult:
                            if (r[0]==h[0][0] and r[1]==h[1][0]):
                                if r[2]>ratioViolations:
                                    nbViewOK=nbViewOK+1
                            if (r[0]==h[1][0] and r[1]==h[0][0]):
                                if r[2]<(-ratioViolations):
                                    nbViewOK=nbViewOK+1
                    dictScore[str(h)]=nbViewOK/len(s1.getCurrentSample())


                #print(dictScore)
                end_time = time.time()
                timings = end_time - start_time
                count = count + 1
                #dfTimes.loc[len(dfTimes)] = [nr, generateIndex, count, timings, ratioCuboidOK, sampleRatio, test]
                dfTimes.loc[len(dfTimes)] = [nr, generateIndex, count, timings, ratioCuboidOK, ratioViolations, sampleRatio, test]

        dfTimes.to_csv(fileResultsTimes,mode='a',header=False)
    s1.clean()


if __name__ == "__main__":

    # The user
    USER = "PM"
    # The DB we want
    theDB = 'F100K'
    #theDB = 'health'
    match theDB:
        case 'F9K': config = Config('configs/flightsDolap.ini', USER)
        case 'FDEBUG': config = Config('configs/flights.ini', "AC")
        case 'F100K': config = Config('configs/flights100k.ini', USER)
        case 'F600K': config = Config('configs/flightsquarterDolap.ini', USER)
        case 'F3M' : config = Config('configs/flights1923Dolap.ini', USER)
        case 'SSB': config = Config('configs/ssbDolap.ini', USER)
        case 'health': config = Config('configs/healthinsurance.ini', USER)

    # exporting results to csv
    current_time = time.localtime()
    formatted_time = time.strftime("%d-%m-%y:%H:%M:%S", current_time)
    fileResultsError = 'results/error_' + formatted_time + '_' + theDB + '.csv'
    column_namesError = ['Runs', 'Initial Sample', 'Query Sample', 'Pair', 'Error on materialized', 'Error on lattice', 'Prediction']
    fileResultsF1 = 'results/f1-r@k_' + formatted_time + '_' + theDB + '.csv'
    column_namesF1 = ['Runs', 'Initial Sample', 'Query Sample', 'Precision on Lattice', 'Recall on Lattice', 'F1 on Lattice', 'Recall@k on Lattice', 'Precision on Queries', 'Recall on Queries', 'F1 on Queries', 'Recall@k on Queries', 'k','Number of Comparisons','Number of Welch','Number of permutation']

    fileResultsTimes = 'results/times-' + formatted_time + '_' + theDB + '.csv'
    column_namesTimes = ['Runs', 'Index',  'count', 'Time', 'Ratio cuboid','sample ratio','Test']
    # testing ratio cuboids and ration violations
    #column_namesTimes = ['Runs', 'Index', 'count', 'Time', 'Ratio cuboids', 'Ratio violations', 'sample ratio','Test']
    column_namesError = ['Runs', 'Initial Sample', 'Query Sample', 'Pair', 'Error on materialized', 'Error on lattice',
                         'Prediction','Ratio cuboids', 'Ratio violations']
    column_namesF1 = ['Runs', 'Initial Sample', 'Query Sample', 'Precision on Lattice', 'Recall on Lattice',
                      'F1 on Lattice', 'Recall@k on Lattice', 'Precision on Queries', 'Recall on Queries',
                      'F1 on Queries', 'Recall@k on Queries', 'k', 'Number of Comparisons', 'Number of Welch',
                      'Number of permutation','Ratio cuboids', 'Ratio violations']

    # Create an empty DataFrame with the specified columns
    dfError = pd.DataFrame(columns=column_namesError)
    dfF1 = pd.DataFrame(columns=column_namesF1)
    dfTimes = pd.DataFrame(columns=column_namesTimes)

    conn = connect_to_db(config.dbname, config.user, config.password, config.host, config.port)

    sel = config.groupbyAtt[0]
    groupbyAtt = config.groupbyAtt[1:]

    nbpairs = 90
    paramTested = list(range(nbpairs))

    pairs = generateAllPairs(conn, sel, config.table, nbpairs)

    #s1=Sample(conn, groupbyAtt, sel, meas, measBase, function, table)
    #s1.generateRandomMC(0.4)

    # Demo code for running pairwise comparison on sample
    """
    from dbStuff import getSample_new
    #conn, measBase, table, sel, sampleSize,
    stest = getSample_new(conn, 1000)
    l = Lattice(stest)
    testing = l.pairwise(["departure_airport", "date"], "UA", "DL", "sum")
    print(testing)
    exit()
    """

    sizeOfR = getSizeOf(conn, config.table)

    initsampleRatio=0.25

    ratioCuboidOK = 0.4
    ratioViolations=0.4

    nbruns=5

    # for Recall @ k
    k = 10

    #dictGTLattice = s1.getGTallLattice(pairs, sizeOfR,ratioViolations)
    #dictGTMC = s1.getGTQueriesOverMC(pairs, sizeOfR,ratioViolations)

    #tabTest = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    #tabTest=[0.1,0.25,0.5,0.75,1]
    tabTest = [0.1, 0.25, 0.5]
    #tabTest=[0.001,0.01,0.1,0.25,0.5,0.75,1]

    comparison=True

    if comparison:
        runComparisons()
    else:
        runTimings()
        #runTimingsByCuboids()
