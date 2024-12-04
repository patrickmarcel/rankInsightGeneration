import configparser
import json
import random
import time
import pandas as pd
from tqdm import tqdm, trange
import utilities

from dolap import hypothesisGeneration,countViolations
from dbStuff import dropAllMVs,getAggQueriesOverMV,createMV,getMVnames,connect_to_db,generateAllPairs,execute_query,getSizeOf,dropAllIndexOnMVs,generateIndexesOnMVs
from bounders import generateAllqueriesOnMVs

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
        createMV(conn, attInGB, selAtt, measBase, function, table, 1, generateIndex)
        self.allMVs=getMVnames(conn)
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
        remaining = self.aggOverMC.copy()
        chosen = []
        for j in range(sizeOfSample):
            nb = random.randint(0, len(remaining)-1)
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
        dictGT={}
        for p in pairs:
            sampleSize = sizeOfR
            hypothesis, hypothesisGenerationTime, samplingTime,pvalue = hypothesisGeneration(self.conn, p, self.sel, self.measBase, self.meas,
                                                                                      self.table, sampleSize,
                                                                                      allComparison=True)
            if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
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

                nbInconclusive = 0
                nbViewOK = 0
                for i in range(len(queryCountviolations)):

                    v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
                    c = execute_query(conn, queryCountCuboid[i])[0][0]

                    if c != 0:
                        if ratio < ratioViolations:
                            nbViewOK = nbViewOK + 1
                    else:
                        nbMVs = nbMVs - 1
                        nbInconclusive = nbInconclusive + 1

                realRatio = nbViewOK / nbMVs

                dictGT[p] = [False, realRatio]

        dictGT = utilities.sort_dict_by_second_entry_desc(dictGT)
        #print(dictGT)
        return dictGT

    def getGTQueriesOverMC(self, pairs,sizeOfR,ratioViolations):
        dictGT={}
        for p in pairs:
            sampleSize = sizeOfR
            hypothesis, hypothesisGenerationTime, samplingTime,pvalue = hypothesisGeneration(self.conn, p, self.sel, self.measBase, self.meas,
                                                                                      self.table, sampleSize,
                                                                                      allComparison=True)
            if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
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

                nbInconclusive = 0
                nbViewOK = 0
                for i in range(len(queryCountviolations)):

                    v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
                    c = execute_query(conn, queryCountCuboid[i])[0][0]

                    if c != 0:
                        if ratio < ratioViolations:
                            nbViewOK = nbViewOK + 1
                    else:
                        nbMVs = nbMVs - 1
                        nbInconclusive = nbInconclusive + 1

                realRatio = nbViewOK / nbMVs

                dictGT[p] = [False, realRatio]

        dictGT = utilities.sort_dict_by_second_entry_desc(dictGT)

        return dictGT




def runComparisons():
    s1 = Sample(conn, groupbyAtt, sel, meas, measBase, function, table)
    s1.generateRandomMC(0.4)

    dictGTLattice = s1.getGTallLattice(pairs, sizeOfR, ratioViolations)
    dictGTMC = s1.getGTQueriesOverMC(pairs, sizeOfR, ratioViolations)


    for nr in tqdm(range(nbruns)):
        # for nr in trange(nbruns, desc='Runs'):

        for initsampleRatio in tabTest:

            s1 = Sample(conn, groupbyAtt, sel, meas, measBase, function, table, 'mc')
            s1.generateRandomMC(0.4)

            for inc in tabTest:
                # if cumulate
                s1.increaseSample(inc)
                # else
                # s1.generateSampleOfAggQueries(inc)

                dict = {}

                for p in pairs:

                    # generate candidate
                    sampleSize = initsampleRatio * sizeOfR
                    hypothesis, hypothesisGenerationTime, samplingTime,pvalue = hypothesisGeneration(conn, p, sel, measBase,
                                                                                              meas,
                                                                                              table, sampleSize,
                                                                                              allComparison=True)

                    # only ok if hypothesis is a<b or a>b
                    # turn it off to check if we miss some pairs
                    # if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
                    #if len(hypothesis) == 2 and pvalue <0.01
                    if len(hypothesis) == 2and hypothesis[0][1] != hypothesis[1][1] and pvalue >0.01:
                        #print(hypothesis,pvalue)
                        # if True:

                        valsToSelect = []
                        for h in hypothesis:
                            valsToSelect.append(h[0])

                        ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                            s1.getCurrentSample(), sel, measBase, function, table, tuple(valsToSelect), hypothesis,
                            s1.getMC())

                        # compute violations over sample
                        nbViewOK = 0
                        nbInconclusive = 0
                        sizeofsample = len(s1.getCurrentSample())

                        for i in range(len(ranks)):
                            v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
                            c = execute_query(conn, queryCountCuboid[i])[0][0]

                            if c != 0:
                                if ratio < ratioViolations:
                                    nbViewOK = nbViewOK + 1
                            else:
                                sizeofsample = sizeofsample - 1
                                nbInconclusive = nbInconclusive + 1

                        if sizeofsample == 0:
                            prediction = 0
                            bennetError = 0
                        else:
                            prediction = nbViewOK / sizeofsample

                        # compute violations over ground truth
                        # first over all queries over MC
                        nbMVs = len(s1.getAggOverMC())

                        ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                            s1.getAggOverMC(), sel,
                            measBase,
                            function,
                            table, tuple(valsToSelect),
                            hypothesis,
                            s1.getMC())

                        nbInconclusive = 0
                        nbViewOK = 0
                        for i in range(len(queryCountviolations)):

                            v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
                            c = execute_query(conn, queryCountCuboid[i])[0][0]

                            if c != 0:
                                if ratio < ratioViolations:
                                    nbViewOK = nbViewOK + 1
                            else:
                                nbMVs = nbMVs - 1
                                nbInconclusive = nbInconclusive + 1

                        gtratio = nbViewOK / nbMVs
                        errorOnMC = abs(prediction - (nbViewOK / nbMVs))

                        # then for all lattice
                        nbMVs = len(s1.getAggOverAllLattice())

                        ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                            s1.getAggOverAllLattice(), sel,
                            measBase,
                            function,
                            table,
                            tuple(valsToSelect),
                            hypothesis,
                            s1.getAllLattice())

                        nbInconclusive = 0
                        nbViewOK = 0
                        for i in range(len(queryCountviolations)):

                            v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
                            c = execute_query(conn, queryCountCuboid[i])[0][0]

                            if c != 0:
                                if ratio < ratioViolations:
                                    nbViewOK = nbViewOK + 1
                            else:
                                nbMVs = nbMVs - 1
                                nbInconclusive = nbInconclusive + 1

                        gtratio = nbViewOK / nbMVs
                        errorOnLattice = abs(prediction - (nbViewOK / nbMVs))

                        dfError.loc[len(dfError)] = [nr, initsampleRatio, inc, p, errorOnMC, errorOnLattice, prediction]
                        dict[p] = [errorOnLattice, prediction]
                    # flush csv
                dict = utilities.sort_dict_by_second_entry_desc(dict)
                precisionLattice, recallLattice, f1Lattice = utilities.f_measure_first_k_keys(dict, dictGTLattice, 0)
                precisionQueries, recallQueries, f1Queries = utilities.f_measure_first_k_keys(dict, dictGTMC, 0)

                pkL, rkL, fkL = utilities.f_measure_first_k_keys(dict, dictGTLattice, k)
                pkQ, rkQ, fkQ = utilities.f_measure_first_k_keys(dict, dictGTMC, k)
                dfF1.loc[len(dfF1)] = [nr, initsampleRatio, inc, precisionLattice, recallLattice, f1Lattice, rkL,
                                       precisionQueries, recallQueries, f1Queries, rkQ, k]

    dfError.to_csv(fileResultsError)
    dfF1.to_csv(fileResultsF1)
    s1.clean()


def runTimings():

    for nr in tqdm(range(nbruns)):
        s1 = Sample(conn, groupbyAtt, sel, meas, measBase, function, table)
        s1.generateRandomMC(0.4)
        mvnames = s1.getMC()

        for generateIndex in [False, True, 'mc']:
            dropAllIndexOnMVs(conn, mvnames)
            generateIndexesOnMVs(conn, sel, mvnames, generateIndex)

            s1.generateSampleOfAggQueries(initsampleRatio)
            timings = 0
            count=0
            for p in pairs:
                start_time = time.time()

                # generate candidate
                sampleSize = initsampleRatio * sizeOfR
                hypothesis, hypothesisGenerationTime, samplingTime,pvalue = hypothesisGeneration(conn, p, sel,
                                                                                                  measBase,
                                                                                                  meas,
                                                                                                  table, sampleSize,
                                                                                                  allComparison=True)

                # only ok if hypothesis is a<b or a>b
                # turn it off to check if we miss some pairs
                # if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:
                if len(hypothesis) == 2:
                    # if True:

                    valsToSelect = []
                    for h in hypothesis:
                        valsToSelect.append(h[0])

                    ranks, queryCountviolations, queryCountCuboid, cuboid = generateAllqueriesOnMVs(
                        s1.getCurrentSample(), sel, measBase, function, table, tuple(valsToSelect), hypothesis,
                        s1.getMC())


                    # compute violations over sample
                    nbViewOK = 0
                    nbInconclusive = 0
                    sizeofsample = len(s1.getCurrentSample())

                    for i in range(len(ranks)):
                        v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
                        c = execute_query(conn, queryCountCuboid[i])[0][0]

                        if c != 0:
                            if ratio < ratioViolations:
                                nbViewOK = nbViewOK + 1
                        else:
                            sizeofsample = sizeofsample - 1
                            nbInconclusive = nbInconclusive + 1

                    if sizeofsample == 0:
                        prediction = 0
                        bennetError = 0
                    else:
                        prediction = nbViewOK / sizeofsample
                end_time = time.time()
                timings=timings + (end_time - start_time)
                count=count+1
                dfTimes.loc[len(dfTimes)] = [nr, generateIndex, count, timings]

    dfTimes.to_csv(fileResultsTimes)
    s1.clean()



if __name__ == "__main__":
    # exporting results to csv
    current_time = time.localtime()
    formatted_time = time.strftime("%d-%m-%y:%H:%M:%S", current_time)
    fileResultsError = 'results/error_' + formatted_time + '.csv'
    column_namesError = ['Runs', 'Initial Sample', 'Query Sample', 'Pair', 'Error on materialized', 'Error on lattice', 'Prediction']
    fileResultsF1 = 'results/f1-r@k_' + formatted_time + '.csv'
    column_namesF1 = ['Runs', 'Initial Sample', 'Query Sample', 'Precision on Lattice', 'Recall on Lattice', 'F1 on Lattice', 'Recall@k on Lattice', 'Precision on Queries', 'Recall on Queries', 'F1 on Queries', 'Recall@k on Queries', 'k']

    fileResultsTimes = 'results/times-' + formatted_time + '.csv'
    column_namesTimes = ['Runs', 'Index',  'count', 'Time']

    # Create an empty DataFrame with the specified columns
    dfError = pd.DataFrame(columns=column_namesError)
    dfF1 = pd.DataFrame(columns=column_namesF1)
    dfTimes = pd.DataFrame(columns=column_namesTimes)

    config = configparser.ConfigParser()

    # The DB we want
    config.read('configs/flightsDolap.ini')
    #config.read('configs/flightsquarterDolap.ini')
    #config.read('configs/ssbDolap.ini')
    # config.read('configs/flights1923Dolap.ini')
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

    conn = connect_to_db(dbname, user, password, host, port)

    sel = groupbyAtt[0]
    groupbyAtt = groupbyAtt[1:]

    nbpairs = 90
    paramTested = list(range(nbpairs))

    pairs = generateAllPairs(conn, sel, table, nbpairs)

    #s1=Sample(conn, groupbyAtt, sel, meas, measBase, function, table)
    #s1.generateRandomMC(0.4)

    sizeOfR = getSizeOf(conn, table)

    initsampleRatio=0.4

    ratioViolations=0.4
    nbruns=5

    # for Recall @ k
    k = 10

    #dictGTLattice = s1.getGTallLattice(pairs, sizeOfR,ratioViolations)
    #dictGTMC = s1.getGTQueriesOverMC(pairs, sizeOfR,ratioViolations)

    tabTest = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    #tabTest=[0.6]

    comparison=True

    if comparison:
        runComparisons()
    else:
        runTimings()
