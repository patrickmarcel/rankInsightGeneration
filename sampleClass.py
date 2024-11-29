import configparser
import json
import random

from dolap import hypothesisGeneration,countViolations
from dbStuff import dropAllMVs,getAggQueriesOverMV,createMV,getMVnames,connect_to_db,generateAllPairs,execute_query
from bounders import generateAllqueriesOnMVs

class Sample:

    def __init__(self, conn,attInGB, selAtt, meas, function, table, generateIndex=False):
        # create all MVs and aggregates for ground truth
        self.conn=conn
        self.sel=selAtt
        dropAllMVs(conn)
        createMV(conn, attInGB, selAtt, meas, function, table, 1, generateIndex)
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
            nb = random.randint(0, len(remaining))
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

    def clean(self):
        dropAllMVs(self.conn)


if __name__ == "__main__":
    config = configparser.ConfigParser()

    # The DB we want
    config.read('configs/flightsDolap.ini')
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

    s1=Sample(conn, groupbyAtt, sel, measBase, function, table)
    s1.generateRandomMC(0.4)

    sampleSize=0.4
    ratioViolations=0.4

    for inc in [0.1,0.6,1]:
        s1.increaseSample(inc)
        #print(s1.getCurrentSample())

        for p in pairs:

            #generate candidate
            hypothesis, hypothesisGenerationTime, samplingTime = hypothesisGeneration(conn, prefs, sel, measBase, meas,
                                                                                      table, sampleSize, allComparison=True)

            # only ok if hypothesis is a<b or a>b
            if len(hypothesis) == 2 and hypothesis[0][1] != hypothesis[1][1]:

                ranks, queryCountviolations, queryCountCuboid, cuboid, newpset=generateAllqueriesOnMVs()

                #compute violations over sample

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



