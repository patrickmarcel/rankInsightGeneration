import math
import random
import dbStuff

# delta is the probabiliy of making an error, sigma is the variance and t the error
def sizeOfSample(delta, sigma, t):
    return (math.log(2/delta) * (2* math.pow(sigma,2) + (2*t)/3)) / math.pow(t,2)

def sizeOfSampleHoeffding(delta, t):
    return math.log(2 / delta) / (2* math.pow(t,2))

# for n<=N/2
def sizeOfSampleHoeffdingSerflingFromBardenet(delta, t, N):
    X=(2*N*math.log(1/delta))/(math.pow(t,2))
    return (X*N - 1) / (1+X)



# probability of making a mistake for error t
def bernsteinBound(sigma, t):
    return 2* math.exp(- (math.pow(t,2) /2) / (sigma +t/3) )

#delta the probability of making an error
def bersteinError(delta, sigma):
    return math.sqrt( math.log(2/delta)*sigma/2 ) + ( math.log(2/delta) * (1+math.sqrt(2)) ) /(math.sqrt(2)*6)

#delta the probability of making an error
# n sample size
def bennetErrorOnAvg(delta, sigma, n):
    return math.sqrt( math.log(1/delta)*sigma*2/n ) + ( math.log(1/delta) ) /(3 * n)


def empiricalBennetFromMaurer(delta, sigma, n):
    return math.sqrt(math.log(2 / delta) * sigma * 2 / n) + (math.log(2 / delta)*7) / (3 * (n-1))

#N is size of population
def empiricalBernsteinFromBardenet(delta, sigma, n, N):
    if(n<=N/2):
        return   math.sqrt(sigma *  math.log(1 / delta) *  2 * (1-((n-1)/N)) / n) + ( math.log(1 / delta) * (7/3 + 3/math.sqrt(2))) / n

    else:
        return  math.sqrt(sigma *  math.log(1 / delta) *  2 * ((1-(n/N))*(1+(1/n)) )  / n) + ( math.log(1 / delta) * (7/3 + 3/math.sqrt(2))) / n



# checks if gb in some view names
# returns the closest one or table if not found
def findMV(mvnames, gb, table):
    #mvnames = dbStuff.getMVnames(conn)
    min=10000
    result=table
    tabgb = gb.split(",")
    for n in mvnames:
        ok = True
        for t in tabgb:
            if t not in n[0]:
                ok=False
        if ok and len(n[0])<min:
            result=n[0]
            min=len(n[0])
    return(result)

# so far name of table in from is name of cuboid (convention: attribute sorted + sel attribute last)
# pwrset is the powerset of categorical attributes that include the target attribute
def getSample(pwrset, sel, meas, function, table, valsToSelect, hypo, mvnames, withReplacement=False, withBias=False, sizeOfSample=1):
    #### remove fact table cuboid from power set
    #print("pwrset in getSample: ",pwrset)
    #del pwrset[-1]
    pset=pwrset.copy()
    if withBias:
        tabNb=[]
        tabWeights=[]
        sumLength=0
        for p in pwrset:
            sumLength=sumLength+len(p)
        maxlength=len(pwrset[-1])+1
        i=0
        for p in pwrset:
            w=(maxlength-len(p))/sumLength
            tabWeights.append(w)
            tabNb.append(i)
            i=i+1
        #print(tabWeights)
    #print(math.fsum(tabWeights))
    #n=int(sizeOfSampleHoeffding(delta, t))
    n=sizeOfSample
    tabRanks=[]
    tabQuery=[]
    tabCount=[]
    tabCuboid=[]
    hyp = ""
    for i in range(len(hypo)):
        hyp = hyp + str(hypo[i])
        if i != len(hypo) - 1:
            hyp = hyp + ","
    for j in range(n):
        if withBias:
            nb=max(random.choices(tabNb,tabWeights))
        else:
            #nb = random.randint(0, len(pwrset) - 1)
            nb = random.randint(0, len(pset) - 1)
        gb = pset[nb]
        if withReplacement==False:
            # without replacement: gb is removed from the list so as not to be drawn twice
            pset.remove(gb)
        strgb = ""
        gbwithoutsel=""
        for i in range(len(gb)):
            strgb = strgb + str(gb[i])
            if i != len(gb) - 1:
                strgb = strgb + ","
        for i in range(len(gb)-1):
            gbwithoutsel = gbwithoutsel + str(gb[i])
            if i != len(gb) - 2:
                gbwithoutsel = gbwithoutsel + ","
        materialized = findMV(mvnames, strgb, table)
        #print(materialized)
        # if materialized = table then reject query???
        queryValidGB=''
        queryHyp = (
                "select " + sel + ",rank  from (values " + hyp + ") as t2 (" + sel + ",rank)")
        if strgb == sel:
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 " FROM \"" + str(materialized) + "\"" +
                 " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
            queryExcept = (
                        "select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            queryCountCuboid = ("select count(*) from (" + q + ") t5;")
            queryRank=("select 'all',string_agg(" + sel + "::text,',') from (" + q + "  order by rank);")
        else:
            queryValidGB = ("select " + gbwithoutsel + " from (SELECT " + strgb + " FROM \"" + str(
                materialized) + "\"" + " WHERE " + sel + " in " + str(
                valsToSelect) + " group by " + strgb + " ) t group by " + gbwithoutsel + " having count(*) >1")
            #queryValidGB = ("select " + gbwithoutsel + " from (SELECT " + strgb + " FROM \"" + str(
            #    materialized) + "\"" + " WHERE " + sel + " in " + str(
            #    valsToSelect) + " group by " + strgb + " order by rank) t group by " + gbwithoutsel + " having count(*) >1")
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( partition by " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
             " FROM \"" + str(materialized) + "\"" +
                " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
            queryExcept = (
                        "select * from  (select * from  (" + q + " ) t7 where ("+ gbwithoutsel + ") in ("+ queryValidGB + ")) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            queryRank = ("select " + gbwithoutsel + ",string_agg(" + sel + "::text,',') from (" + q + " ) t7 where ("+ gbwithoutsel + ") in ("+ queryValidGB + ") group by "+gbwithoutsel +";")

            #queryExcept = (
            #        "select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            #print('queryExcept:',queryExcept)
            queryCountCuboid = ("select count(*) from (Select * from(" + q + ") t8 where ("+ gbwithoutsel + ") in ("+ queryValidGB + ")) t5;")
        #print('queryCountCuboid:',queryCountCuboid)
        queryCountExcept = ("select count(*) from (" + queryExcept + ") t6;")
        #print('queryCountExcept:',queryCountExcept)

        #print('query rank:', queryRank)

        tabRanks.append(queryRank)
        tabQuery.append(queryCountExcept)
        tabCount.append(queryCountCuboid)
        tabCuboid.append(strgb)
    return tabRanks,tabQuery,tabCount,tabCuboid,pset

def getMoreRandomQueries(sizeofquerysample,currentSample, sel, meas, function,table,valsToSelect,hypo,mvnames, withReplacement=False, withBias=False):
    currentsize=len(currentSample["ranks"])
    #print("sizeofquerysample - currentsize:",sizeofquerysample - currentsize)

    tabRanks,tabQuery,tabCount,tabCuboid,pset=getSample(currentSample["pset"], sel, meas, function, table, valsToSelect, hypo, mvnames, withReplacement=False, withBias=False, sizeOfSample=sizeofquerysample - currentsize)
    return tabRanks,tabQuery,tabCount,tabCuboid,pset

def generateAllqueriesOnMVs(pwrset, sel, meas, function, table, valsToSelect, hypo, mvnames):
    pset=pwrset
    n=len(pwrset)
    tabRanks=[]
    tabQuery=[]
    tabCount=[]
    tabCuboid=[]
    hyp = ""
    for i in range(len(hypo)):
        hyp = hyp + str(hypo[i])
        if i != len(hypo) - 1:
            hyp = hyp + ","
    for j in range(n):
        #nb = random.randint(0, len(pwrset) - 1)
        gb = pset[j]
        # without replacement: gb is removed from the list so as not to be drawn twice
        #pset.remove(gb)
        strgb = ""
        gbwithoutsel=""
        for i in range(len(gb)):
            strgb = strgb + str(gb[i])
            if i != len(gb) - 1:
                strgb = strgb + ","
        for i in range(len(gb)-1):
            gbwithoutsel = gbwithoutsel + str(gb[i])
            if i != len(gb) - 2:
                gbwithoutsel = gbwithoutsel + ","
        materialized = findMV(mvnames, strgb, table)
        #print("STRGB: ",strgb)
        #print("MATERIALIZED: ",materialized)
        #materialized=strgb
        queryHyp = (
                "select " + sel + ",rank  from (values " + hyp + ") as t2 (" + sel + ",rank)")
        if strgb == sel:
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 " FROM \"" + str(materialized) + "\"" +
                 " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
            queryExcept = (
                    "select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            queryCountCuboid = ("select count(*) from (" + q + ") t5;")
            queryRank = ("select 'all',string_agg(" + sel + "::text,',') from (" + q + "  order by rank);")
        else:
            queryValidGB = ("select " + gbwithoutsel + " from (SELECT " + strgb + " FROM \"" + str(
                materialized) + "\"" + " WHERE " + sel + " in " + str(
                valsToSelect) + " group by " + strgb + " ) t group by " + gbwithoutsel + " having count(*) >1")
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( partition by " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 " FROM \"" + str(materialized) + "\"" +
                 " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
            queryExcept = (
                    "select * from  (select * from  (" + q + " ) t7 where (" + gbwithoutsel + ") in (" + queryValidGB + ")) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            queryRank = (
                        "select " + gbwithoutsel + ",string_agg(" + sel + "::text,',') from (" + q + " ) t7 where (" + gbwithoutsel + ") in (" + queryValidGB + ") group by " + gbwithoutsel + ";")

            # queryExcept = (
            #        "select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            # print('queryExcept:',queryExcept)
            queryCountCuboid = (
                        "select count(*) from (Select * from(" + q + ") t8 where (" + gbwithoutsel + ") in (" + queryValidGB + ")) t5;")
        # print('queryCountCuboid:',queryCountCuboid)
        queryCountExcept = ("select count(*) from (" + queryExcept + ") t6;")
        # print('queryCountExcept:',queryCountExcept)

        # print('query rank:', queryRank)

        tabRanks.append(queryRank)
        tabQuery.append(queryCountExcept)
        tabCount.append(queryCountCuboid)
        tabCuboid.append(strgb)
    return tabRanks, tabQuery, tabCount, tabCuboid


def generateAllqueries(pwrset, sel, meas, function, table, valsToSelect, hypo, mvnames):
    pset=pwrset
    n=len(pwrset)
    tabRanks=[]
    tabQuery=[]
    tabCount=[]
    tabCuboid=[]
    hyp = ""
    for i in range(len(hypo)):
        hyp = hyp + str(hypo[i])
        if i != len(hypo) - 1:
            hyp = hyp + ","
    for j in range(n):
        nb = random.randint(0, len(pwrset) - 1)
        gb = pset[nb]
        # without replacement: gb is removed from the list so as not to be drawn twice
        pset.remove(gb)
        strgb = ""
        gbwithoutsel=""
        for i in range(len(gb)):
            strgb = strgb + str(gb[i])
            if i != len(gb) - 1:
                strgb = strgb + ","
        for i in range(len(gb)-1):
            gbwithoutsel = gbwithoutsel + str(gb[i])
            if i != len(gb) - 2:
                gbwithoutsel = gbwithoutsel + ","
        #materialized = findMV(mvnames, strgb, table)
        materialized=strgb
        queryHyp = (
                "select " + sel + ",rank  from (values " + hyp + ") as t2 (" + sel + ",rank)")
        if strgb == sel:
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 " FROM \"" + str(materialized) + "\"" +
                 " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
            queryExcept = (
                    "select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            queryCountCuboid = ("select count(*) from (" + q + ") t5;")
            queryRank = ("select 'all',string_agg(" + sel + "::text,',') from (" + q + "  order by rank);")
        else:
            queryValidGB = ("select " + gbwithoutsel + " from (SELECT " + strgb + " FROM \"" + str(
                materialized) + "\"" + " WHERE " + sel + " in " + str(
                valsToSelect) + " group by " + strgb + " ) t group by " + gbwithoutsel + " having count(*) >1")
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( partition by " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 " FROM \"" + str(materialized) + "\"" +
                 " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
            queryExcept = (
                    "select * from  (select * from  (" + q + " ) t7 where (" + gbwithoutsel + ") in (" + queryValidGB + ")) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            queryRank = (
                        "select " + gbwithoutsel + ",string_agg(" + sel + "::text,',') from (" + q + " ) t7 where (" + gbwithoutsel + ") in (" + queryValidGB + ") group by " + gbwithoutsel + ";")

            # queryExcept = (
            #        "select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
            # print('queryExcept:',queryExcept)
            queryCountCuboid = (
                        "select count(*) from (Select * from(" + q + ") t8 where (" + gbwithoutsel + ") in (" + queryValidGB + ")) t5;")
        # print('queryCountCuboid:',queryCountCuboid)
        queryCountExcept = ("select count(*) from (" + queryExcept + ") t6;")
        # print('queryCountExcept:',queryCountExcept)

        # print('query rank:', queryRank)

        tabRanks.append(queryRank)
        tabQuery.append(queryCountExcept)
        tabCount.append(queryCountCuboid)
        tabCuboid.append(strgb)
    return tabRanks, tabQuery, tabCount, tabCuboid


def runSampleQueries(tabQuery):
    tabViolations=[]
    for tupleOfQueries in tabQuery:
        #getMVfor(q)
        #execon(q,MV)
        #for now just run queries over R
        # compute violations of q
        r=dbStuff.execute_query(tupleOfQueries[3])
        tabViolations.append(r[0])

