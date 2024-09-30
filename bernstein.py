import math
import random
import dbStuff

# delta is the probabiliy of making an error, sigma is the variance and t the error
def sizeOfSample(delta, sigma, t):
    return (math.log(2/delta) * (2* math.pow(sigma,2) + (2*t)/3)) / math.pow(t,2)

def sizeOfSampleHoeffding(delta, t):
    return math.log(2 / delta, 10) / (2* math.pow(t,2))

# probability of making a mistake for error t
def bernsteinBound(delta, sigma, t):
    return 2* math.exp(-(math.pow(t,2) /2)/(math.pow(sigma,2) +t/3))

#delta the probability of making an error
def bersteinError(delta, sigma):
    return math.sqrt(math.ln(2/delta)*sigma/2) + (math.ln(2/delta)*(1+math.sqrt(2)))/math.sqrt(2)*6

# so far name of table in from is name of cuboid (convention: attribute sorted + sel attribute last)
# todo search for materialized cuboids that is the closest to the one drawn (using substring should be enough)
# pwrset is the powerset of categorical attributes that include the target attribute
def getSample(delta, t, pwrset, sel, meas, function, table, valsToSelect, hypothesis):
    pset=pwrset
    n=int(sizeOfSampleHoeffding(delta, t))
    tabQuery=[]
    tabCount=[]
    tabCuboid=[]
    hyp = ""
    for i in range(len(hypothesis)):
        hyp = hyp + str(hypothesis[i])
        if i != len(hypothesis) - 1:
            hyp = hyp + ","
    for i in range(n):
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
        if strgb == sel:
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 # " FROM " + table +
                 " FROM \"" + strgb + "\"" +
                 " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
        else:
            q = ("SELECT " + strgb + ", " + function + '(' + meas + "), "
                 + " rank () over ( partition by " + gbwithoutsel + " order by " + function + '(' + meas + ") desc ) as rank" +
                 #" FROM " + table +
             " FROM \"" + strgb + "\"" +
                " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + " ")
        queryHyp = (
                "select " + sel + ",rank  from (values " + hyp + ") as t2 (" + sel + ",rank)")
        queryExcept = ("select * from  (" + q + " ) t3 , (" + queryHyp + ") t4 where t3." + sel + "=t4." + sel + " and t3.rank!=t4.rank")
        queryCountCuboid= ("select count(*) from (" + q + ") t5;")
        queryCountExcept = ("select count(*) from (" + queryExcept + ") t6;")

        """
        q = ("SELECT " + strgb + "," + sel + "," + meas + ", "
                 + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
                 #" FROM " + table +
             " FROM \"" + strgb + "\"" +
                " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ")
        queryHyp = (
                "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")
        queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + q + " ) t3 except all " + queryHyp + " ")
        queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")
        """
        tabQuery.append(queryCountExcept)
        tabCount.append(queryCountCuboid)
        tabCuboid.append(strgb)
    return tabQuery,tabCount,tabCuboid


def runSampleQueries(tabQuery):
    tabViolations=[]
    for tupleOfQueries in tabQuery:
        #getMVfor(q)
        #execon(q,MV)
        #for now just run queries over R
        # compute violations of q
        r=dbStuff.execute_query(tupleOfQueries[3])
        tabViolations.append(r[0])

