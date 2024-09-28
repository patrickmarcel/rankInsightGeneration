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

# pwrset is the powerset of categorical attributes that include the target attribute
def getSample(delta, t, pwrset, sel, meas, table, valsToSelect, hypothesis):
    pset=pwrset
    n=sizeOfSampleHoeffding(delta, t)
    tabQuery=[]
    hyp = ""
    for i in range(len(hypothesis)):
        hyp = hyp + str(hypothesis[i])
        if i != len(hypothesis) - 1:
            hyp = hyp + ","
    for i in range(n):
        nb = random.randint(0, len(pwrset) - 1)
        gb = pset[nb]
        pset.remove(nb)
        strgb = ""
        for i in range(len(gb)):
            strgb = strgb + str(gb[i])
            if i != len(gb) - 1:
                strgb = strgb + ","
        q = ("SELECT " + strgb + "," + sel + "," + meas + ", "
                 + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
                 " FROM " + table + " WHERE " + sel + " in " + str(
                    valsToSelect) + " group by " + strgb + "," + sel + " ")
        queryHyp = (
                "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")
        queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + q + " ) t3 except all " + queryHyp + " ")
        queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")

        tabQuery.append((q,queryHyp,queryExcept,queryCountExcept))
    return tabQuery


def runSampleQueries(tabQuery):
    tabViolations=[]
    for tupleOfQueries in tabQuery:
        #getMVfor(q)
        #execon(q,MV)
        #for now just run queries over R
        # compute violations of q
        r=dbStuff.execute_query(tupleOfQueries[3])
        tabViolations.append(r[0])

