import math
from math import hypot

from dbStuff import getSample_new, connect_to_db, execute_query, getCuboidsOfAtt
import Config
from Lattice import Lattice
from DataSampler import DataSampler
from RankingFromPairwise import RankingFromPairwise
import time
from sampleClassRanking import SampleRanking

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

if __name__ == '__main__':

#    numpy.seterr(all='raise')
    cfg = Config.Config('configs/flights100k.ini','PM')
    conn = connect_to_db(cfg.dbname, cfg.user, cfg.password, cfg.host, cfg.port)

    adom = [x[0] for x in execute_query(conn, "select distinct  " + cfg.sel + " from " + cfg.table + ";")]
    #table_size = execute_query(conn, "select count(1) from " + cfg.table + ";")[0][0]

    start_time = time.time()
    ds=DataSampler(conn, cfg)
    sample=ds.getSample(10358,samplingMethod = 'naive')
    l = Lattice(sample)
    #testing = l.pairwise(["departure_airport", "date"], "UA", "DL", "sum")
    #print(testing)
    #print(l.getVal('DL',cfg.measBase))
    #print(testing.columns)
    #valA=testing[(testing['airline'] == 'DL')]
    #print(valA)

    r=int(math.pow(2,len(cfg.groupbyAtt)-1))
    #to increase the chances of gaps in deltak, enabling drawing with replacement
    r=r*1
    #print('r:',r)
    p=1
    #ranking=RankingFromPairwise(cfg.prefs, r,p)
    ranking=RankingFromPairwise(adom, r,p, 'Welch', True)
    ranking.run(l,method='withoutTest')
    #print('Delta:',ranking.delta)
    print('F:',ranking.F)
    #print('Tau:',ranking.tau)
    #print('M',ranking.M)
    end_time = time.time()
    timings = end_time - start_time
    print('Completed in ',timings, 'seconds')
    hypothesis=ranking.getHypothesis()
    print('Hypothesis:',hypothesis)

    groupbyAtt = cfg.groupbyAtt[1:]
    #sampleOfLattice=SampleRanking(conn, groupbyAtt, cfg.sel, cfg.meas, cfg.measBase, cfg.function, cfg.table, generateIndex=False)

