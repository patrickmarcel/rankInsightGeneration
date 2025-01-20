import math

from dbStuff import getSample_new, connect_to_db, execute_query, getCuboidsOfAtt
import Config
from Lattice import Lattice
from DataSampler import DataSampler
from RankingFromPairwise import RankingFromPairwise
import numpy

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

    ds=DataSampler(conn, cfg)
    sample=ds.getSample(60000,samplingMethod = 'naive')
    l = Lattice(sample)
    #testing = l.pairwise(["departure_airport", "date"], "UA", "DL", "sum")
    #print(testing)
    #print(l.getVal('DL',cfg.measBase))
    #print(testing.columns)
    #valA=testing[(testing['airline'] == 'DL')]
    #print(valA)

    r=int(math.pow(2,len(cfg.groupbyAtt)-1))
    #to increase the chances of gaps in deltak, enabling drawing with replacement
    r=r*10
    print('r:',r)
    p=1
    #ranking=RankingFromPairwise(cfg.prefs, r,p)
    ranking=RankingFromPairwise(adom, r,p, 'Welch', True)
    ranking.run(l)
    print('Delta:',ranking.delta)
    print('F:',ranking.F)
    #print('Tau:',ranking.tau)
    #print('M',ranking.M)

