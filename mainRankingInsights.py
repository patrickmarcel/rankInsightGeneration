from dbStuff import getSample_new, connect_to_db
import Config
from Lattice import Lattice
from DataSampler import DataSampler
from RankingFromPairwise import RankingFromPairwise

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
    cfg = Config.Config('configs/flights.ini','PM')
    conn = connect_to_db(cfg.dbname, cfg.user, cfg.password, cfg.host, cfg.port)
    ds=DataSampler(conn, cfg)
    sample=ds.getSample(1000,samplingMethod = 'naive')
    l = Lattice(sample)
    #testing = l.pairwise(["departure_airport", "date"], "UA", "DL", "sum")
    #print(testing)
    #print(l.getVal('DL',cfg.measBase))
    #print(testing.columns)
    #valA=testing[(testing['airline'] == 'DL')]
    #print(valA)

    r=len(cfg.groupbyAtt)-1
    p=1
    ranking=RankingFromPairwise(cfg.prefs, r,p)
    ranking.run(l)