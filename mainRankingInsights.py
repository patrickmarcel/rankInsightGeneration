import math
from math import hypot

from tqdm import tqdm

from dbStuff import getSample_new, connect_to_db, execute_query, getCuboidsOfAtt, getSizeOf
import Config
from Lattice import Lattice
from DataSampler import DataSampler
from RankingFromPairwise import RankingFromPairwise
import time
from sampleClassRanking import SampleRanking
from statStuff import transform_to_rankings, compute_kendall_tau
import pandas as pd

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

    #parameters

    #samplingMethod='naive'
    samplingMethod='congressional'
    #groundTruth = ['WN', 'AA', 'DL', 'OO', 'UA', 'NK', '9E', 'YX', 'MQ', 'YV', 'OH', 'B6', 'F9', 'G4', 'AS', 'HA']
    groundTruth=['G4', 'YV', 'AA', 'OO', 'NK', 'MQ', 'UA', 'WN', 'DL', 'OH', 'B6', '9E', 'YX', 'HA', 'F9', 'AS']
    computeGT = False
    computeHyp = True
    #samplesize=10358
    #config='configs/flights100k.ini'
    user='PM'
    theDB = 'F100K'

    #comparison method
    #method = 'withoutTest'
    method = 'withTest'

    match theDB:
        case 'F9K': cfg = Config.Config('configs/flightsDolap.ini', user)
        case 'FDEBUG': cfg = Config.Config('configs/flights.ini', user)
        case 'F100K': cfg = Config.Config('configs/flights100k.ini', user)
        case 'F600K': cfg = Config.Config('configs/flightsquarterDolap.ini', user)
        case 'F3M' : cfg = Config.Config('configs/flights1923Dolap.ini', user)
        case 'SSB': cfg = Config.Config('configs/ssbDolap.ini', user)


    # exporting results to csv
    current_time = time.localtime()
    formatted_time = time.strftime("%d-%m-%y:%H:%M:%S", current_time)
    fileResultsError = 'results/error_' + formatted_time + '_' + theDB + '.csv'
    column_namesError = ['r','samplesize','k','tau']
    #fileResultsF1 = 'results/f1-r@k_' + formatted_time + '_' + theDB + '.csv'
    #column_namesF1 = ['Runs', 'Initial Sample', 'Query Sample', 'Precision on Lattice', 'Recall on Lattice', 'F1 on Lattice', 'Recall@k on Lattice', 'Precision on Queries', 'Recall on Queries', 'F1 on Queries', 'Recall@k on Queries', 'k','Number of Comparisons','Number of Welch','Number of permutation']
    #fileResultsTimes = 'results/times-' + formatted_time + '_' + theDB + '.csv'
    #column_namesTimes = ['Runs', 'Index',  'count', 'Time', 'Ratio cuboid','sample ratio','Test']

    # Create an empty DataFrame with the specified columns
    dfError = pd.DataFrame(columns=column_namesError)
    #dfF1 = pd.DataFrame(columns=column_namesF1)
    #dfTimes = pd.DataFrame(columns=column_namesTimes)

#    numpy.seterr(all='raise')
#    cfg = Config.Config(config,user)
    conn = connect_to_db(cfg.dbname, cfg.user, cfg.password, cfg.host, cfg.port)

    adom = [x[0] for x in execute_query(conn, "select distinct  " + cfg.sel + " from " + cfg.table + ";")]
    #table_size = execute_query(conn, "select count(1) from " + cfg.table + ";")[0][0]

    sizeOfR = getSizeOf(conn, cfg.table)

    tabR=[1,2,5,10]
    tabR=[1]
    tabSampleSize=[0.01,0.1,0.3,0.5,1]
    #tabSampleSize=[1]

    for coef in tqdm(tabR, desc='coef for r'):
        for percentSize in tqdm(tabSampleSize, desc='sample size'):
            samplesize=sizeOfR*percentSize


            start_time = time.time()
            ds=DataSampler(conn, cfg)
            if samplingMethod == 'naive':
                sample=ds.getSample(samplesize, adom, samplingMethod, adom)
                l = Lattice(sample,conn)
            else:
                newadom,congress = ds.getSample(samplesize, adom, samplingMethod)
                l = Lattice(congress,conn)

            if computeHyp:
                r=int(math.pow(2,len(cfg.groupbyAtt)-1))
                #to increase the chances of gaps in deltak, enabling drawing with replacement
                r=r*coef
                p=1
                #ranking=RankingFromPairwise(cfg.prefs, r,p)
                ranking=RankingFromPairwise(adom, r,p, 'Welch', True)
                ranking.run(l,method)
                #print('Delta:',ranking.delta)
                print('F:',ranking.F)
                #print('Tau:',ranking.tau)
                #print('M',ranking.M)
                end_time = time.time()
                timings = end_time - start_time
                print('Completed in ',timings, 'seconds')
                hypothesis=ranking.getHypothesis()
                print('Hypothesis:',hypothesis)
                tauHypothesis=ranking.getTauHypothesis()
                print('TauHypothesis:',tauHypothesis)
                print('OrderedN:',ranking.orderedN)
                l1, l2 = transform_to_rankings(hypothesis, tauHypothesis)
                print('Kendall tau between N and Tau: ', compute_kendall_tau(l1, l2))

            if computeGT:
                groupbyAtt = cfg.groupbyAtt[1:]
                sampleOfLattice=SampleRanking(conn, groupbyAtt, cfg.sel, cfg.meas, cfg.measBase, cfg.function, cfg.table, generateIndex=False)
                groundTruth=sampleOfLattice.getGTallLattice(adom,method)
                print('Ground truth: ',groundTruth)
                print('orderedN:',sampleOfLattice.orderedN)
                tauGT=sampleOfLattice.tauGT
                print('ordered tau:', tauGT)
                print('F:',sampleOfLattice.F)
                l1, l2 = transform_to_rankings(groundTruth, tauGT)
                print('Kendall tau between N and Tau: ', compute_kendall_tau(l1, l2))

            #for k in [3,5,10,16]:
            #l1,l2=transform_to_rankings(hypothesis[:k],groundTruth[:k])
            k=16
            l1,l2=transform_to_rankings(hypothesis,groundTruth)
            tau, pval = compute_kendall_tau(l1,l2)
            #print('Kendall tau between hypothesis and ground truth: ', compute_kendall_tau(l1,l2))
            dfError.loc[len(dfError)] = [r,samplesize,k,tau]
            print('kendall tau between hypothesis and ground truth:',tau)
    dfError.to_csv(fileResultsError, mode='a', header=True)

