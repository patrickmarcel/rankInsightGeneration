import dbStuff
from dbStuff import  getCuboidsOfAtt
from Config import Config
from pandas import DataFrame
import numpy as np
from statStuff import claireStat, welch_ttest, permutation_test, compute_skewness
import pandas as pd


class Lattice:
    """
    A simple Singleton class.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Lattice, cls).__new__(cls)
        return cls._instance

    def __init__(self, sample,conn):
        #Fetch db info
        cfg = Config()

        self.colNames = []
        self.colNames.extend(cfg.groupbyAtt)
        self.colNames.append(cfg.measBase)

        self.function=cfg.function

        self.measure = cfg.measBase
        self.selection = cfg.groupbyAtt[0]

        groupbyAtt = cfg.groupbyAtt[1:]
        pwset2 = getCuboidsOfAtt(groupbyAtt, cfg.sel)
        # remove last
        #self.cuboids= pwset2[:-1]
        self.cuboids = pwset2
        self.data = DataFrame(sample, columns=self.colNames)

        self.conn=conn
        #print(self.data)

    def getCuboids(self):
        return self.cuboids

    def agg(self, gb_set, function):
        return self.data.groupby([self.selection] + gb_set).agg({self.measure: function})

    def gb_agg(self, gb_set, function):
        if function=='avg':
            function='mean'
        return self.data.groupby(gb_set, as_index=False).agg({self.measure: function})

    def pairwise(self, gb_set, val1, val2, function):
        return self.data[(self.data[self.selection] == val1) | (self.data[self.selection] == val2)].groupby([self.selection] + gb_set).agg({self.measure: function})

    def getVal(self,val1,meas):
        return self.data[(self.data[self.selection] == val1)].loc[:,meas]

    def getValInGb(self,val1, val2, meas,gb):
        #cuboid = self.data.groupby(list(gb), group_keys=True).agg({self.measure: self.function})
        cuboid=self.gb_agg(list(gb),self.function)
        cuboid1 = cuboid[(cuboid[self.selection] == val1)]
        cuboid2 = cuboid[(cuboid[self.selection] == val2)]
        cuboid1=cuboid1.rename(columns={meas:val1})
        cuboid2=cuboid2.rename(columns={meas:val2})
        bgtmp=list(gb)
        bgtmp.remove(self.selection)
        bgtmp.append(val1)
        cuboid1=cuboid1[bgtmp]
        bgtmp.remove(val1)
        bgtmp.append(val2)
        cuboid2=cuboid2[bgtmp]
        #join=pd.concat([cuboid1, cuboid2], axis=1,  join="inner")
        bgtmp.remove(val2)
        if len(gb)>1:
            join=cuboid1.merge(cuboid2, how='inner', on=bgtmp)
            #aonmv, bonmv, query = self.checkOnMVs(val1, val2, gb)
            valA=join.loc[:,val1]
            valB=join.loc[:,val2]
            #if len(set(valA).difference(set(aonmv))) != 0 or len(set(valB).difference(set(bonmv))) != 0:
            #    print(query)
            #    print(len(set(valA)))
            #    print(len(set(aonmv)))
            return valA, valB
        else:
            valA = cuboid1.loc[:,val1]
            valB = cuboid2.loc[:, val2]
            return valA, valB

    def checkOnMVs(self, a, b, gb):
        #groupby=gb[0].split(',')[:-1]
        groupby=list(gb)
        gb2=list(gb)[:-1]
        strgb=''
        for s in groupby:
            strgb=strgb+s+','
        strgb=strgb[:-1]
        strgb2 = ''
        for s in gb2:
            strgb2 = strgb2 + s + ','
        strgb2 = strgb2[:-1]
        if strgb2=='':
            queryGetGroupByVals = ('select \"'+ a + '\",\"' + b +'\" from (select ' + self.measure + ' as \"' + a + '\" from \"' + strgb + '\" where ' + self.selection + '= \'' + a + '\' ) natural join (' +
                             'select '  + self.measure + ' as \"' + b + '\" from \"' + strgb + '\" where ' + self.selection + '= \'' + b + '\' );')
        else:
            queryGetGroupByVals=('select \"'+ a + '\",\"' + b +'\" from (select ' + strgb2  + ', ' + self.measure + ' as \"' + a + '\" from \"' + strgb + '\" where ' + self.selection + '= \'' + a + '\' ) natural join (' +
                             'select ' + strgb2  + ', ' + self.measure + ' as \"' + b + '\" from \"' + strgb + '\" where ' + self.selection + '= \'' + b + '\' );')

        res=dbStuff.execute_query(self.conn, queryGetGroupByVals)
        return [t[0] for t in res],[t[1] for t in res],queryGetGroupByVals


    # return 0 if no comparison, 1 if a>b, -1 if b>a using statistical tests
    def compare(self, a, b, gb, test ='stat',method='withTest'):
        S = []
        if method=='withTest':
            valsA=np.array(self.getVal(a,self.measure))
            valsB=np.array(self.getVal(b,self.measure))
            #valsA, valsB = np.array(self.getValInGb(a, b, self.measure, gb))
            nA = len(valsA)
            nB = len(valsB)
            skewA = compute_skewness(valsA)
            skewB = compute_skewness(valsB)
            S.append((a, nA, skewA, valsA))
            S.append((b, nB, skewB, valsB))
            return self.runStatisticalTest(S, test)
        else:
            #valsA = np.array(self.getValInGb(a, self.measure, gb))
            #valsB = np.array(self.getValInGb(b, self.measure, gb))
            return self.compareWithoutTest(a,b,gb)

    # compare without test returning probabilities of a wins over b
    def compareWithoutTest(self, a, b, gb):
        S = []
        #valsA=np.array(self.getVal(a,self.measure))
        #valsB=np.array(self.getVal(b,self.measure))
        valsA,valsB = np.array(self.getValInGb(a, b, self.measure,gb))
        #print("group by:",gb)
        #print('a',a)
        #print('b',b)
        #print("valsA:",valsA)
        #print("valsB:",valsB)
        #valsB = np.array(self.getValInGb(b, self.measure,gb))
        nA = len(valsA)
        nB = len(valsB)
        #skewA = compute_skewness(valsA)
        #skewB = compute_skewness(valsB)
        #S.append((a, nA, skewA, valsA))
        #S.append((b, nB, skewB, valsB))
        #seriesA = S[0][3]
        #seriesB = S[1][3]
        seriesA=valsA
        seriesB=valsB
        nbWonA = 0
        nbWonB = 0
        for i in range(len(seriesA)):
            if seriesA[i] > seriesB[i]:
                nbWonA = nbWonA + 1
            if seriesA[i] < seriesB[i]:
                nbWonB = nbWonB + 1
        if (nbWonA+nbWonB)!=0:
            probaWonA = nbWonA / (nbWonA+nbWonB)
            probaWonB = nbWonB / (nbWonA+nbWonB)
            return nbWonA, probaWonA, nbWonB, probaWonB
        else:
            return 0,0,0,0


    # legacy code
    def runStatisticalTest(self, S, test='stat'):
        b = claireStat(S[0][2], S[1][2], S[0][1], S[1][1])
        pairwiseComparison = []
        #claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))
        match test:
            case 'stat':
                if b:
                    # print("Welch test can be used")
                    #self.nbWelch = self.nbWelch + 1
                    comp = 0  # not significant
                    if len(S[0][3]) <= 1 or len(S[1][3]) <= 1:
                        comp = -2
                    else:
                        t_stat, p_value, conclusion = welch_ttest(S[0][3], S[1][3])
                        if p_value < 0.05 and t_stat < 0:
                            comp = -1
                        if p_value < 0.05 and t_stat > 0:
                            comp = 1
                    #pairwiseComparison.append((S[0][0], S[1][0], comp, t_stat, float(p_value)))
                else:
                    # print("Permutation test is used")
                    #self.nbPerm = self.nbPerm + 1
                    comp = 0  # not significant
                    if len(S[0][3] <= 1) or len(S[1][3] <= 1):
                        comp = -2
                    else:
                        observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[0][3], S[1][3])
                        if p_value < 0.05 and observed_t_stat < 0:
                            comp = -1
                        if p_value < 0.05 and observed_t_stat > 0:
                            comp = 1
                    #pairwiseComparison.append((S[0][0], S[1][0], comp, observed_t_stat, float(p_value)))
            case 'Welch':
                # only Welch
                #self.nbWelch = self.nbWelch + 1
                comp = 0  # not significant
                if len(S[0][3])<=1 or len(S[1][3])<=1:
                    comp=-2
                else:
                    t_stat, p_value, conclusion = welch_ttest(S[0][3], S[1][3])
                #if len(S[0][3])==0 or len(S[1][3])==0:
                #    print("empty sample !!")
                    if p_value < 0.05 and t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and t_stat > 0:
                        comp = 1
#                pairwiseComparison.append((S[0][0], S[1][0], comp, t_stat, float(p_value)))
            case 'Permutation':
                # only permutation
                #self.nbPerm = self.nbPerm + 1
                comp = 0  # not significant
                if len(S[0][3]) <= 1 or len(S[1][3]) <= 1:
                    comp = -2
                else:
                    observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[0][3], S[1][3])
                    if p_value < 0.05 and observed_t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and observed_t_stat > 0:
                        comp = 1
                #pairwiseComparison.append((S[0][0], S[1][0], comp, observed_t_stat, float(p_value)))
        return comp




def test():
    cfg = Config()
    # Enumerate all cuboids to construct
    sel = cfg.groupbyAtt[0]
    groupbyAtt = cfg.groupbyAtt[1:]
    pwset2 = getCuboidsOfAtt(groupbyAtt, sel)
    # remove last
    pwset2 = pwset2[:-1]

    for i in range(len(pwset2)):
        gb = pwset2[i]
        gbs = ''
        for s in gb:
            gbs = gbs + s + ','
        gbs = gbs[:-1]

        print(gb, cfg.function + "(" + cfg.meas + ")")