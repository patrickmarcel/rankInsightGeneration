from dbStuff import  getCuboidsOfAtt
from Config import Config
from pandas import DataFrame
import numpy as np
from statStuff import claireStat, welch_ttest, permutation_test, compute_skewness


class Lattice:
    """
    A simple Singleton class.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Lattice, cls).__new__(cls)
        return cls._instance

    def __init__(self, sample):
        #Fetch db info
        cfg = Config()

        self.colNames = []
        self.colNames.extend(cfg.groupbyAtt)
        self.colNames.append(cfg.measBase)

        self.measure = cfg.measBase
        self.selection = cfg.groupbyAtt[0]

        self.data = DataFrame(sample, columns=self.colNames)

        print(self.data)

    def agg(self, gb_set, function):
        return self.data.groupby([self.selection] + gb_set).agg({self.measure: function})

    def pairwise(self, gb_set, val1, val2, function):
        return self.data[(self.data[self.selection] == val1) | (self.data[self.selection] == val2)].groupby([self.selection] + gb_set).agg({self.measure: function})

    def getVal(self,val1,meas):
        return self.data[(self.data[self.selection] == val1)].loc[:,meas]

    # return 0 if no comparison, 1 if a>b, -1 if b>a
    def compare(self, a, b, test ='stat'):

        S = []
        valsA=np.array(self.getVal(a,self.measure))
        valsB=np.array(self.getVal(b,self.measure))
        nA = len(valsA)
        nB = len(valsB)
        skewA = compute_skewness(valsA)
        skewB = compute_skewness(valsB)
        S.append((a, nA, skewA, valsA))
        S.append((b, nB, skewB, valsB))
        return self.comparison(S,test)

    def comparison(self,S,test='stat'):
        b = claireStat(S[0][2], S[1][2], S[0][1], S[1][1])
        pairwiseComparison = []
        #claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))
        match test:
            case 'stat':
                if b:
                    # print("Welch test can be used")
                    #self.nbWelch = self.nbWelch + 1
                    t_stat, p_value, conclusion = welch_ttest(S[0][3], S[1][3])

                    comp = 0  # not significant
                    if p_value < 0.05 and t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and t_stat > 0:
                        comp = 1
                    pairwiseComparison.append((S[0][0], S[1][0], comp, t_stat, float(p_value)))
                else:
                    # print("Permutation test is used")
                    #self.nbPerm = self.nbPerm + 1
                    observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[0][3], S[1][3])

                    comp = 0  # not significant
                    if p_value < 0.05 and observed_t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and observed_t_stat > 0:
                        comp = 1
                    pairwiseComparison.append((S[0][0], S[1][0], comp, observed_t_stat, float(p_value)))
            case 'Welch':
                # only Welch
                #self.nbWelch = self.nbWelch + 1
                t_stat, p_value, conclusion = welch_ttest(S[0][3], S[1][3])

                comp = 0  # not significant
                if p_value < 0.05 and t_stat < 0:
                    comp = -1
                if p_value < 0.05 and t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[0][0], S[1][0], comp, t_stat, float(p_value)))
            case 'Permutation':
                # only permutation
                #self.nbPerm = self.nbPerm + 1
                observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[0][3], S[1][3])

                comp = 0  # not significant
                if p_value < 0.05 and observed_t_stat < 0:
                    comp = -1
                if p_value < 0.05 and observed_t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[0][0], S[1][0], comp, observed_t_stat, float(p_value)))
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