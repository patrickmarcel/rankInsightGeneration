from dbStuff import  getCuboidsOfAtt
from Config import Config
from pandas import DataFrame

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