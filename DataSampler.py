from random import randint

from dbStuff import getSample, execute_query, getSample_new
from statStuff import welch_ttest, permutation_test, compute_skewness, benjamini_hochberg, claireStat
from statsmodels.stats.multitest import fdrcorrection


# so far two sampling methods available, naive and congressional
class DataSampler():

    def __init__(self, conn, cfg):
        self.conn = conn
        self.cfg=cfg

    def getSample(self, sampleSize, adom, samplingMethod = 'naive', method="SYSTEM_ROWS", repeatable=False):
        match samplingMethod:
            case 'naive':
                sample=self.getNaive(sampleSize, method, repeatable)
                return list(map(lambda x: x[1:], sample))
            case 'congressional': return self.getCongressional(sampleSize, adom)

    def getNaive(self, sampleSize, method="SYSTEM_ROWS", repeatable=False):
        #resultVals = getSample(self.conn, self.cfg.measBase, self.cfg.table, self.cfg.sel, sampleSize, method, repeatable, None)
        resultVals=getSample_new(self.conn, sampleSize)
        #print('size of sample:',len(resultVals))
        return resultVals

    def getCongressional(self, sampleSize, adom):
        #adom, congress = self.fetchCongressionalSample(self.conn, self.cfg.sel, self.cfg.table, self.cfg.measBase, sampleSize, adom_restr=self.cfg.prefs)
        adom, congress = self.fetchCongressionalSample(self.conn, self.cfg.sel, self.cfg.table, self.cfg.measBase, sampleSize, adom)
        return adom, congress

    def get_state_sample(self, conn, measBase, table, sel, sampleSize, state):

        querySample = "SELECT ctid, " + sel + ", " + measBase + " FROM " + table + " where " + sel + " = '" + str(state) + "' limit " + str(sampleSize) + ";"

        # modifs patrick
        gbs = ''
        for s in self.cfg.groupbyAtt:
            gbs = gbs + s + ','
        gbs = gbs[:-1]
        querySample = ("SELECT " + gbs + ", " + measBase + " FROM " + table + " TABLESAMPLE " + "SYSTEM_ROWS" + " (" + str(sampleSize) + ")" + " where " + sel + " = '" + str(state) + "';")

        #print('state query:', querySample)
        resultVals = execute_query(conn, querySample)
        return resultVals

    def fetchCongressionalSample(self, conn, sel, table, measBase, sampleSize, adom, adom_restr=None):
        # fetch the congressional sample
        if adom_restr:
            adom = adom_restr
        #else:
        #    adom = [x[0] for x in execute_query(conn, "select distinct  " + sel + " from " + table + ";")]
        #table_size = execute_query(conn, "select count(1) from " + table + ";")[0][0]

        #sample_size = int(table_size * sampleSize)
        sample_size = sampleSize
        alpha = 0.5
        house_size = sample_size * alpha
        senate_size = sample_size * (1 - alpha)

        #house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)
        house = getSample_new(conn, house_size + 0.1*senate_size, method="SYSTEM_ROWS", repeatable=False)


        senate = set()
        state_sample_size = int(senate_size / len(adom))
        for state in adom:
            senate.update(self.get_state_sample(conn, measBase, table, sel, state_sample_size, state))

        if adom_restr:
            house = list(filter(lambda x: x[0] in adom_restr, house))

        #remove duplicates
        house_clean = []
        for sample in house:
            if sample not in senate:
                house_clean.append(sample)

        while len(house_clean)+len(senate) > house_size + senate_size:
            house_clean.pop(randint(0, len(house_clean) - 1))
        congress = house_clean + list(senate)
        #print(congress)
        #print("sampler", len(house_clean), "/", len(house))

        # END - fetch the congressional sample
        #print("sampler", len(congress), "/", house_size+senate_size)

        #return adom, self.getSample(sampleSize, adom, samplingMethod = 'naive', method="SYSTEM_ROWS", repeatable=False)
        return adom, list(map(lambda x : x[1:],congress))

