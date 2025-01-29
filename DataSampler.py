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
            case 'naive': return self.getNaive(sampleSize, method, repeatable)
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

        #querySample = "SELECT " + sel + ", " + measBase + " FROM " + table + " where " + sel + " = '" + str(state) + "' limit " + str(sampleSize) + ";"

        querySample = ("SELECT " + sel + ", " + measBase + " FROM " + table + " TABLESAMPLE " + "SYSTEM_ROWS" + " (" + str(sampleSize) + ")" + " where " + sel + " = '" + str(state) + "';")

        ##print('stat query:', querySample)
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
        alpha = 0.10
        house_size = sample_size * alpha
        senate_size = sample_size * (1 - alpha)

        #house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)
        house = getSample_new(conn, house_size, method="SYSTEM_ROWS", repeatable=False)


        senate = []
        state_sample_size = int(senate_size / len(adom))
        for state in adom:
            senate.extend(self.get_state_sample(conn, measBase, table, sel, state_sample_size, state))

        if adom_restr:
            house = list(filter(lambda x: x[0] in adom_restr, house))
        congress = house + senate
        # END - fetch the congressional sample
        return adom, congress

