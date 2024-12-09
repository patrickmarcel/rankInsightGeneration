import time
import numpy as np
from dbStuff import getSample, execute_query
from dolap import countViolations
from statStuff import welch_ttest, permutation_test, compute_skewness, benjamini_hochberg, claireStat
from statsmodels.stats.multitest import fdrcorrection


class Hypothesis:
    def __init__(self ):
        self.nbWelch=0
        self.nbPerm=0

    def getNbWelch(self):
            return self.nbWelch
    def getNbPerm(self):
            return self.nbPerm

    def computeRanksForAll(self, pairwiseComparison, Sels):
        """
        For each item in adom will count how many times it has 'won' a pairwise comparison
        :param pairwiseComparison:
        :param Sels:
        :return:
        """
        # this is borda count style
        ranks = {}
        for s in Sels:
            ranks[s] = 0
        for left, right, result, _, _ in pairwiseComparison:
            if result == 1:
                ranks.update({left: ranks[left] + 1})
            if result == -1:
                ranks.update({right: ranks[right] + 1})

        # print("ranks:", ranks)
        return ranks

    def computeBHcorrection(self, pairwiseComparison, alpha=0.05):
        tabPValues = []
        for p in pairwiseComparison:
            tabPValues.append(p[4])

        corrected = benjamini_hochberg(tabPValues, alpha)

        pairwiseComp2 = []
        i = 0
        nbChanges = 0
        for c in pairwiseComparison:
            comp = 0  # not significant
            if corrected[i] < 0.05 and c[3] < 0:
                comp = -1
            if corrected[i] < 0.05 and c[3] > 0:
                comp = 1
            if comp != c[2]:
                nbChanges = nbChanges + 1
            pairwiseComp2.append((c[0], c[1], comp, c[2], corrected[i]))
            i = i + 1

        return pairwiseComp2

    def generateAllComparisons(self, Sels, S):
        # tabPValues = []
        pairwiseComparison = []
        # tabStat = []

        # compute Claire statistics for all pairs
        claireTab = []
        for i in range(1, len(S)):
            for j in range(i, len(S)):
                b = claireStat(S[i - 1][2], S[j][2], S[i - 1][1], S[j][1])
                claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))

                if b:
                    # print("Welch test can be used")
                    self.nbWelch=self.nbWelch + 1
                    t_stat, p_value, conclusion = welch_ttest(S[i - 1][3], S[j][3])

                    comp = 0  # not significant
                    if p_value < 0.05 and t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and t_stat > 0:
                        comp = 1
                    pairwiseComparison.append((S[i - 1][0], S[j][0], comp, t_stat, float(p_value)))
                else:
                    # print("Permutation test is used")
                    self.nbPerm=self.nbPerm + 1
                    observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[i - 1][3], S[j][3])

                    comp = 0  # not significant
                    if p_value < 0.05 and observed_t_stat < 0:
                        comp = -1
                    if p_value < 0.05 and observed_t_stat > 0:
                        comp = 1
                    pairwiseComparison.append((S[i - 1][0], S[j][0], comp, observed_t_stat, float(p_value)))

        pairwiseComparison = self.computeBHcorrection(pairwiseComparison, 0.05)
        # print('paiwise after BH:', pairwiseComparison)
        return pairwiseComparison

    def generateHypothesisTest(self, conn, meas, measBase, table, sel, sampleSize, method, valsToSelect=None):
        # sampling
        start_time = time.time()
        resultVals = getSample(conn, measBase, table, sel, sampleSize, method, False, valsToSelect)
        end_time = time.time()
        samplingTime = end_time - start_time


        start_time = time.time()

        # get adom values
        Sels = list(set([x[0] for x in resultVals]))
        # print('Sels:',Sels)
        # analyse sample for each adom value: value, nb of measures, skewness, and tuples
        S = []
        for v in Sels:

            data = []
            for row in resultVals:
                if row[0] == v:
                    data.append(float(row[1]))

            nvalues = len(data)
            data = np.array(data)
            skewness = compute_skewness(data)
            S.append((v, nvalues, skewness, data))


        pairwiseComparison = self.generateAllComparisons(Sels, S)

        ranks = self.computeRanksForAll(pairwiseComparison, Sels)

        sorted_items = sorted(ranks.items(), key=lambda item: item[1], reverse=True)

        # Construct a rank from the number of comparison won for each adom values
        hypothesis = []
        rank = 0
        for s in sorted_items:
            if rank == 0:
                rank = 1
                hypothesis.append((s[0], rank))
                val = s[1]
            else:
                if s[1] == val:
                    hypothesis.append((s[0], rank))
                    val = s[1]
                else:
                    rank = rank + 1
                    hypothesis.append((s[0], rank))
                    val = s[1]

        end_time = time.time()
        hypothesisGenerationTime = end_time - start_time
        # print('Hypothesis generation time:', hypothesisGenerationTime)
        return hypothesis, samplingTime, hypothesisGenerationTime

    def hypothesisGeneration(self,conn, prefs, sel, measBase, meas, table, sampleSize, allComparison):
        if allComparison == False:
            # sampling
            start_time = time.time()
            adom, congress = self.fetchCongressionalSample(conn, sel, table, measBase, sampleSize, adom_restr=prefs)
            end_time = time.time()
            samplingTime = end_time - start_time

            # compute hypothesis
            start_time = time.time()
            hypothesis = self.getHypothesisCongressionalSampling(adom, congress)
            end_time = time.time()
            hypothesisGenerationTime = end_time - start_time
        else:
            # sampling and hypothesis
            hypothesis, samplingTime, hypothesisGenerationTime, pvalue = self.getHypothesisAllComparisons(conn, meas,
                                                                                                     measBase, table,
                                                                                                     sel, tuple(prefs),
                                                                                                     sampleSize,
                                                                                                     method='SYSTEM_ROWS')
        return hypothesis, hypothesisGenerationTime, samplingTime, pvalue

    def getHypothesisAllComparisons(self, conn, meas, measBase, table, sel, valsToSelect, sampleSize, method='SYSTEM_ROWS'):
        # checking all comparisons
        correctHyp, samplingTime, hypothesisGenerationTime, pvalue = self.generateHypothesisTestDolap(conn, meas, measBase,
                                                                                                 table, sel, sampleSize,
                                                                                                 method, valsToSelect)
        ##print('all comp. hypothesis:', correctHyp)
        return correctHyp, samplingTime, hypothesisGenerationTime, pvalue

    def get_state_sample(self, conn, measBase, table, sel, sampleSize, state):

        querySample = "SELECT " + sel + ", " + measBase + " FROM " + table + " where " + sel + " = '" + str(
            state) + "' limit " + str(sampleSize) + ";"
        ##print('stat query:', querySample)
        resultVals = execute_query(conn, querySample)
        return resultVals

    def fetchCongressionalSample(self, conn, sel, table, measBase, sampleSize, adom_restr=None):
        # fetch the congressional sample
        if adom_restr:
            adom = adom_restr
        else:
            adom = [x[0] for x in execute_query(conn, "select distinct  " + sel + " from " + table + ";")]
        table_size = execute_query(conn, "select count(1) from " + table + ";")[0][0]

        sample_size = int(table_size * sampleSize)
        alpha = 0.10
        house_size = sample_size * alpha
        senate_size = sample_size * (1 - alpha)

        house = getSample(conn, measBase, table, sel, house_size, method="SYSTEM_ROWS", repeatable=False)

        senate = []
        state_sample_size = int(senate_size / len(adom))
        for state in adom:
            senate.extend(self.get_state_sample(conn, measBase, table, sel, state_sample_size, state))

        if adom_restr:
            house = list(filter(lambda x: x[0] in adom_restr, house))
        congress = house + senate
        # END - fetch the congressional sample
        return adom, congress

    def getHypothesisCongressionalSampling(self, adom, congress):

        buckets = {s: [] for s in adom}
        skews = dict()
        for item in congress:
            buckets[item[0]].append(item[1])
        for k in buckets.keys():
            skews[k] = compute_skewness(buckets[k])

        # do all welch tests
        param_budget = 20
        param_budget = int(param_budget / 2)

        from scipy.stats import ttest_ind

        raw_comparisons = []

        for i in range(len(adom)):
            for j in range(i + 1, len(adom)):
                left = adom[i]
                right = adom[j]
                res = ttest_ind(buckets[left], buckets[right], equal_var=False)
                stat_c = claireStat(skews[left], skews[right], len(buckets[left]), len(buckets[right]))
                if res.statistic < 0:
                    raw_comparisons.append((left, right, stat_c, res.pvalue))
                else:
                    raw_comparisons.append((right, left, stat_c, res.pvalue))

        w_comparisons = []
        w_comparisons_rej = []
        ##print(raw_comparisons)
        rejected, corrected = fdrcorrection([x[3] for x in raw_comparisons], alpha=0.05)
        for i in range(len(raw_comparisons)):
            if rejected[i]:
                w_comparisons_rej.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))
            else:
                w_comparisons.append((raw_comparisons[i][0], raw_comparisons[i][1], raw_comparisons[i][2]))

        # print("NB de comparaisons significatives (welch)", len(w_comparisons))
        # #print_comp_list(sorted(w_comparisons, key=lambda x: x[0] + x[1]))
        by_prox_to_threshold = sorted(w_comparisons, key=lambda x: abs(0.05 - x[2]), reverse=True)
        # #print(by_prox_to_threshold)

        final = by_prox_to_threshold[param_budget:]
        to_redo = by_prox_to_threshold[:param_budget]

        to_redo.extend(sorted(w_comparisons_rej, key=lambda x: abs(0.05 - x[2]), reverse=True)[:param_budget])

        for left, right, _ in to_redo:
            res = permutation_test(buckets[left], buckets[right])
            if res[3].startswith("Reject"):
                if res[1] > 0:
                    final.append((left, right, -1))
                else:
                    final.append((right, left, -1))

        # print("NB de comparaisons significatives (welch + X param)", len(final))
        # #print_comp_list(sorted(final, key=lambda x: x[0] + x[1]))

        # borda hypothesis
        patrick_format = [(a, b, 1, None, None) for (a, b, c) in final]
        ##print('alex pairwise:',patrick_format)
        hypothesis = self.computeRanksForAll(patrick_format, adom).items()
        ##print(hypothesis)
        hypothesis = sorted(
            hypothesis,
            key=lambda x: x[1],
            reverse=True
        )
        ##print(hypothesis)
        # hypothesis = [(a, b + 1) for (a, b) in hypothesis]
        correctHyp = []
        i = 1
        prevb = -1
        for (a, b) in hypothesis:
            if prevb == -1:
                currentRank = 1
                correctHyp.append((a, currentRank))
                prevb = b
            else:
                if b == prevb:
                    correctHyp.append((a, currentRank))
                else:
                    currentRank = currentRank + 1
                    correctHyp.append((a, currentRank))
                    prevb = b

        ##print('Alex hypothesis:',correctHyp)
        return correctHyp

    def generateHypothesisTestDolap(self,conn, meas, measBase, table, sel, sampleSize, method, valsToSelect=None):
        # sampling
        start_time = time.time()
        resultVals = getSample(conn, measBase, table, sel, sampleSize, method, False, valsToSelect)
        end_time = time.time()
        samplingTime = end_time - start_time
        # print('sampling time:', samplingTime)

        # resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=DEBUG_FLAG)
        # print(resultVals)

        start_time = time.time()

        # get adom values
        Sels = list(set([x[0] for x in resultVals]))
        # print('Sels:',Sels)
        # analyse sample for each adom value: value, nb of measures, skewness, and tuples
        S = []
        for v in Sels:

            data = []
            for row in resultVals:
                if row[0] == v:
                    data.append(float(row[1]))

            nvalues = len(data)
            data = np.array(data)
            skewness = compute_skewness(data)
            S.append((v, nvalues, skewness, data))

        # print('S:',S)

        # nlog(n) comparisons enough for recovering the true ranking when comparisons are certain (not noisy)
        # we should try less
        # print(len(Sels))
        # nbOfComparisons = len(Sels) * math.log(len(Sels), 2)
        # print("Number of comparisons to make: " + str(nbOfComparisons))

        pairwiseComparison = self.generateAllComparisons(Sels, S)

        # separation=computeSeparationJMLR18(pairwiseComparison,len(valsToSelect))

        # print("pairwise comparisons:")
        # for p in pairwiseComparison:
        #    print("p: ", p)

        # pairwiseComparison = generateComparisonsWithMergeSort(Sels, S)

        # ranking
        # ranks = balanced_rank_estimation(pairwiseComparison)
        # print("Balanced Rank Estimation:", ranks)
        ranks = self.computeRanksForAll(pairwiseComparison, Sels)

        sorted_items = sorted(ranks.items(), key=lambda item: item[1], reverse=True)

        # Construct a rank from the number of comparison won for each adom values
        hypothesis = []
        rank = 0
        for s in sorted_items:
            if rank == 0:
                rank = 1
                hypothesis.append((s[0], rank))
                val = s[1]
            else:
                if s[1] == val:
                    hypothesis.append((s[0], rank))
                    val = s[1]
                else:
                    rank = rank + 1
                    hypothesis.append((s[0], rank))
                    val = s[1]

        end_time = time.time()
        hypothesisGenerationTime = end_time - start_time
        # print('Hypothesis generation time:', hypothesisGenerationTime)
        if pairwiseComparison != []:
            pvalue = float(pairwiseComparison[0][4])
        else:
            pvalue = 1000  # change me
        return hypothesis, samplingTime, hypothesisGenerationTime, pvalue

    def checkViolationsOpt(self, ratioViolations, ratioCuboidOK, conn, ranks, hypothesis, queryCountCuboid,sizeofsample):
        # cut here
        #  if ratio is reached, or if sure ratio won't be reached given nb of cuboids remaining, then stop don't test more cuboids
        nbViewOK = 0
        nbInconclusive = 0
        stop = False
        i = 0
        while not stop and i < len(ranks):
            v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
            c = execute_query(conn, queryCountCuboid[i])[0][0]
            if c != 0:
                if ratio < ratioViolations:
                    nbViewOK = nbViewOK + 1
            else:
                sizeofsample = sizeofsample - 1
                nbInconclusive = nbInconclusive + 1
            i = i + 1
            if sizeofsample != 0:
                prediction = nbViewOK / sizeofsample
                if prediction >= ratioCuboidOK or (nbViewOK + len(ranks) - i) / sizeofsample < ratioCuboidOK:
                    stop = True
        # cut here
        return nbViewOK,sizeofsample,nbInconclusive

    def checkViolations(self, ratioViolations, conn, ranks, hypothesis, queryCountCuboid, queryCountviolations, nbMVs):
        nbInconclusive = 0
        nbViewOK = 0
        for i in range(len(queryCountviolations)):

            v, ratio, qtime = countViolations(conn, ranks[i], hypothesis)
            c = execute_query(conn, queryCountCuboid[i])[0][0]

            if c != 0:
                if ratio < ratioViolations:
                    nbViewOK = nbViewOK + 1
            else:
                nbMVs = nbMVs - 1
                nbInconclusive = nbInconclusive + 1

        if nbMVs == 0:
            prediction = 0
            # bennetError = 0
        else:
            prediction = nbViewOK / nbMVs

        return prediction,nbViewOK,nbMVs,nbInconclusive