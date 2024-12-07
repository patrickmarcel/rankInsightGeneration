import time
import numpy as np
from dbStuff import getSample
from statStuff import welch_ttest, permutation_test, compute_skewness, benjamini_hochberg, benjamini_hochberg_statmod, claireStat
from dolap import fetchCongressionalSample, getHypothesisCongressionalSampling, getHypothesisAllComparisons


class Hypothesis:
    def __init__(self ):
        self.nbWelch=0
        self.nbPerm=0

        def getNbWelch():
            return self.nbWelch
        def getNbPerm():
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
            adom, congress = fetchCongressionalSample(conn, sel, table, measBase, sampleSize, adom_restr=prefs)
            end_time = time.time()
            samplingTime = end_time - start_time

            # compute hypothesis
            start_time = time.time()
            hypothesis = getHypothesisCongressionalSampling(adom, congress)
            end_time = time.time()
            hypothesisGenerationTime = end_time - start_time
        else:
            # sampling and hypothesis
            hypothesis, samplingTime, hypothesisGenerationTime, pvalue = getHypothesisAllComparisons(conn, meas,
                                                                                                     measBase, table,
                                                                                                     sel, tuple(prefs),
                                                                                                     sampleSize,
                                                                                                     method='SYSTEM_ROWS')
        return hypothesis, hypothesisGenerationTime, samplingTime, pvalue


