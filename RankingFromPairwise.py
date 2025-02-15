import math
import time
import random

from scipy.sparse import dok_matrix
import numpy as np
from tqdm import tqdm

from utilities import sort_dict_descending
from dbStuff import getSample
from statStuff import welch_ttest, permutation_test, compute_skewness, benjamini_hochberg, benjamini_hochberg_statmod, claireStat
from Lattice import Lattice

# global variable
pairwiseComparison = []





class RankingFromPairwise:

    #values: tab of values to rank
    def __init__(self, values, r, p, test='Welch', replacement=True):
        self.values = values
        self.r=r
        self.p=p
        self.test=test
        self.replacement=replacement
        self.n=len(values)
        self.threshold=8 * math.sqrt(math.log( self.n )/ (self.n*self.p*self.r))
        print('threshold:', self.threshold)
        self.M = dok_matrix((self.n, self.n),dtype=np.float32)
        for i in range(self.n):
            self.M[i,i]=1/2
        self.tau=dict()
        self.delta=[]
        self.F=[]
        #initialize N
        self.N=dict()
        for i in range(len(self.values)):
            a=self.values[i]
            self.N[a]=0

    def getHypothesis(self):
        return self.hypothesis

    def getTauHypothesis(self):
        return self.tauHypothesis

    def get_rp(self):
        # r= nb of trials = nb of cuboids?
        # p= probability of each trial
        return self.r,self.p
        #return 64, 1

    #val is empirical probabilty that a beats b
    def updateM(self,a,b,val):
        i=self.values.index(a)
        j=self.values.index(b)
        self.M[i,j]=val

    def computeTau(self):
        self.tau=(1/self.n) * self.M.sum(axis=1)
        taudict={}
        for v in self.values:
            indexintau = self.values.index(v)
            tauv = self.tau[indexintau][0]
            taudict[v] = tauv
        orderedTau = sort_dict_descending(taudict)
        self.orderedTau = orderedTau
        orderedTautmp = sort_dict_descending(self.orderedTau)
        # print("ordered N", orderedN)
        self.tauHypothesis = list(orderedTautmp.keys())
        #print('Ordered Tau:',orderedTau)

    def binomialForPair(self,r,p):
        #return np.random.binomial(r,p)
        return r

    #compares a, b on in-memory cuboids
    def performComparisons(self, L, nb, a,b,method='withTest'):
        nbWon=0
        nbLost=0
        nbZeros=0
        nbFailedTest=0
        setOfCuboidsOnSample=L.getCuboids().copy()
        remaining=len(setOfCuboidsOnSample)-1

        match method:
            case 'withTest':
                for i in range(nb):
                    nbr = random.randint(0, remaining)
                    gb = setOfCuboidsOnSample[nbr]
                    if not self.replacement:
                        remaining = remaining - 1
                        setOfCuboidsOnSample.remove(gb)
                    res=L.compare(a,b,gb,self.test,method)
                    if res==1:
                        nbWon=nbWon+1
                    if res==0:
                        nbZeros=nbZeros+1
                    if res==-1:
                        nbLost=nbLost+1
                    if res==-2:
                        nbFailedTest=nbFailedTest+1
                self.N[a] = self.N[a] + nbWon
                self.N[b] = self.N[b] + nbLost
                if nbZeros==nb or nbZeros+nbFailedTest==nb:
                    self.updateM(a, b, .5)
                    self.updateM(b, a, .5)
                else:
                    self.updateM(a,b,nbWon/(nb-(nbZeros+nbFailedTest)))
                    self.updateM(b, a, 1- (nbWon / (nb - (nbZeros+nbFailedTest))))
            case 'withoutTest':
                for i in range(nb):
                    nbr = random.randint(0, remaining)
                    gb = setOfCuboidsOnSample[nbr]
                    if not self.replacement:
                        remaining = remaining - 1
                        setOfCuboidsOnSample.remove(gb)
                    #nbA,pA,nbB,pB=L.compare(a,b,gb,self.test)
                    nbA, pA, nbB, pB = L.compareWithoutTest(a, b, gb)
                    nbWon = nbWon + nbA
                    nbLost = nbLost + nbB
                #CHECK HERE
                self.N[a] = self.N[a] + nbWon
                self.N[b] = self.N[b] + nbLost
                self.updateM(a, b, nbWon / (nb))
                self.updateM(b, a, 1 - (nbWon / (nb)))
            case 'onlyFacts':
                res = L.compare(a, b, 'None', self.test, method)
                if res == 1:
                    nbWon = nbWon + 1
                if res == 0:
                    nbZeros = nbZeros + 1
                if res == -1:
                    nbLost = nbLost + 1
                if res == -2:
                    nbFailedTest = nbFailedTest + 1
                self.N[a] = self.N[a] + nbWon
                self.N[b] = self.N[b] + nbLost
                if nbZeros == nb or nbZeros + nbFailedTest == nb:
                    self.updateM(a, b, .5)
                    self.updateM(b, a, .5)
                else:
                    self.updateM(a, b, nbWon / (nb - (nbZeros + nbFailedTest)))
                    self.updateM(b, a, 1 - (nbWon / (nb - (nbZeros + nbFailedTest))))


    def run(self,L,method='withTest'):
        #for each pair in self.values:

        for i in range(len(self.values)):
            for j in range(i+1,len(self.values)):
                a,b=self.values[i],self.values[j]
                # get r and p
                r,p=self.get_rp()
                # draw number of comparison to make for each pair
                nbOfComp=self.binomialForPair(r,p)
                # make comparison and update M, N
                self.performComparisons(L,nbOfComp,a,b,method)
        # compute tau
        self.computeTau()
        # compute deltak
        self.computeDeltak()
        #check lowest k for which theorem 1 holds

    def computeDeltak(self):
        #order N desc
        orderedN=sort_dict_descending(self.N)
        self.orderedN=orderedN
        #print("ordered N", orderedN)
        self.hypothesis=list(orderedN.keys())
        orderedN=self.orderedTau
        #print("ordered keys", orderedN.keys())
        for k in range(len(orderedN.keys())-1):
            valk=list(orderedN.keys())[k]
            valNext=list(orderedN.keys())[k+1]
            indexkintau=self.values.index(valk)
            indexnextintaux=self.values.index(valNext)
            tauk=self.tau[indexkintau][0]
            tauNext = self.tau[indexnextintaux][0]
            self.delta.append(tauk-tauNext)
        self.checkGoodk()

    def checkGoodk(self):
        #r, p = self.get_rp()
        for k in range(len(self.delta)):
            if self.delta[k] >= self.threshold:
                self.F.append(True)
            else:
                self.F.append(False)


if __name__ == '__main__':
    test=RankingFromPairwise(['a', 'b', 'c', 'd', 'e', 'f'])

    test.run()

    print(test.M)
    print(test.tau)
    print(test.N)
    print(test.delta)
    print(test.F)




def computeRanksForAll(pairwiseComparison, Sels):
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

    #print("ranks:", ranks)
    return ranks


def balanced_rank_estimation(pairwise_comparisons, max_iterations=1000, tol=1e-6):
    """
    Perform Balanced Rank Estimation based on pairwise comparisons.

    Parameters:
    - pairwise_comparisons: list of tuples (item1, item2, result)
      where result is 1 if item1 > item2, -1 if item1 < item2, and 0 if item1 == item2.
    - max_iterations: maximum number of iterations for the algorithm (default is 1000).
    - tol: tolerance for convergence (default is 1e-6).

    Returns:
    - ranks: dictionary with items as keys and their estimated ranks as values.
    """
    # Extract unique items from the pairwise comparisons
    items = set()
    for item1, item2, result in pairwise_comparisons:
        items.add(item1)
        items.add(item2)

    items = list(items)
    n = len(items)
    item_index = {item: i for i, item in enumerate(items)}

    # Initialize ranks to zero
    ranks = np.zeros(n)

    # Create a matrix to keep track of comparisons and results
    comparison_matrix = np.zeros((n, n))
    for item1, item2, result in pairwise_comparisons:
        i, j = item_index[item1], item_index[item2]
        comparison_matrix[i, j] = result
        comparison_matrix[j, i] = -result

    for iteration in range(max_iterations):
        old_ranks = ranks.copy()

        for i in range(n):
            sum_comparisons = 0
            sum_ranks = 0
            for j in range(n):
                if comparison_matrix[i, j] != 0:
                    sum_comparisons += 1
                    sum_ranks += old_ranks[j] + comparison_matrix[i, j]

            if sum_comparisons > 0:
                ranks[i] = sum_ranks / sum_comparisons

        if np.linalg.norm(ranks - old_ranks, ord=1) < tol:
            break

    return {items[i]: ranks[i] for i in range(n)}




def findTuple(a, b, claireTab):
    for c in claireTab:
        if (c[0] == a and c[1] == b) or (c[1] == a and c[0] == b):
            return c


def compare(a, b):
    cl = findTuple(a, b)

    if tuple[2]:
        # print("Welch test can be used")
        #nbWelch = nbWelch + 1
        t_stat, p_value, conclusion = welch_ttest(cl[3], cl[4])

        comp = 0  # not significant
        if p_value < 0.05 and t_stat < 0:
            comp = -1
        if p_value < 0.05 and t_stat > 0:
            comp = 1
        pairwiseComparison.append((cl[0], cl[1], comp, t_stat, float(p_value)))
    else:
        # print("Permutation test is used")
        #nbPermut = nbPermut + 1
        observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(cl[3], cl[4])

        comp = 0  # not significant
        if p_value < 0.05 and observed_t_stat < 0:
            comp = -1
        if p_value < 0.05 and observed_t_stat > 0:
            comp = 1
        pairwiseComparison.append((cl[0], cl[1], comp, observed_t_stat, float(p_value)))
    result = False
    if comp == -1 or comp == 0:
        result = True
    return result, pairwiseComparison


def compare_cost(a, b, claireTab):
    """
    This function should be implemented to return the cost of comparing two objects a and b.
    """
    # Placeholder implementation; replace with the actual cost function.
    result = 0
    for c in claireTab:
        if (c[0] == a and c[1] == b) or (c[1] == a and c[0] == b):
            if c[2]:
                result = 0
            else:
                result = 1
    return result


def split_with_min_cost(arr, claireTab):
    """
    Splits the array into two parts with the minimum comparison cost.
    """
    n = len(arr)
    min_cost = float('inf')
    best_split = None

    for i in range(1, n):
        left = arr[:i]
        right = arr[i:]
        cost = sum(compare_cost(x, y, claireTab) for x in left for y in right)

        if cost < min_cost:
            min_cost = cost
            best_split = (left, right)

    return best_split


def merge(left, right, claireTab):
    """
    Merge two sorted lists into one sorted list.
    """
    sorted_list = []
    i = j = 0

    while i < len(left) and j < len(right):
        bool = compare(left[i], right[j], claireTab)  # True if left[i] <= right[j]
        if bool:  #left[i] <= right[j]
            sorted_list.append(left[i])
            i += 1
        else:
            sorted_list.append(right[j])
            j += 1

    sorted_list.extend(left[i:])
    sorted_list.extend(right[j:])

    return sorted_list


def merge_sort(arr, claireTab):
    """
    Modified merge sort that splits the array to minimize the cost of comparisons.
    """
    if len(arr) <= 1:
        return arr

    left, right = split_with_min_cost(arr, claireTab)
    sorted_left = merge_sort(left, claireTab)
    sorted_right = merge_sort(right, claireTab)

    return merge(sorted_left, sorted_right, claireTab)

# Example usage
#arr = [10, 20, 15, 5, 25, 30]
#sorted_arr = merge_sort(arr)
#print("Sorted array:", sorted_arr)

# computes separation as in JMLR 2018 By Shah and Wainwright
def computeSeparationJMLR18(pairwiseComparison,nbItems):
    dict={}
    for (c1,c2,comp,stat,proba) in pairwiseComparison:
        dict[c1]=[0.5]
        dict[c2]=[0.5]
    for (c1,c2,comp,stat,proba) in pairwiseComparison:
        if comp==1:
            list=dict[c1]
            list.append(1)
            dict[c1]=list
        if comp==-1:
            list=dict[c2]
            list.append(1)
            dict[c2]=list
        if comp==0:
            list1=dict[c1]
            list1.append(0.5)
            dict[c1]=list1
            list2=dict[c2]
            list2.append(0.5)
            dict[c2]=list2
    print(dict)
    separation={}
    for k in dict.keys():
        list=dict[k]
        tau=0
        if len(list)!=0:
            tau=sum(list)/nbItems
        separation[k]=tau
    print(separation)
    n=nbItems
    p=1
    r=1
    alpha=8
    probaError=1/(math.pow(n,14))
    bound=alpha * (math.sqrt( (math.log(n)) / (n*p*r) ) )
    print("bound:",bound)
    result=[]
    for s1 in separation.keys():
        scorek=separation[s1]
        for s2 in separation.keys():
            scorekplus1=separation[s2]
            if scorek-scorekplus1 >= bound:
                result.append((s1,s2))
    print("result:",result)
    return result

def generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, method, valsToSelect=None):
    # sampling
    start_time = time.time()
    resultVals = getSample(conn, measBase, table, sel, sampleSize, method, False, valsToSelect)
    end_time = time.time()
    samplingTime = end_time - start_time
    #print('sampling time:', samplingTime)

    #resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=DEBUG_FLAG)
    #print(resultVals)

    start_time = time.time()

    # get adom values
    Sels = list(set([x[0] for x in resultVals]))
    #print('Sels:',Sels)
    #analyse sample for each adom value: value, nb of measures, skewness, and tuples
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

    #print('S:',S)

    # nlog(n) comparisons enough for recovering the true ranking when comparisons are certain (not noisy)
    # we should try less
    #print(len(Sels))
    #nbOfComparisons = len(Sels) * math.log(len(Sels), 2)
    #print("Number of comparisons to make: " + str(nbOfComparisons))

    pairwiseComparison=generateAllComparisons(Sels, S)

    #separation=computeSeparationJMLR18(pairwiseComparison,len(valsToSelect))

    #print("pairwise comparisons:")
    #for p in pairwiseComparison:
    #    print("p: ", p)

    #pairwiseComparison = generateComparisonsWithMergeSort(Sels, S)

    # ranking
    #ranks = balanced_rank_estimation(pairwiseComparison)
    #print("Balanced Rank Estimation:", ranks)
    ranks = computeRanksForAll(pairwiseComparison, Sels)

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
    #print('Hypothesis generation time:', hypothesisGenerationTime)
    return hypothesis, samplingTime, hypothesisGenerationTime


def generateAllComparisons(Sels, S):
    #tabPValues = []
    pairwiseComparison = []
    #tabStat = []

    # compute Claire statistics for all pairs
    claireTab = []
    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i - 1][2], S[j][2], S[i - 1][1], S[j][1])
            claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))

            if b:
                # print("Welch test can be used")
                t_stat, p_value, conclusion = welch_ttest(S[i - 1][3], S[j][3])
                # print(t_stat, p_value, conclusion)
                #tabStat.append(t_stat)
                #tabPValues.append(float(p_value))
                comp = 0  # not significant
                if p_value < 0.05 and t_stat < 0:
                    comp = -1
                if p_value < 0.05 and t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[i - 1][0], S[j][0], comp, t_stat, float(p_value)))
            else:
                # print("Permutation test is used")
                observed_t_stat, p_value, permuted_t_stats, conclusion = permutation_test(S[i - 1][3], S[j][3])
                # print(f"Observed Welch's t-statistic: {observed_t_stat}")
                # print(f"P-value: {p_value}")
                # print(f"conclusion: {conclusion}")
                # print(observed_t_stat, p_value, conclusion)
                #tabStat.append(observed_t_stat)
                #tabPValues.append(float(p_value))
                comp = 0  # not significant
                if p_value < 0.05 and observed_t_stat < 0:
                    comp = -1
                if p_value < 0.05 and observed_t_stat > 0:
                    comp = 1
                pairwiseComparison.append((S[i - 1][0], S[j][0], comp, observed_t_stat, float(p_value)))

    pairwiseComparison = computeBHcorrection(pairwiseComparison, 0.05)
    #print('paiwise after BH:', pairwiseComparison)
    return pairwiseComparison


def computeBHcorrection(pairwiseComparison, alpha=0.05):
    # rejected, corrected_p_values = benjamini_hochberg_gpt(tabPValues, alpha)
    # print("Rejected hypotheses:", rejected)
    # print("raw p-values:", tabPValues)
    # print("Corrected p-values (gpt):", corrected_p_values)
    tabPValues = []
    for p in pairwiseComparison:
        tabPValues.append(p[4])

    corrected = benjamini_hochberg(tabPValues, alpha)
    rejected, corrected2 = benjamini_hochberg_statmod(tabPValues, alpha)

    #print("nb of True in rejected: ", utilities.nbTrueInList(rejected))

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

    #print("Number of BH corrections: ", nbChanges)

    #print("nb non zeros after corrections: ", utilities.countNonZeros(pairwiseComp2))

    return pairwiseComp2