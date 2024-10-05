def generateRandomQuery(pwsert, valsToSelect, hypothesis):
    nb = random.randint(0, len(pwsert) - 1)
    gb = pwsert[nb]
    strgb = ""
    for i in range(len(gb)):
        strgb = strgb + str(gb[i])
        if i != len(gb) - 1:
            strgb = strgb + ","

    print("group by is: " + strgb)

    # for debugging
    # strgb = "departure_airport"

    #print("vals in gen queries:", valsToSelect)
    hyp = ""
    for i in range(len(hypothesis)):
        hyp = hyp + str(hypothesis[i])
        if i != len(hypothesis) - 1:
            hyp = hyp + ","
    queryHyp = (
            "select * from (select  " + strgb + " from  " + table + ") t1  cross join (values " + hyp + ") as t2 ")

    query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
             + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
             " FROM " + table + " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ")

    queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "
                   + table + " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ) x;")

    queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + " ")

    queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
    queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")

    #return query, queryHyp, queryValues, queryExcept, strgb
    return queryValues, queryCountGb, queryCountExcept


def getValues(queryValues, vals, v, conn):
    queryVal = queryValues.replace(str(vals), "('" + v + "')")
    resultValues = execute_query(conn, queryVal)
    data = []
    for row in resultValues:
        data.append(float(row[0]))

    np.array(data)
    return data



def generateHypothesisTest(conn, meas, measBase, table, sel, sampleSize, method):
    resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=False)
    #resultVals = getSample(conn, measBase, table, sel, sampleSize, method=method, repeatable=DEBUG_FLAG)

    # get adom values
    Sels = list(set([x[0] for x in resultVals]))

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

    #print(S)

    # nlog(n) comparisons enough for recovering the true ranking when commarisons are certain (not noisy)
    # we should try less
    nbOfComparisons = len(Sels) * math.log(len(Sels), 2)
    #print("Number of comparisons to make: " + str(nbOfComparisons))

    pairwiseComparison=generateAllComparisons(Sels, S, nbOfComparisons)

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

    return hypothesis

#This function estimates the number of violations in all the cube of R
#by randomly drawing tuples from the materialized cuboids (R included)
#It uses Hoeffding concentration inequality for bounding the number of draws according to a confidence interval
def estimateViolations(conn, meas, measBase, table, sel, cuboids, ranking, epsilon = 0.1, alpha = 0.1):
    #n is number of draws
    n = math.log(2 / alpha, 10) / (2 * epsilon * epsilon)
    n = math.ceil(n)

    estimates=0

    for i in range(n):
        nCuboid = random.randint(1, len(cuboids)) #+1 is for R itself
        print("nCuboid: ",nCuboid)
        #if nCuboid == len(cuboids):
        #    #draw in R
        #    #TODO draw tuples where only diff is on sel attribute!
        #    tuples=getSample(conn, measBase, table, sel, 5)
        #else:
        #draw in cuboid nCuboid
        view=cuboids[nCuboid][0]
        tuples=getSample(conn, "avg", view, sel, 5)
        if checkViolation(tuples, ranking) == True:
            estimates=estimates+1

    return n, estimates

#checks if tuples violate ranking, return True if this is the case
def checkViolation(tuples, ranking):
    print("tuples: ", tuples)
    print("ranking: ", ranking)

    meast1 = tuples[0][1]
    meast2 = tuples[1][1]

    valuet1 = tuples[0][0]
    valuet2 = tuples[1][0]

    rankt1 = getRank(valuet1,ranking)
    rankt2 = getRank(valuet2, ranking)

    if meast1<meast2 and rankt1<rankt2:
        return True
    else:
        return False

def hoeffdingForRank(groupbyAtt, n, valsToSelect, limitedHyp):
    print("Size of confidence interval around p: " + str(epsilon))
    print("Probability is of making a mistake: " + str(alpha))

    # n queries enough according to Hoeffding
    print("n: " + str(n))

    # compute powerset of categorical attributes
    pwset = powerset(groupbyAtt)

    # empty group by set removed from powerset
    # since it WAS used to generate the hypothesis

    pwset.remove(())

    #print("Hypothesis is:" + str(hypothesis))

    nbTests = 0
    for i in range(n):

        # generate the random query
        # query, queryHyp, queryValues, queryExcept, strgb = generateRandomQuery(pwsert,hypothesis)
        queryValues, queryCountGb, queryCountExcept = generateRandomQuery(pwset, valsToSelect, limitedHyp)

        # strategy: use the db engine to check whether the hypothesis holds
        # cross join hypothesis to chosen group by set
        # compute the actual ranks of vals in the hypothesis for each group
        # except all the actual ranks with the hypothesis
        # remaining only the groups where hypothesis does not hold

        resultCountGb = execute_query(conn, queryCountGb)
        resultCountExcept = execute_query(conn, queryCountExcept)

        print("number of tuples checked: " + str(resultCountGb[0][0]))
        print("number of exceptions: " + str(resultCountExcept[0][0]))
        print("ratio is: " + str(resultCountExcept[0][0] / resultCountGb[0][0]))

        ratio = resultCountExcept[0][0] / resultCountGb[0][0]

        # keep actual ratio
        # could use Bernouilli random variables instead: 1 if less than 10% errors
        if ratio > threshold:
            H.append(ratio)
        #    H.append(0)
        else:
            H.append(ratio)
        #    H.append(1)

        #nbTests = nbTests + 1

    expectedValue = sum(H) / len(H)
    print("Expected value is: " + str(sum(H) / len(H)))
    return expectedValue

# returns true if tuple violates ranking
# TODO check!
def countViolations(conn, viewName, ranking):
    viewDef=dbStuff.getDefOfMV(conn, viewName)
    strgb=viewDef.split("GROUP BY ")[1].split(";")[0]
    queryHyp = (
            "select * from (select  * from " + viewName +  ") t1  cross join (values " + ranking + ") as t2 ")

    query = ("SELECT " + strgb + "," + sel + "," + meas + ", "
             + " rank () over ( partition by " + strgb + " order by " + meas + " desc ) as rank" +
             " FROM " + table +  " group by " + strgb + "," + sel + " ")

    #queryValues = ("SELECT measure FROM (SELECT " + strgb + "," + sel + "," + meas + " as measure FROM "+ table + " WHERE " + sel + " in " + str(valsToSelect) + " group by " + strgb + "," + sel + " ) x;")

    queryExcept = ("select " + strgb + "," + sel + ", rank from  (" + query + " ) t3 except all " + queryHyp + " ")

    #queryCountGb = ("select count(*) from (" + queryHyp + ") t4;")
    queryCountExcept = ("select count(*) from (" + queryExcept + ") t5;")
    return dbStuff.execute_query(conn, queryCountExcept)




#draws n views, return those having less than threshold violations
# TODO check!
def azuma(conn, n, threshold, ranking):
    print(n + " draws, you have "+ (100-n) +"% of chances to get " + math.sqrt(2*n*math.log(n))+ " cuboids with acceptable violations")
    res = dbStuff.getMVnames(conn)
    #print(res)
    for i in range(n):
        nb = random.randint(0, len(res) - 1)
        viewName = res[nb][0]
        print(viewName)
        res.remove(nb)
        view=dbStuff.execute_query(conn, "select * from "+viewName + ";")
        print(view)
        nbViolations=countViolations(view, ranking)
        tabView=[]
        if nbViolations < threshold:
            tabView.append(view)
        return tabView




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

    print("nb of True in rejected: ", utilities.nbTrueInList(rejected))

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

    print("Number of BH corrections: ", nbChanges)

    print("nb non zeros after corrections: ", utilities.countNonZeros(pairwiseComp2))

    return pairwiseComp2


def generateComparisonsWithMergeSort(Sels, S):
    # compute Claire statistics for all pairs
    claireTab = []
    for i in range(1, len(S)):
        for j in range(i, len(S)):
            b = claireStat(S[i - 1][2], S[j][2], S[i - 1][1], S[j][1])
            claireTab.append((S[i - 1][0], S[j][0], b, S[i - 1][3], S[j][3]))

    #claireComp = [(x[0],x[1],x[2]) for x in claireTab]
    #print("claireComp: ", claireTab)
    print("merge: ", merge_sort(Sels, claireTab))
    print("pairwise: ", rankingFromPairwise.pairwiseComparison)

    return rankingFromPairwise.pairwiseComparison


def generateAllComparisons(Sels, S, nbOfComparisons):
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
    return pairwiseComparison


# returns the rank of value in ranking
# returns 0 if value not found
def getRank(value, ranking):
    for r in ranking:
        if r[0] == value:
            rank=r[1]
    return rank

