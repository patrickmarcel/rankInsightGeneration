import statistics
import plotStuff
from main import materializeViews,test



def testTimingsQuerySampleSize(nbruns,conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice, groupbyAtt, sel, measBase, meas, function,table, comparison, generateIndex,
                                                                           allComparisons, initsampleSize, sizeOfR, ratioOfQuerySample, cumulate):


    dictQuery={}
    dictSamp={}
    dictHyp={}
    dictVal={}

    paramTested = 'Query sample size'

    tabTest=(0.1, 0.25, 0.5)
    #tabTest = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)


    for j in range(nbruns):
        print("-----RUN: ", j)

        mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex)
        currentSample = {}


        for ratioOfQuerySample in tabTest:
        #for percentOfLattice in tabTest:
        # for initsampleSize in tabTest:
        # for nbAdomVals in range(2,10):

            print("--- TESTING VALUE:", ratioOfQuerySample)

            sampleSize = initsampleSize * sizeOfR

            #benTab = []
            #samplingTab = []
            #hypoTab = []
            #validTab = []

            #for i in range(nbOfRuns):
            #    print("-----RUN: ", i)

            queryTime, samplingTime, hypothesisTime, validationTime = test(conn, nbAdomVals, prefs, ratioViolations,
                                                                               proba, error,
                                                                               percentOfLattice, groupbyAtt, sel, measBase, meas,
                                                                               function,
                                                                               table, sampleSize, comparison, generateIndex,
                                                                               allComparisons, ratioOfQuerySample, mvnames,
                                                                               aggQueries,
                                                                               currentSample, cumulate=True)

            if str(ratioOfQuerySample) in dictQuery:
                dictQuery[str(ratioOfQuerySample)].extend([queryTime])
            else:
                dictQuery[str(ratioOfQuerySample)]=[queryTime]
            if str(ratioOfQuerySample) in dictSamp:
                dictSamp[str(ratioOfQuerySample)].extend([samplingTime])
            else:
                dictSamp[str(ratioOfQuerySample)]=[samplingTime]
            if str(ratioOfQuerySample) in dictHyp:
                dictHyp[str(ratioOfQuerySample)].extend([hypothesisTime])
            else:
                dictHyp[str(ratioOfQuerySample)]=[hypothesisTime]
            if str(ratioOfQuerySample) in dictVal:
                dictVal[str(ratioOfQuerySample)].extend([validationTime])
            else:
                dictVal[str(ratioOfQuerySample)]=[validationTime]
            #benTab.append(queryTime)
            #samplingTab.append(samplingTime)
            #hypoTab.append(hypothesisTime)
            #validTab.append(validationTime)

        #tabres[j]=dictTest


    meanQ=[]
    stdevQ=[]
    meanSamp=[]
    stdevSamp=[]
    meanHypo=[]
    stdevHypo=[]
    meanValid=[]
    stdevValid=[]

    for i in tabTest:
        # query time
        meanQ.append(statistics.mean(dictQuery[str(i)]))
        stdevQ.append(statistics.stdev(dictQuery[str(i)]))
        # sampling time
        meanSamp.append(statistics.mean(dictSamp[str(i)]))
        stdevSamp.append(statistics.stdev(dictSamp[str(i)]))
        # hypothesis time
        meanHypo.append(statistics.mean(dictHyp[str(i)]))
        stdevHypo.append(statistics.stdev(dictHyp[str(i)]))
        # validation time
        meanValid.append(statistics.mean(dictVal[str(i)]))
        stdevValid.append(statistics.stdev(dictVal[str(i)]))

        #meanQ = statistics.mean(qtab)
        #meanSamp = statistics.mean(samplingTab)
        #meanHypo = statistics.mean(hypoTab)
        #meanValid = statistics.mean(validTab)


    data = [
        {'x': tabTest, 'y': meanQ, 'yerr': stdevQ, 'label': 'Aggregate queries time'},
        {'x': tabTest, 'y': meanSamp, 'yerr': stdevSamp, 'label': 'Sampling time'},
        {'x': tabTest, 'y': meanHypo, 'yerr': stdevHypo, 'label': 'Hypothesis time'},
        {'x': tabTest, 'y': meanValid, 'yerr': stdevValid, 'label': 'Validation time'}
    ]

    plotStuff.plot_curves_with_error_bars(data, x_label=paramTested, y_label='Time (s)',title='Times',scale='log')


def testAccuracyQuerySampleSize(nbruns,conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice, groupbyAtt, sel, measBase, meas, function,table, comparison, generateIndex,
                                                                           allComparisons, initsampleSize, sizeOfR, ratioCuboidOK, ratioOfQuerySample, cumulate):
    dictPred = {}
    dictBennet = {}
    dictError = {}
    dictWR = {}


    paramTested = 'Query sample size'
    # paramTested = 'Sample size'
    # paramTested = 'Percent of lattice'
    # tabTest=(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9, 1)
    tabTest = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
    # tabTest=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
    # tabTest=(5,10,20,50,75,100

    for j in range(nbruns):
        print("-----RUN: ", j)

        mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,
                                               generateIndex)
        currentSample = {}

        # for percentOfLattice in tabTest:
        # for initsampleSize in tabTest:
        for ratioOfQuerySample in tabTest:
            # for nbAdomVals in range(2,10):

            print("--- TESTING VALUE:", ratioOfQuerySample)

            sampleSize = initsampleSize * sizeOfR

            prediction, bennetError, realError, gtratio = test(conn, nbAdomVals, prefs, ratioViolations, proba, error,
                                                           percentOfLattice, groupbyAtt,
                                                           sel, measBase, meas, function, table, sampleSize, comparison,
                                                           generateIndex, allComparisons, ratioOfQuerySample, mvnames,
                                                           aggQueries, currentSample, cumulate=True)

            if str(ratioOfQuerySample) in dictPred:
                dictPred[str(ratioOfQuerySample)].extend([prediction])
            else:
                dictPred[str(ratioOfQuerySample)] = [prediction]
            if str(ratioOfQuerySample) in dictBennet:
                dictBennet[str(ratioOfQuerySample)].extend([bennetError])
            else:
                dictBennet[str(ratioOfQuerySample)] = [bennetError]
            if str(ratioOfQuerySample) in dictError:
                dictError[str(ratioOfQuerySample)].extend([realError])
            else:
                dictError[str(ratioOfQuerySample)] = [realError]

            print("Desired cuboid ratio is:", ratioCuboidOK, ". We predicted ratio of: ", prediction, ". Real ratio is: ",
                  gtratio)
            if (gtratio < ratioCuboidOK and prediction > ratioCuboidOK) or (
                    gtratio > ratioCuboidOK and prediction < ratioCuboidOK):
                nbWrongRanking = 1
            else:
                nbWrongRanking = 0
            #nbWrongRankingTab.append(nbWrongRanking)

            print("interval: [", prediction - bennetError, ",", prediction + bennetError, "]")
            print("user threshold:", ratioCuboidOK)
            if ratioCuboidOK >= prediction - bennetError and ratioCuboidOK <= prediction + bennetError:
                print("continue")
            else:
                print("WE CAN STOP")

            if str(ratioOfQuerySample) in dictWR:
                dictWR[str(ratioOfQuerySample)].extend([nbWrongRanking])
            else:
                dictWR[str(ratioOfQuerySample)] = [nbWrongRanking]

    meanPred = []
    stdevPred = []
    meanBennet = []
    stdevBennet = []
    meanError = []
    stdevError = []
    meanWRW = []
    stdevWR = []

    for i in tabTest:
        # prediction
        meanPred.append(statistics.mean(dictPred[str(i)]))
        stdevPred.append(statistics.stdev(dictPred[str(i)]))
        # Bennet
        meanBennet.append(statistics.mean(dictBennet[str(i)]))
        stdevBennet.append(statistics.stdev(dictBennet[str(i)]))
        # error
        meanError.append(statistics.mean(dictError[str(i)]))
        stdevError.append(statistics.stdev(dictError[str(i)]))
        # WR
        meanWRW.append(statistics.mean(dictWR[str(i)]))
        stdevWR.append(statistics.stdev(dictWR[str(i)]))

    data = [
        {'x': tabTest, 'y': meanBennet, 'yerr': stdevBennet, 'label': 'Bennet theoretical error'},
        {'x': tabTest, 'y': meanError, 'yerr': stdevError, 'label': 'real error'},
        {'x': tabTest, 'y': meanPred, 'yerr': stdevPred, 'label': 'prediction'},
        {'x': tabTest, 'y': meanWRW, 'yerr': stdevWR, 'label': 'unvalidated prediction'}
    ]

    plotStuff.plot_curves_with_error_bars(data, x_label=paramTested, y_label='Error',
                                title='prediction and errors')


def testTimingsLattice(conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice,
                groupbyAtt, sel, measBase, meas, function,table, comparison, generateIndex,
                allComparisons,
                initsampleSize,sizeOfR,nbOfRuns,
                ratioOfQuerySample, cumulate):
    listBennet = []
    devBennet = []
    listSampling = []
    devSampling = []
    listHypo = []
    devHypo = []
    listValid = []
    devValid = []

    paramTested='Percent of Lattice'
    #paramTested = 'Query sample size'

    tabTest=(0.1, 0.25, 0.5, 0.75)
    #tabTest = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)



    #for ratioOfQuerySample in tabTest:
    for percentOfLattice in tabTest:
        # for initsampleSize in tabTest:
        # for nbAdomVals in range(2,10):

        print("--- TESTING VALUE:", percentOfLattice)



        sampleSize = initsampleSize * sizeOfR

        benTab = []
        samplingTab = []
        hypoTab = []
        validTab = []

        for i in range(nbOfRuns):
            print("-----RUN: ", i)

            mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice,generateIndex)
            currentSample = {}

            queryTime, samplingTime, hypothesisTime, validationTime = test(conn, nbAdomVals, prefs, ratioViolations,
                                                                           proba, error,
                                                                           percentOfLattice, groupbyAtt, sel, measBase, meas,
                                                                           function,
                                                                           table, sampleSize, comparison, generateIndex,
                                                                           allComparisons, ratioOfQuerySample, mvnames,
                                                                           aggQueries,
                                                                           currentSample, cumulate=False)

            benTab.append(queryTime)
            samplingTab.append(samplingTime)
            hypoTab.append(hypothesisTime)
            validTab.append(validationTime)

        # resultRuns.append((percentOfLattice, bennetError, hypothesisTime, validationTime))
        meanBen = statistics.mean(benTab)
        meanSamp = statistics.mean(samplingTab)
        meanHypo = statistics.mean(hypoTab)
        meanValid = statistics.mean(validTab)

        if nbOfRuns == 1:
            stdevBen = 0
            stdevSamp = 0
            stdevHypo = 0
            stdevValid = 0
        else:
            stdevBen = statistics.stdev(benTab)
            stdevSamp = statistics.stdev(samplingTab)
            stdevHypo = statistics.stdev(hypoTab)
            stdevValid = statistics.stdev(validTab)

        listBennet.append(meanBen)
        devBennet.append(stdevBen)
        listSampling.append(meanSamp)
        devSampling.append(stdevSamp)
        listHypo.append(meanHypo)
        devHypo.append(stdevHypo)
        listValid.append(meanValid)
        devValid.append(stdevValid)

    data = [
        {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Aggregate queries time'},
        {'x': tabTest, 'y': listSampling, 'yerr': devSampling, 'label': 'Sampling time'},
        {'x': tabTest, 'y': listHypo, 'yerr': devHypo, 'label': 'Hypothesis time'},
        {'x': tabTest, 'y': listValid, 'yerr': devValid, 'label': 'Validation time'}
    ]

    plotStuff.plot_curves_with_error_bars(data, x_label=paramTested, y_label='Time (s)',title='Times',scale='log')





def testAccuracyInitSampleSize(conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice,
                groupbyAtt, sel, measBase, meas, function,table, comparison, generateIndex,
                allComparisons,
                initsampleSize,sizeOfR,nbOfRuns,ratioCuboidOK,
                ratioOfQuerySample, cumulate):
        listPred=[]
        devPred=[]
        listError=[]
        devError=[]
        listWR=[]
        devWR=[]
        listBennet=[]
        devBennet=[]

        #paramTested = 'Query sample size'
        paramTested = 'Sample size'
        #paramTested = 'Percent of lattice'
        #tabTest=(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9, 1)
        tabTest=(0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1)
        #tabTest=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
        #tabTest=(5,10,20,50,75,100)

        mvnames,aggQueries=materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex)
        currentSample={}

        #for percentOfLattice in tabTest:
        for initsampleSize in tabTest:
        #for ratioOfQuerySample in tabTest:
        #for nbAdomVals in range(2,10):

            print("--- TESTING VALUE:",ratioOfQuerySample)

            sampleSize = initsampleSize * sizeOfR

            predictionTab=[]
            realErrorTab=[]
            nbWrongRankingTab=[]
            bennetTab = []

            for i in range(nbOfRuns):

                print("-----RUN: ",i)
                prediction,bennetError,realError,gtratio=test(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice, groupbyAtt,
                                                              sel, measBase, meas, function,table, sampleSize, comparison,generateIndex,allComparisons,ratioOfQuerySample,mvnames,aggQueries,currentSample,cumulate=False)
                #resultRuns.append((percentOfLattice,prediction,bennetError,realError))

                predictionTab.append(prediction)
                bennetTab.append(bennetError)
                realErrorTab.append(realError)
                print("Desired cuboid ratio is:",ratioCuboidOK,". We predicted ratio of: ",prediction,". Real ratio is: ",gtratio)
                #if gtratio <ratioCuboidOK:
                if (gtratio < ratioCuboidOK and prediction > ratioCuboidOK) or (gtratio > ratioCuboidOK and prediction < ratioCuboidOK):
                    nbWrongRanking=1
                else:
                    nbWrongRanking = 0
                nbWrongRankingTab.append(nbWrongRanking)

                print("interval: [",prediction-bennetError,",",prediction+bennetError,"]")
                print("user threshold:",ratioCuboidOK)
                if ratioCuboidOK >= prediction-bennetError and ratioCuboidOK <= prediction+bennetError:
                    print("continue")
                else:
                    print("WE CAN STOP")

            meanPred=statistics.mean(predictionTab)
            meanBen= statistics.mean(bennetTab)
            meanError=statistics.mean(realErrorTab)
            meanWRTab = statistics.mean(nbWrongRankingTab)

            if nbOfRuns==1:
                stdevPred = 0
                stdevBen = 0
                stdevError = 0
                stdevWRTab = 0
            else:
                stdevPred = statistics.stdev(predictionTab)
                stdevBen = statistics.stdev(bennetTab)
                stdevError = statistics.stdev(realErrorTab)
                stdevWRTab = statistics.stdev(nbWrongRankingTab)


            listPred.append(meanPred)
            devPred.append(stdevPred)
            listBennet.append(meanBen)
            devBennet.append(stdevBen)
            listError.append(meanError)
            devError.append(stdevError)
            listWR.append(meanWRTab)
            devWR.append(stdevWRTab)

        # Example usage:
        data = [
            {'x': tabTest, 'y':listPred,  'yerr': devPred, 'label': 'prediction'},
            {'x': tabTest, 'y': listError, 'yerr': devError, 'label': 'real error'},
            {'x': tabTest, 'y': listWR, 'yerr': devWR, 'label': 'unvalidated prediction'},
            {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Bennet theoretical error'}
        ]

        plotStuff.plot_curves_with_error_bars(data, x_label=paramTested, y_label='Error',
                                    title='prediction and errors')



def testTimingsOLD(conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice,
                groupbyAtt, sel, measBase, meas, function,table, comparison, generateIndex,
                allComparisons,
                initsampleSize,sizeOfR,nbOfRuns,
                ratioOfQuerySample, cumulate):
    listBennet = []
    devBennet = []
    listSampling = []
    devSampling = []
    listHypo = []
    devHypo = []
    listValid = []
    devValid = []

    paramTested='Percent of Lattice'
    #paramTested = 'Query sample size'

    tabTest=(0.1, 0.25, 0.5, 0.75)
    #tabTest = (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

    mvnames, aggQueries = materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex)
    currentSample = {}

    #for ratioOfQuerySample in tabTest:
    for percentOfLattice in tabTest:
        # for initsampleSize in tabTest:
        # for nbAdomVals in range(2,10):

        print("--- TESTING VALUE:", ratioOfQuerySample)

        sampleSize = initsampleSize * sizeOfR

        benTab = []
        samplingTab = []
        hypoTab = []
        validTab = []

        for i in range(nbOfRuns):
            print("-----RUN: ", i)

            queryTime, samplingTime, hypothesisTime, validationTime = test(conn, nbAdomVals, prefs, ratioViolations,
                                                                           proba, error,
                                                                           percentOfLattice, groupbyAtt, sel, measBase, meas,
                                                                           function,
                                                                           table, sampleSize, comparison, generateIndex,
                                                                           allComparisons, ratioOfQuerySample, mvnames,
                                                                           aggQueries,
                                                                           currentSample, cumulate=True)

            benTab.append(queryTime)
            samplingTab.append(samplingTime)
            hypoTab.append(hypothesisTime)
            validTab.append(validationTime)

        # resultRuns.append((percentOfLattice, bennetError, hypothesisTime, validationTime))
        meanBen = statistics.mean(benTab)
        meanSamp = statistics.mean(samplingTab)
        meanHypo = statistics.mean(hypoTab)
        meanValid = statistics.mean(validTab)

        if nbOfRuns == 1:
            stdevBen = 0
            stdevSamp = 0
            stdevHypo = 0
            stdevValid = 0
        else:
            stdevBen = statistics.stdev(benTab)
            stdevSamp = statistics.stdev(samplingTab)
            stdevHypo = statistics.stdev(hypoTab)
            stdevValid = statistics.stdev(validTab)

        listBennet.append(meanBen)
        devBennet.append(stdevBen)
        listSampling.append(meanSamp)
        devSampling.append(stdevSamp)
        listHypo.append(meanHypo)
        devHypo.append(stdevHypo)
        listValid.append(meanValid)
        devValid.append(stdevValid)

    data = [
        {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Aggregate queries time'},
        {'x': tabTest, 'y': listSampling, 'yerr': devSampling, 'label': 'Sampling time'},
        {'x': tabTest, 'y': listHypo, 'yerr': devHypo, 'label': 'Hypothesis time'},
        {'x': tabTest, 'y': listValid, 'yerr': devValid, 'label': 'Validation time'}
    ]

    plotStuff.plot_curves_with_error_bars(data, x_label=paramTested, y_label='Time (s)',title='Times',scale='log')




def testAccuracyOLD(conn, nbAdomVals, prefs, ratioViolations,proba, error, percentOfLattice,
                groupbyAtt, sel, measBase, meas, function,table, comparison, generateIndex,
                allComparisons,
                initsampleSize,sizeOfR,nbOfRuns,ratioCuboidOK,
                ratioOfQuerySample, cumulate):
        listPred=[]
        devPred=[]
        listError=[]
        devError=[]
        listWR=[]
        devWR=[]
        listBennet=[]
        devBennet=[]

        paramTested = 'Query sample size'
        #paramTested = 'Sample size'
        #paramTested = 'Percent of lattice'
        #tabTest=(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9, 1)
        tabTest=(0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1)
        #tabTest=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
        #tabTest=(5,10,20,50,75,100)

        mvnames,aggQueries=materializeViews(conn, groupbyAtt, sel, measBase, function, table, percentOfLattice, generateIndex,nbOfRuns)
        currentSample={}

        #for percentOfLattice in tabTest:
        #for initsampleSize in tabTest:
        for ratioOfQuerySample in tabTest:
        #for nbAdomVals in range(2,10):

            print("--- TESTING VALUE:",ratioOfQuerySample)

            sampleSize = initsampleSize * sizeOfR

            predictionTab=[]
            realErrorTab=[]
            nbWrongRankingTab=[]
            bennetTab = []

            for i in range(nbOfRuns):

                print("-----RUN: ",i)
                prediction,bennetError,realError,gtratio=test(conn, nbAdomVals, prefs, ratioViolations, proba, error, percentOfLattice, groupbyAtt,
                                                              sel, measBase, meas, function,table, sampleSize, comparison,generateIndex,allComparisons,ratioOfQuerySample,mvnames,aggQueries,currentSample,cumulate=True)
                #resultRuns.append((percentOfLattice,prediction,bennetError,realError))

                predictionTab.append(prediction)
                bennetTab.append(bennetError)
                realErrorTab.append(realError)
                print("Desired cuboid ratio is:",ratioCuboidOK,". We predicted ratio of: ",prediction,". Real ratio is: ",gtratio)
                #if gtratio <ratioCuboidOK:
                if (gtratio < ratioCuboidOK and prediction > ratioCuboidOK) or (gtratio > ratioCuboidOK and prediction < ratioCuboidOK):
                    nbWrongRanking=1
                else:
                    nbWrongRanking = 0
                nbWrongRankingTab.append(nbWrongRanking)

                print("interval: [",prediction-bennetError,",",prediction+bennetError,"]")
                print("user threshold:",ratioCuboidOK)
                if ratioCuboidOK >= prediction-bennetError and ratioCuboidOK <= prediction+bennetError:
                    print("continue")
                else:
                    print("WE CAN STOP")

            meanPred=statistics.mean(predictionTab)
            meanBen= statistics.mean(bennetTab)
            meanError=statistics.mean(realErrorTab)
            meanWRTab = statistics.mean(nbWrongRankingTab)

            if nbOfRuns==1:
                stdevPred = 0
                stdevBen = 0
                stdevError = 0
                stdevWRTab = 0
            else:
                stdevPred = statistics.stdev(predictionTab)
                stdevBen = statistics.stdev(bennetTab)
                stdevError = statistics.stdev(realErrorTab)
                stdevWRTab = statistics.stdev(nbWrongRankingTab)


            listPred.append(meanPred)
            devPred.append(stdevPred)
            listBennet.append(meanBen)
            devBennet.append(stdevBen)
            listError.append(meanError)
            devError.append(stdevError)
            listWR.append(meanWRTab)
            devWR.append(stdevWRTab)

        # Example usage:
        data = [
            {'x': tabTest, 'y':listPred,  'yerr': devPred, 'label': 'prediction'},
            {'x': tabTest, 'y': listError, 'yerr': devError, 'label': 'real error'},
            {'x': tabTest, 'y': listWR, 'yerr': devWR, 'label': 'unvalidated prediction'},
            {'x': tabTest, 'y': listBennet, 'yerr': devBennet, 'label': 'Bennet theoretical error'}
        ]

        plotStuff.plot_curves_with_error_bars(data, x_label=paramTested, y_label='Error',
                                    title='prediction and errors')
