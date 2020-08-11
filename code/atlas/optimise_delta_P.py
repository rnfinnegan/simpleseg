import SimpleITK as sitk
import numpy as np
import os
import gc
import pandas as pd

import matplotlib.pyplot as plt

from robstools import LabelComparison, LabelFusion

def quickOptimiseCost(metricFunction, referenceMetric, probabilityImage, pOptimal=0.5, tolerance = 0.01, metricType='min', show=False, doseMap=False):

    n=0

    xList  = list(np.linspace(0, 0.3, 5))
    delta = xList[1]-xList[0]
    if doseMap:
        fList = [(referenceMetric - metricFunction(probabilityImage, pOptimal, p, doseMap))**2 for p in xList]
    else:
        fList = [(referenceMetric - metricFunction(probabilityImage, pOptimal, p))**2 for p in xList]

    if metricType=='min':
        xn = xList[np.argmin(fList)]
    elif metricType=='max':
        xn = xList[np.argmax(fList)]

    fn = np.min(fList)
    f0 = fn

    print("  n = {0} | DP = {1:.3f}, Metric = {2:.3f}".format(n, xn, fn))

    while np.abs( (fn-f0) ) > tolerance or n==0:

        n += 1
        f0 = fn

        xNew = [xn-3*delta/4, xn-delta/2, xn-delta/4, xn+delta/4, xn+delta/2, xn+3*delta/4]

        if doseMap:
            fNew = [(referenceMetric - metricFunction(probabilityImage, pOptimal, p, doseMap))**2 for p in xNew]
        else:
            fNew = [(referenceMetric - metricFunction(probabilityImage, pOptimal, p))**2 for p in xNew]

        xList = xList + xNew
        fList = fList + fNew

        if metricType=='min':
            xn = xList[np.argmin(fList)]
            fn = np.min(fList)
        elif metricType=='max':
            xn = xList[np.argmax(fList)]
            fn = np.max(fList)

        delta /= 4

        print("  n = {0} | DP = {1:.3f}, Metric = {2:.3f}".format(n, xn, fn))

    if show:
        fig, ax = plt.subplots(1,1)
        ax.scatter(xList,fList, c='k')
        ax.set_xlim(0,1)
        ax.set_xlabel("Probability Difference (from Optimal)")
        ax.set_ylabel("Metric Value")
        ax.grid()
        plt.show()

    return xn, fn

def metricFunction(probabilityImage, pThresh, pRange):
    lower = LabelFusion.processProbabilityImage(probabilityImage, max([0.05,pThresh-pRange]))
    upper = LabelFusion.processProbabilityImage(probabilityImage, min([0.95,pThresh+pRange]))

    MASD = LabelComparison.metric_MASD(lower, upper)

    return MASD

def doseVariation(probabilityImage, pThresh, pRange, doseMap):
    lower = LabelFusion.processProbabilityImage(probabilityImage, max([0.05,pThresh-pRange]))
    upper = LabelFusion.processProbabilityImage(probabilityImage, min([0.95,pThresh+pRange]))

    dArr = sitk.GetArrayFromImage(doseMap)
    lArr = sitk.GetArrayFromImage(lower)
    uArr = sitk.GetArrayFromImage(upper)

    dLower = np.mean( dArr[np.where(lArr)] )
    dUpper = np.mean( dArr[np.where(uArr)] )

    return dLower - dUpper

idList = ['A_2CF7FA45B81784','A_8E9603E0367472','A_3790D87A0FA1CD',
          'A_9611A8A3B59DA0','A_395474114B4C70','A_CBD6D1241EDB54',
          'A_5C7C972446070F','A_D9E41C39A4A9A2','A_77ADFE6E74BAB8',
          'A_E18A73126AE10A','A_F97ADCA9BE58CC','A_85E34D88B2C2D2',
          'A_8020C64E85ABF6','A_EA6432A5A43375','A_84DEDFEC1CA427']

f = open('ProbabilityThresholdData-DOSE-VARIATION.csv','w')
f.write('TARGET_ID POPT PMIN PMAX METRIC_MANUAL METRIC_AUTO\n')

df = pd.read_csv('../../../Processing/AustralianAtlas/MeanOptimalThresholds.csv', delimiter=',')

for targetId in idList:
    print(targetId)

    """
    MASD
    """
    # manualProbImage = sitk.ReadImage("../../../Data/Case_{0}/Structures/Case_{0}_COR_MANUAL_PROBABILITY_CROP.nii.gz".format(targetId))
    # manualInner = LabelFusion.processProbabilityImage(manualProbImage, 1)
    # manualOuter = LabelFusion.processProbabilityImage(manualProbImage, 0.01)
    #
    # referenceMetric = LabelComparison.metric_MASD(manualInner, manualOuter)

    """
    Dose variation
    """
    doseMap = sitk.ReadImage(f"../../../Data/Case_{targetId}/Doses/Case_{targetId}_dose_FRACTION_RESAMPLED_CROP.nii.gz")
    doseMapArr = sitk.GetArrayFromImage(doseMap)

    doseList = []

    for obs in range(9):
        s = sitk.ReadImage(f"../../../Data/Case_{targetId}/Structures/Case_{targetId}_COR_{obs}_CROP.nii.gz")
        sArr = sitk.GetArrayFromImage(s)
        d = np.mean(doseMapArr[np.where(sArr)])

        print(f"  {obs}: {d:.2f}")

        doseList.append(d)

    referenceMetric = np.subtract(*np.percentile(doseList,[97.5,2.5]))
    referenceMetric = 1.96*np.std(doseList)

    print('  METRIC = {0:.3f}'.format(referenceMetric))

    autoProbImage   = sitk.ReadImage("../../../Processing/AustralianAtlas/LabelFusion/LeaveOut{0}/Case_{0}_WHOLEHEART_PROBABILITY_CORR.nii.gz".format(targetId))
    pOptimal = df.loc[df['STRUCTURE']=='WHOLEHEART']['P_OPTIMAL'].values[0]

    DP, metricOpt = quickOptimiseCost(doseVariation, referenceMetric, autoProbImage, pOptimal=pOptimal, metricType='min', show=False, doseMap=doseMap)
    comparisonMetric = doseVariation(autoProbImage, pOptimal, DP, doseMap)
    # metricFunction(autoProbImage, pOptimal, DP)

    print('    P_OPT = {0:.3f} [+/- {1:.3f}]\n    METRIC (M/A) = {2:.3f}/{3:.3f}'.format(pOptimal, DP, referenceMetric, comparisonMetric))

    f.write('{0} {1} {2} {3} {4} {5}\n'.format(targetId, pOptimal, pOptimal-DP, pOptimal+DP, referenceMetric, comparisonMetric))
    f.flush()

    gc.collect()

f.close()
