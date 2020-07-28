"""
Service to run calculate label similarity metrics

Rob Finnegan
December 9, 2019
"""

import SimpleITK as sitk
import numpy as np

settings = {
    "labels": ["HEART", "LUNG_L", "LUNG_R", "SPINALCORD", "ESOPHAGUS"],
    "manual_dir" : "./data",
    "auto_dir" : "./output"
}

def surfaceMetrics(imFixed, imMoving, verbose=False):
    """
    HD, meanSurfDist, medianSurfDist, maxSurfDist, stdSurfDist
    """
    hausdorffDistance = sitk.HausdorffDistanceImageFilter()
    hausdorffDistance.Execute(imFixed, imMoving)
    HD = hausdorffDistance.GetHausdorffDistance()

    meanSDList = []
    maxSDList = []
    stdSDList = []
    medianSDList = []
    numPoints = []
    for (imA, imB) in ((imFixed, imMoving), (imMoving, imFixed)):

        labelIntensityStat = sitk.LabelIntensityStatisticsImageFilter()
        referenceDistanceMap = sitk.Abs(sitk.SignedMaurerDistanceMap(imA, squaredDistance=False, useImageSpacing=True))
        movingLabelContour = sitk.LabelContour(imB)
        labelIntensityStat.Execute(movingLabelContour,referenceDistanceMap)

        meanSDList.append(labelIntensityStat.GetMean(1))
        maxSDList.append(labelIntensityStat.GetMaximum(1))
        stdSDList.append(labelIntensityStat.GetStandardDeviation(1))
        medianSDList.append(labelIntensityStat.GetMedian(1))

        numPoints.append(labelIntensityStat.GetNumberOfPixels(1))

    if verbose:
        print("        Boundary points:  {0}  {1}".format(numPoints[0], numPoints[1]))

    meanSurfDist = np.dot(meanSDList, numPoints)/np.sum(numPoints)
    maxSurfDist = np.max(maxSDList)
    stdSurfDist = np.sqrt( np.dot(numPoints,np.add(np.square(stdSDList), np.square(np.subtract(meanSDList,meanSurfDist)))) )
    medianSurfDist = np.mean(medianSDList)

    resultDict = {}
    resultDict['hausdorffDistance'] = HD
    resultDict['meanSurfaceDistance'] = meanSurfDist
    resultDict['medianSurfaceDistance'] = medianSurfDist
    resultDict['maximumSurfaceDistance'] = maxSurfDist
    resultDict['sigmaSurfaceDistance'] = stdSurfDist

    return resultDict

def volumeMetrics(imFixed, imMoving):
    """
    DSC, VolOverlap, FracOverlap, TruePosFrac, TrueNegFrac, FalsePosFrac, FalseNegFrac
    """
    arrFixed = sitk.GetArrayFromImage(imFixed).astype(bool)
    arrMoving = sitk.GetArrayFromImage(imMoving).astype(bool)

    arrInter = arrFixed & arrMoving
    arrUnion = arrFixed | arrMoving

    voxVol = np.product(imFixed.GetSpacing())/1000. # Conversion to cm^3

    # 2|A & B|/(|A|+|B|)
    DSC =  (2.0*arrInter.sum())/(arrFixed.sum()+arrMoving.sum())

    #  |A & B|/|A | B|
    FracOverlap = arrInter.sum()/arrUnion.sum().astype(float)
    VolOverlap = arrInter.sum() * voxVol

    TruePos = arrInter.sum()
    TrueNeg = (np.invert(arrFixed) & np.invert(arrMoving)).sum()
    FalsePos = arrMoving.sum()-TruePos
    FalseNeg = arrFixed.sum()-TruePos

    #
    TruePosFrac = (1.0*TruePos)/(TruePos+FalseNeg)
    TrueNegFrac = (1.0*TrueNeg)/(TrueNeg+FalsePos)
    FalsePosFrac = (1.0*FalsePos)/(TrueNeg+FalsePos)
    FalseNegFrac = (1.0*FalseNeg)/(TruePos+FalseNeg)

    resultDict = {}
    resultDict['DSC'] = DSC
    resultDict['volumeOverlap'] = VolOverlap
    resultDict['fractionOverlap'] = FracOverlap
    resultDict['truePositiveFraction'] = TruePosFrac
    resultDict['trueNegativeFraction'] = TrueNegFrac
    resultDict['falsePositiveFraction'] = FalsePosFrac
    resultDict['falseNegativeFraction'] = FalseNegFrac

    return resultDict

if __name__ == '__main__':

    labels = settings['labels']
    manual_dir = settings['manual_dir']
    auto_dir = settings['auto_dir']

    test_id = 'Train-S1-009'

    # Read manual & auto structures
    manual_struct = {i:sitk.ReadImage(f'{manual_dir}/{test_id}/Structures/{test_id}_{i}.nii.gz') for i in labels}
    auto_struct = {i:sitk.ReadImage(f'{auto_dir}/{test_id}/{test_id}_{i}_AUTO.nii.gz') for i in labels}

    data = []

    # Compare labels
    for i,struct in enumerate(labels):
        volume_metrics = volumeMetrics(manual_struct[struct], auto_struct[struct])
        surface_metrics = surfaceMetrics(manual_struct[struct], auto_struct[struct])

        data.append([struct, volume_metrics['DSC'], surface_metrics['meanSurfaceDistance'], surface_metrics['hausdorffDistance']])


    # Save to file
    np.savetxt(f'Output-{test_id}.csv', data, fmt='%s', delimiter=',', newline='\n', header='struct,dsc,masd,hd')
