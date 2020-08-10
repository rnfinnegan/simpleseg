"""
Service to display output of automatic thorax CT segmentation

Rob Finnegan
July, 2020
"""

import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt

### Read in the data

manual_directory = '../data/'



### Step 1: Some nice figures

    image = sitk.ReadImage(imageBase.format(caseId))
    ndaImage = sitk.GetArrayFromImage(image)

    manHeartImage = sitk.ReadImage(structBaseMan.format(caseId))
    manHeartArr = sitk.GetArrayFromImage(manHeartImage)

    z,y,x = np.array(measurements.center_of_mass(manHeartArr), dtype=np.int)

    (AxSize, CorSize, SagSize) = ndaImage.shape
    spPlane, _, spSlice = image.GetSpacing()
    asp = (spSlice/spPlane)

    fSize = (6, 6.0*(asp*AxSize+CorSize)/(1.0*SagSize+CorSize))

    fig, ((axAx, blank), (axCor, axSag)) = plt.subplots(2, 2, figsize=fSize,gridspec_kw = {'height_ratios':[CorSize,asp*AxSize],'width_ratios': [SagSize,CorSize]})
    blank.axis('off')

    imAx = axAx.imshow(ndaImage[z,:,:], cmap=plt.cm.gray)
    axAx.axis('off')

    imCor = axCor.imshow(ndaImage[:,y,:], origin='lower', aspect=asp, cmap=plt.cm.gray)
    axCor.axis('off')

    imSag = axSag.imshow(ndaImage[:,:,x], origin='lower', aspect=asp, cmap=plt.cm.gray)
    axSag.axis('off')

    vmin = -250; window=500
    imAx.set_clim(vmin,vmin+window)
    imCor.set_clim(vmin,vmin+window)
    imSag.set_clim(vmin,vmin+window)

    fig.subplots_adjust(left=0, right=1, wspace=0.01, hspace=0.01, top=1, bottom=0)
    fig.savefig(name)
    plt.close('all')



### Step 2: Some quantitative results