"""
Service to run multi-atlas based thorax CT segmentation.

Rob Finnegan
July, 2020
"""

# Read modules
import os
import SimpleITK as sitk
from loguru import logger

from simpleseg.code.cardiac.run import CARDIAC_SETTINGS_DEFAULTS, run_cardiac_segmentation

# Set the settings
settings = {
    "atlasSettings": {
        "atlasIdList": ["002","003","005","006","007","008"],
        "atlasStructures": ["Heart","Lung-Left","Lung-Right","Esophagus","Spinal-Cord"],
        "atlasPath": '../data/NIFTI_CONVERTED',
        "atlasImageFormat": "Study_{0}/Study_{0}.nii.gz",
        "atlasLabelFormat": "Study_{0}/Study_{0}_{1}.nii.gz",
    },
    "autoCropSettings": {"expansion": [0,0,0],},
    "intialRegSettings": {
        "initialReg": "Similarity",
        "options": {
            "shrinkFactors": [16, 8, 4],
            "smoothSigmas": [0, 0, 0],
            "samplingRate": 0.75,
            "defaultValue": -1024,
            "optimiser": "gradient_descent",
            "numberOfIterations": 50,
            "finalInterp": sitk.sitkBSpline,
            "metric": "mean_squares",
        },
        "trace": False,
        "guideStructure": False,
    },
    "deformableSettings": {
        "isotropicResample": True,
        "resolutionStaging": [16, 8, 2],  # specify voxel size (mm) since isotropic_resample is set
        "iterationStaging": [5, 5, 5],
        "smoothingSigmas": [0, 0, 0],
        "ncores": 8,
        "trace": False,
    },
    "IARSettings": {
        "referenceStructure": "Heart",
        "smoothDistanceMaps": True,
        "smoothSigma": 1,
        "zScoreStatistic": "MAD",
        "outlierMethod": "IQR",
        "outlierFactor": 1.5,
        "minBestAtlases": 5,
        "project_on_sphere": False,
    },
    "labelFusionSettings": {
        "voteType": "unweighted",
        "voteParams": {},
        "optimalThreshold":    {"Heart":0.5,
                                "Lung-Left":0.5,
                                "Lung-Right":0.5,
                                "Esophagus":0.5,
                                }
    },
    "vesselSpliningSettings": {
        "vesselNameList": ['Spinal-Cord'],
        "vesselRadius_mm": {'Spinal-Cord': 6},
        "spliningDirection": {'Spinal-Cord': 'z'},
        "stopCondition": {'Spinal-Cord': 'count'},
        "stopConditionValue": {'Spinal-Cord': 2}
    },
}

# Run the segmentation pipeline
if __name__ == "__main__":
    test_id = '010'
    output_name = f'./output/Study_{test_id}/Study_{test_id}_{{0}}_AUTO.nii.gz'
    test_image = sitk.ReadImage(f'../data/NIFTI_CONVERTED/Study_{test_id}/Study_{test_id}.nii.gz')

    try:
        os.mkdir('./output')
    except:
        None

    try:
        os.mkdir(f'./output/Study_{test_id}/')
    except:
        None

    auto_structs = run_cardiac_segmentation(image=test_image, settings=settings)

    logger.info('Writing outputs')
    for struct_name in list(auto_structs.keys()):
        logger.info(f'  > {struct_name}')
        sitk.WriteImage(auto_structs[struct_name], output_name.format(struct_name))