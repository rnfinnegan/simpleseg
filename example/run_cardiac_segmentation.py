"""
Service to run multi-atlas based thorax CT segmentation.

Rob Finnegan
July, 2020
"""

# Read modules
import SimpleITK as sitk
from loguru import logger

from simpleseg.code.cardiac.run import CARDIAC_SETTINGS_DEFAULTS, run_cardiac_segmentation

# Set the settings using the default
settings = CARDIAC_SETTINGS_DEFAULTS

# Structures we might like
structure_list = ["HEART","LUNG_L","LUNG_R","ESOPHAGUS","SPINALCORD"]

settings['atlasSettings']['atlasPath'] = '../data'
settings['atlasSettings']['atlasStructures'] = structure_list
settings['labelFusionSettings']['optimalThreshold'] = {s:0.5 for s in structure_list}
settings['IARSettings']['referenceStructure'] = 'HEART'

settings["deformableSettings"]["isotropicResample"] =  True
settings["deformableSettings"]["resolutionStaging"] =  [32,8,4]  # specify voxel size (mm) since isotropic_resample is set
settings["deformableSettings"]["iterationStaging"] =  [50,25,25]

# Run the segmentation pipeline
if __name__ == "__main__":
    test_id = 'Train-S1-002'
    test_image = sitk.ReadImage(f'../data/{test_id}/Images/{test_id}.nii.gz')

    output_name = './output_{0}.nii.gz'

    auto_structs = run_cardiac_segmentation(image=test_image, settings=settings)

    logger.info('Writing outputs')
    for struct_name in list(auto_structs.keys()):
        logger.info(f'  > {struct_name}')
        sitk.WriteImage(auto_structs[struct_name], output_name.format(struct_name))