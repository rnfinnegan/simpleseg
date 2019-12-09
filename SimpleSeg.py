"""
Service to run multi-atlas based cardiac segmentation.

Rob Finnegan
December 9, 2019
"""
import os

import SimpleITK as sitk
import numpy as np
from loguru import logger

from simpleseg.code.atlas.registration import (
    initial_registration,
    transform_propagation,
    fast_symmetric_forces_demons_registration,
    apply_field,
)

from simpleseg.code.atlas.label import (
    compute_weight_map,
    combine_labels,
    process_probability_image,
)

from simpleseg.code.atlas.iterative_atlas_removal import run_iar

from simpleseg.code.thorax.thoracic_utils import (
    AutoLungSegment,
    CropImage,
    vesselSplineGeneration,
)

settings = {
    "outputSettings": {
        "outputDir": "./output"
    },
    "atlasSettings": {
        "atlasIdList": ["Train-S1-001", "Train-S1-002", "Train-S1-003", "Train-S1-004", "Train-S1-005", "Train-S1-006", "Train-S1-007", "Train-S1-008"],
        "atlasStructures": ["HEART", "LUNG_L", "LUNG_R", "SPINALCORD", "ESOPHAGUS"],
        "atlasPath": "./data",
    },
    "autoCropSettings": {
        "expansion": 25
    },
    "rigidSettings": {
        "initialReg": "Rigid",
        "options": {
            "shrinkFactors": [8, 2, 1],
            "smoothSigmas": [8, 2, 1],
            "samplingRate": 0.25,
            "finalInterp": sitk.sitkBSpline,
        },
        "trace": True,
        "guideStructure": False,
    },
    "deformableSettings": {
        "resolutionStaging": [8,4,2,1],
        "iterationStaging": [150,100,50,25],
        "ncores": 8,
        "trace": True,
    },
    "IARSettings": {
        "referenceStructure": "WHOLEHEART",
        "smoothDistanceMaps": True,
        "smoothSigma": 1,
        "zScoreStatistic": "MAD",
        "outlierMethod": "IQR",
        "outlierFactor": 1.5,
        "minBestAtlases": 3,
        "project_on_sphere":False,
    },
    "labelFusionSettings": {"voteType": "local",
                            "optimalThreshold":
                                {
                                "HEART": 0.5,
                                "LUNG_L": 0.5,
                                "LUNG_R": 0.5,
                                "ESOPHAGUS": 0.5,
                                "SPINALCORD": 0.5
                                }
                            },
    "vesselSpliningSettings": {
        "vesselNameList": [],
        "vesselRadius_mm": {},
        "spliningDirection": {},
        "stopCondition": {},
        "stopConditionValue": {},
    }
}

def SimpleSeg(image, output_name, settings=settings):
    """
    Implements the simpleseg framework to provide cardiac atlas based segmentation.

    """

    logger.info("Running Cardiac Segmentation")
    logger.info("Using settings: " + str(settings))

    """
    Initialisation - Read in atlases
    - image files
    - structure files

        Atlas structure:
        'ID': 'Original': 'CT Image'    : sitk.Image
                          'Struct A'    : sitk.Image
                          'Struct B'    : sitk.Image
              'RIR'     : 'CT Image'    : sitk.Image
                          'Transform'   : transform parameter map
                          'Struct A'    : sitk.Image
                          'Struct B'    : sitk.Image
              'DIR'     : 'CT Image'    : sitk.Image
                          'Transform'   : displacement field transform
                          'Weight Map'  : sitk.Image
                          'Struct A'    : sitk.Image
                          'Struct B'    : sitk.Image

    """

    logger.info("")
    # Settings
    atlas_path = settings["atlasSettings"]["atlasPath"]
    atlas_id_list = settings["atlasSettings"]["atlasIdList"]
    atlas_structures = settings["atlasSettings"]["atlasStructures"]

    atlas_set = {}
    for atlas_id in atlas_id_list:
        atlas_set[atlas_id] = {}
        atlas_set[atlas_id]["Original"] = {}

        atlas_set[atlas_id]["Original"]["CT Image"] = sitk.ReadImage(
            "{0}/{1}/Images/{1}_CROP.nii.gz".format(atlas_path, atlas_id)
        )

        for struct in atlas_structures:
            atlas_set[atlas_id]["Original"][struct] = sitk.ReadImage(
                "{0}/{1}/Structures/{1}_{2}_CROP.nii.gz".format(
                    atlas_path, atlas_id, struct
                )
            )

    """
    Step 1 - Automatic cropping using a translation transform
    - Registration of atlas images (maximum 5)
    - Potential expansion of the bounding box to ensure entire volume of interest is enclosed
    - Target image is cropped
    """
    # Settings
    quick_reg_settings = {"shrinkFactors": [8, 2],
                          "smoothSigmas": [8, 2],
                          "samplingRate": 0.2
                         }

    registered_crop_images = []

    for atlas_id in atlas_id_list[:max([5, len(atlas_id_list)])]:
        # Register the atlases
        atlas_set[atlas_id]["RIR"] = {}
        atlas_image = atlas_set[atlas_id]["Original"]["CT Image"]

        reg_image, crop_tfm = initial_registration(
            image,
            atlas_image,
            moving_structure=False,
            fixed_structure=False,
            options=quick_reg_settings,
            trace=False,
            reg_method="Rigid",
        )

        registered_crop_images.append(reg_image)

    combined_image_extent = (sum(registered_crop_images) > 0)

    shape_filter = sitk.LabelShapeStatisticsImageFilter()
    shape_filter.Execute(combined_image_extent)
    bounding_box = np.array(shape_filter.GetBoundingBox(1))

    expansion = settings["autoCropSettings"]["expansion"]

    # Avoid starting outside the image
    crop_box_index = np.max([bounding_box[:3]-expansion, np.array([0,0,0])], axis=0)

    # Avoid ending outside the image
    crop_box_size = np.min([np.array(image.GetSize())-crop_box_index,  bounding_box[3:]+2*expansion], axis=0)

    crop_box_size = [int(i) for i in crop_box_size]
    crop_box_index = [int(i) for i in crop_box_index]

    img_crop = sitk.RegionOfInterest(image, size=crop_box_size, index=crop_box_index)

    #sitk.WriteImage(img_crop, "image_crop.nii.gz")

    """
    Step 2 - Rigid registration of target images
    - Individual atlas images are registered to the target
    - The transformation is used to propagate the labels onto the target
    """
    initial_reg = settings["rigidSettings"]["initialReg"]
    rigid_options = settings["rigidSettings"]["options"]
    trace = settings["rigidSettings"]["trace"]
    guide_structure = settings["rigidSettings"]["guideStructure"]

    for atlas_id in atlas_id_list:
        # Register the atlases
        atlas_set[atlas_id]["RIR"] = {}
        atlas_image = atlas_set[atlas_id]["Original"]["CT Image"]

        if guide_structure:
            atlas_struct = atlas_set[atlas_id]["Original"][guide_structure]
        else:
            atlas_struct = False

        rigid_image, initial_tfm = initial_registration(
            img_crop,
            atlas_image,
            moving_structure=atlas_struct,
            options=rigid_options,
            trace=trace,
            reg_method=initial_reg,
        )

        # Save in the atlas dict
        atlas_set[atlas_id]["RIR"]["CT Image"] = rigid_image
        atlas_set[atlas_id]["RIR"]["Transform"] = initial_tfm

        # sitk.WriteImage(rigidImage, f'./RR_{atlas_id}.nii.gz')

        for struct in atlas_structures:
            input_struct = atlas_set[atlas_id]["Original"][struct]
            atlas_set[atlas_id]["RIR"][struct] = transform_propagation(
                img_crop, input_struct, initial_tfm, structure=True, interp=sitk.sitkLinear
            )

    """
    Step 3 - Deformable image registration
    - Using Fast Symmetric Diffeomorphic Demons
    """
    # Settings
    resolution_staging = settings["deformableSettings"]["resolutionStaging"]
    iteration_staging = settings["deformableSettings"]["iterationStaging"]
    ncores = settings["deformableSettings"]["ncores"]
    trace = settings["deformableSettings"]["trace"]

    for atlas_id in atlas_id_list:
        # Register the atlases
        atlas_set[atlas_id]["DIR"] = {}
        atlas_image = atlas_set[atlas_id]["RIR"]["CT Image"]
        deform_image, deform_field = fast_symmetric_forces_demons_registration(
            img_crop,
            atlas_image,
            resolution_staging=resolution_staging,
            iteration_staging=iteration_staging,
            ncores=ncores,
            trace=trace,
        )

        # Save in the atlas dict
        atlas_set[atlas_id]["DIR"]["CT Image"] = deform_image
        atlas_set[atlas_id]["DIR"]["Transform"] = deform_field

        # sitk.WriteImage(deformImage, f'./DIR_{atlas_id}.nii.gz')

        for struct in atlas_structures:
            input_struct = atlas_set[atlas_id]["RIR"][struct]
            atlas_set[atlas_id]["DIR"][struct] = apply_field(
                input_struct, deform_field, structure=True, interp=sitk.sitkLinear
            )

    """
    Step 4 - Iterative atlas removal
    - This is an automatic process that will attempt to remove inconsistent atlases from the entire set

    Not included in this example - stay tuned for an update!
    Code is available in this repository, example:

    # Use simple GWV as this minises the potentially negative influence of mis-registered atlases
    for atlas_id in atlas_id_list:
        atlas_image = atlas_set[atlas_id]["DIR"]["CT Image"]
        weight_map = compute_weight_map(img_crop, atlas_image, vote_type='global')
        atlas_set[atlas_id]["DIR"]["Weight Map"] = weight_map

    reference_structure = settings["IARSettings"]["referenceStructure"]
    smooth_distance_maps = settings["IARSettings"]["smoothDistanceMaps"]
    smooth_sigma = settings["IARSettings"]["smoothSigma"]
    z_score_statistic = settings["IARSettings"]["zScoreStatistic"]
    outlier_method = settings["IARSettings"]["outlierMethod"]
    outlier_factor = settings["IARSettings"]["outlierFactor"]
    min_best_atlases = settings["IARSettings"]["minBestAtlases"]
    project_on_sphere = settings["IARSettings"]["project_on_sphere"]

    atlas_set = run_iar(
        atlas_set=atlas_set,
        structure_name=reference_structure,
        smooth_maps=smooth_distance_maps,
        smooth_sigma=smooth_sigma,
        z_score=z_score_statistic,
        outlier_method=outlier_method,
        min_best_atlases=min_best_atlases,
        n_factor=outlier_factor,
        iteration=0,
        single_step=False,
        project_on_sphere=project_on_sphere
    )

    """

    """
    Step 4 - Vessel Splining

    Not included in this example - stay tuned for an update!
    Code is available in this repository, example:

    vessel_name_list = settings["vesselSpliningSettings"]["vesselNameList"]
    vessel_radius_mm = settings["vesselSpliningSettings"]["vesselRadius_mm"]
    splining_direction = settings["vesselSpliningSettings"]["spliningDirection"]
    stop_condition = settings["vesselSpliningSettings"]["stopCondition"]
    stop_condition_value = settings["vesselSpliningSettings"]["stopConditionValue"]

    segmented_vessel_dict = vesselSplineGeneration(
        atlas_set,
        vessel_name_list,
        vesselRadiusDict=vessel_radius_mm,
        stopConditionTypeDict=stop_condition,
        stopConditionValueDict=stop_condition_value,
        scanDirectionDict=splining_direction,
    )

    """

    """
    Step 5 - Label Fusion
    """
    # Compute weight maps

    vote_type = settings['labelFusionSettings']['voteType']

    for atlas_id in list(atlas_set.keys()):
        atlas_image = atlas_set[atlas_id]["DIR"]["CT Image"]
        weight_map = compute_weight_map(img_crop, atlas_image, vote_type=vote_type)
        atlas_set[atlas_id]["DIR"]["Weight Map"] = weight_map

    combined_label_dict = combine_labels(atlas_set, atlas_structures)

    """
    Step 6 - Paste the cropped structure into the original image space
    """

    output_dir = settings['outputSettings']['outputDir']
    auto_structs = {}

    template_im = sitk.Cast((image * 0), sitk.sitkUInt8)
    vote_structures = settings["labelFusionSettings"]["optimalThreshold"].keys()

    for structure_name in vote_structures:
        optimal_threshold = settings["labelFusionSettings"]["optimalThreshold"][structure_name]
        binary_struct = process_probability_image(
            combined_label_dict[structure_name], optimal_threshold
        )
        paste_img = sitk.Paste(
            template_im, binary_struct, binary_struct.GetSize(), (0, 0, 0), crop_box_index
        )

        # Write the mask to a file in the working_dir
        mask_file = os.path.join(output_dir, f'{output_name}/{output_name}_{structure_name}_AUTO.nii.gz')

        try:
            os.mkdir(os.path.join(output_dir, f'{output_name}'))
        except:
            print('Output directory already exists')

        sitk.WriteImage(paste_img, mask_file)

        # Append the data to return
        auto_structs[structure_name] = paste_img

    return auto_structs


if __name__ == "__main__":
    test_id = 'Train-S1-009'
    test_image = sitk.ReadImage(f'./data/{test_id}/Images/{test_id}.nii.gz')

    output_name = test_id

    auto_structs = SimpleSeg(image=test_image, output_name=output_name, settings=settings)
