#/usr/bin/env python

"""
Name:           convert-dicom.py
Author:         Rob Finnegan
Date:           August 2020
Description:
    This script demonstrates a systematic approach to converting DICOM files
    It uses code in this repository, which can also be run from the command line.
    Run "python simpleseg_directory/code/rtstruct_to_nifti/convert.py --help" for more info.
"""

import os
from simpleseg.code.rtstruct_to_nifti import convert

#convert_rtstruct(dcm_img, dcm_rt_file, prefix=prefix, output_dir=output_dir, output_img=output_img, spacing=spacing)

# Step 1: Definitions

input_directory = "../data/DICOM/NSCLC-Radiomics"
output_directory = "../data/NIFTI_CONVERTED"

# Step 2: Search directory for data

try:
    os.mkdir(output_directory)
except:
    None

for input_file in os.listdir(input_directory):
    case_number = input_file[-3:]
    output_name = f"Study_{case_number}"
    output_subdirectory = '/'.join([output_directory, output_name])

    print(f"Running conversion on dataset: {case_number}")

    study_directory = '/'.join([input_directory, input_file])
    
    for study_file in os.listdir(study_directory):
        study_location = '/'.join([study_directory, study_file])
        dcm_dir = os.listdir(study_location)

        if 'CTLUNG' not in study_file:
           dcm_filename = '/'.join([study_location, dcm_dir[0]])

        if 'CTLUNG' in study_file:
            for dcm_rt_test in dcm_dir:
                if 'Segmentation' not in dcm_rt_test:
                    dcm_rt_filename = '/'.join([study_location, dcm_rt_test, '000000.dcm'])

    try:
        os.mkdir(output_subdirectory)
    except:
        None

    convert.convert_rtstruct(dcm_filename, dcm_rt_filename, prefix=output_name+'_', output_dir=output_subdirectory, output_img=output_name)

    all_struct = sum(structs.values())
    all_struct = all_struct > 0 # In case there are any overlaps in the lungs

    shape_filter = sitk.LabelShapeStatisticsImageFilter()
    shape_filter.Execute(all_struct)
    bounding_box = np.array(shape_filter.GetBoundingBox(1))
