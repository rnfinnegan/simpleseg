# SimpleSeg

An open-source, Python-based atlas-based segmentation framework. If you use this code in a publication please cite [this paper](https://doi.org/10.1088/1361-6560/ab0ea6):


> Finnegan, R., Dowling, J., Koh, E.-S., Tang, S., Otton, J., Delaney, G., Batumalai, V., Luo, C., Atluri, P., Satchithanandha, A., Thwaites, D., Holloway, L. (2019). Feasibility of multi-atlas cardiac segmentation from thoracic planning CT in a probabilistic framework. Physics in Medicine & Biology, 64(8), 085006. https://doi.org/10.1088/1361-6560/ab0ea6



## Design Aims

This software is designed with flexibility in mind. To this end, the code is quite generalisable and modular. We encourage using this code to build your own image segmentation pipelines, and have included a demonstration of how this code could be run for simple (single script!) image segmentation.

In particular, all intermediate data objects (deformation vector fields, probabilistic labels, etc.) are available. These contain a wealth of information, and can be used in a wide variety of ways.

We would love to hear how you have used this code, as well as any suggestions, comments or feedback!

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. This was originally developed on Linux (Ubuntu), although getting it to work on Windows/Mac should be fairly straightforward.

The deformable image registration is multi-threaded and can utilise as many CPU cores as you have available. Currently, GPU-acceleration is under development.

### Prerequisites

Before starting, make sure you have Python (version 3) and a Python package manager (e.g. Pip, Conda). Next, you can install all of the required packages using the included list of Python packages.

```
pip install -r requirements.txt
```

### Included dataset

We have included a dataset in this repository, accessed from The Cancer Imaging Archive:

```
Aerts, H. J. W. L., Wee, L., Rios Velazquez, E., Leijenaar, R. T. H., Parmar, C., Grossmann, P., ... Lambin, P. (2019). Data From NSCLC-Radiomics [Data set]. The Cancer Imaging Archive. https://doi.org/10.7937/K9/TCIA.2015.PF0M9REI
```

This dataset consists of non-contrast CT imaging of 60 lung cancer patients, with five organ-at-risk contours (heart, left lung, right lung, spinal cord, esophagus). We have only included 36 of these in this repository (from the *training* subset).

This collection may not be used for commercial purposes. This collection is freely available to browse, download, as outlined in the Attribution-NonCommercial 3.0 Unported (CC BY-NC 3.0) https://creativecommons.org/licenses/by-nc/3.0/. It is able to be shared and adapted. We have converted the original DICOM images and RT-STRUCT files into compressed NifTI files. We have also cropped the images to contain to the extent of the contours (to make deformable image registration more efficient).


### Demonstration

The SimpleSeg codebase is included (in the *code* directory), and a demonstration of how this code can be used for end-to-end segmentation is available in the *SimpleSeg.py* script.

Have a read of the *SimpleSeg.py* script, most of the settings are made to be user-understandable. At the bottom of the script, the target image is defined.

The script can be run in a terminal:

```
python SimpleSeg.py
```

### Measuring segmentation accuracy

We include a simple script to generate quantitative label similarity metrics, this is useful to evaluate the performance of an automatic segmentation algorithm. To run this script, in a terminal:

```
python ComputeMetrics.py
```

This will create a CSV file containing results.

# Additional information

## Authors

* **Robert Finnegan** - robert *dot* finnegan *at* sydney.edu.au
* **Philip Chlap** - philip *dot* chlap *at* unsw.edu.au

## Contributing

You are welcome to contribute to this project!

## Versioning

This software is currently under development, and is subject to change that may cause compatibility issues.

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE, see the [LICENSE](LICENSE) file.

## Acknowledgments

* [SimpleITK](http://www.simpleitk.org) - A fantastic abstraction layer/wrapper for [ITK](http://www.itk.org). This is used for much of the image processing and registration.
