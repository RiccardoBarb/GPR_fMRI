# GPR_fMRI

This code belongs to the paper *"Encoding of continuous perceptual choices in human early visual cortex"* by Barbieri, Töpfer et al. (2023), currently available as a preprint. It consists of a collection MATLAB scripts and functions for reproducing the main analyses described in the paper. The necessary data to be used for performing the analyses will soon be available on OSF. 

- Preprint: [Encoding of continuous perceptual choices in human early visual cortex](https://www.biorxiv.org/content/10.1101/2023.02.10.527876v)
- Data: TBA

## Requirements

This code was developed and run in [MATLAB R2019b](https://de.mathworks.com/help/matlab/release-notes.html) and [MATLAB R2021b](https://de.mathworks.com/help/matlab/release-notes.html). The analyses of data also requires additional Toolboxes:

- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- [GPRML](http://gaussianprocess.org/gpml/code/matlab/doc/)
- [Fieldtrip-Toolbox](https://www.fieldtriptoolbox.org/download/)

## Data

***TBA***

## Description

The [main](https://github.com/RiccardoBarb/GPR_fMRI/tree/main) directory includes the scripts for performing the main analyses, as well as a [tutorial](https://github.com/RiccardoBarb/GPR_fMRI/blob/main/tutorial_gprML.m) (`tutorial_gprML.m`). 
The [/functions](https://github.com/RiccardoBarb/GPR_fMRI/tree/main/functions) directory includes the entire collection of functions used by the entry scripts, to facilitate the exploration of the codebase they have been grouped by analysis type.
The [/other_scripts](https://github.com/RiccardoBarb/GPR_fMRI/tree/main/other_scripts) directory contains additional scripts that have been used to produce intermediate results such as individual ROIs, the fMRI preprocessing pipeline and the eye-tracking analyses (note that some of this script *cannot be used with the available data*, however they were included in this repository for the sake of transparency).

The following section provides a brief overview of the main analysis scripts. Further details can be found in the scripts comments. 

***N.B. Make sure to download the data, install the required libraries and add their path to MATLAB before running the analyses.***

### Tutorial

The script `tutorial_gprML.m` provides an introduction to the methodology employed for the analysis of fMRI data through a simulation. This is a good starting point for getting familiar with the concepts of voxel-wise response profile estimation and searchlight reconstruction, as well as the functions included in [/functions/gpr_analysis](https://github.com/RiccardoBarb/GPR_fMRI/tree/main/functions/gpr_analysis). The user can experiment with different parameters and observe the corresponding results. No additional data are required for running the tutorial.

### Behavioral Data Analysis

The script `behavioral_data_analysis.m` can be used to reproduce **Figure 2 A, B** and **C** as well as for producing single-subjects scatterplots of responses and stimulus directions like those included in **Appendix 2** of the manuscript. The final part of the script includes the code to run checks on the amount of discarded trials.

### GPR Estimation Pipeline

The script `GPR_estimation_pipeline.m` performs the estimation of voxel-wise response profiles as described in the paper. The analysis is set-up to automatically run on all subjects and for all the 6 conditions (*high_stimulus, mid_stimulus, low_stimulus, high_report, mid_report, low_report*). 
***N.B. This analysis is time consuming (the exact timing depends on specific machine setup) and will produce large output files (between 4GB and 6GB), we recommend testing the analysis on a single subject and for a single condition before running it for the entire dataset.***

### GPR Reconstruction Pipeline

The script `GPR_reconstruction_pipeline.m` performs searchlight-based multivariate stimulus/report reconstruction based on the voxel response profiles estimated with the `GPR_estimation_pipeline.m` script. The method is described in details in the manuscript. The analysis is set-up to automatically run on all subjects and for all the 6 conditions (*high_stimulus, mid_stimulus, low_stimulus, high_report, mid_report, low_report*). 
***N.B. This analysis is time consuming (the exact timing depends on specific machine setup), we recommend testing the analysis on a single subject and for a single condition before running it for the entire dataset. Before running this script is necessary to  run the `GPR_estimation_pipeline.m` since the outputs of that analysis will be used to perform the reconstruction.***

### GPR Generalization Pipeline

The script `GPR_generalization_pipeline.m` uses the GPR Reconstruction Pipeline to predict the reported/stimulus directions based on models estimated in a different coherence-condition (model generalization and model consistency). The method is described in details in the manuscript. The analysis is set-up to automatically run on all subjects but the user is required to manually specify the conditions on which the generalization should be performed.  Since the analysis is constrained to a specific ROI, the computations required to perform this analysis will be faster than those of the standard GPR Reconstruction Pipeline. The final part of the script can be used to directly plot the results displayed on **Figure 8** of the manuscript by using the `xclass_results_figure_8.mat` file included in the dataset.
***N.B. Before running this script is necessary to run the `GPR_estimation_pipeline.m`  since the outputs of that analysis will be used to perform the reconstruction. Alternatively, the third part of the script can be run by using the appropriate file including ready-to-be-plotted results*** 

### Pipeline to Reproduce Figure 3

The script `pipeline_to_reproduce_figure_3.m` produces a plot of 16 voxel response profiles estimated with the `GPR_estimation_pipeline.m` and the corresponding beta estimate of the measured brain activity. The voxels are randomly selected from the highest informative searchlight of a single subject. The final part of the script can be used to simply plot the results displayed on figure 3 by using the `s318_tuning_curves_figure_3.mat` file included in the dataset.
***N.B. Before running this script is necessary to run the `GPR_estimation_pipeline.m`. Alternatively, the third part of the script can be run by using the appropriate file including ready-to-be-plotted results*** 