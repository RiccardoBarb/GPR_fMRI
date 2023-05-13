% This template script can be used to extract information (i.e. balanced accuracy 
% values) from a ROI and store them in a matrix. This enables the use of
% other software (i.e. JASP) to compute statistics.

% N.B. this script is meant as a tool to extract information from brain
% maps containing results of an analysis. We adopted this method to perform
% the statistical analyses illustrated in figure 6 of the manuscript. More 
% precisely:
% - for the analyses of results in the parietal areas we used SPM
%   Anatomy Toolbox to generate ROIs (see manusccripts for details). The
%   parietal ROIs are included in "data/fmri_data/ROIs/parietal".
% - For the analysis of the early visual cortex we combined the results of
%   the 2nd level localizer GLM specifying the contrast [1,1] (stimulus and
%   noise) and the anatomical areas V1, V2, V3 identified with SPM Anatomy
%   Toolbox). For details on how the rois have been combined see the
%   function "src/functions/roi_analyses/combine_roi.m" and
%   "src/functions/roi_analyses/convert_roi.m". The early visual ROI is
%   included in "data/fmri_data/ROIs/early_visual".
% - For the analysis of results in MT+ we used a ROI defined with the
%   procedure described in the script:
%   "src/other_scripts/create_single_subjects_roi_mt.m". please refer to
%   the script extract_data_from_individual_roi.m to analyze results of the
%   results from MT+.
%
clear all
SPM_path = ''
addpath(genpath(SPM_path))
results_dir = '../../data/fmri_data/gpr/';
%the path_to_roi specifies the location of the roi from which to extract
%the data to analyze - in the example we specified the early visual ROI.
path_to_roi = '../../data/fmri_data/ROIs/early_visual/ROI_v1v2v3_MNI_2mm_inter.nii'

%output_measure_tag refers to the name of the image from which we want to
%extract information. In this example we want to extract individual balanced accuracy
%from the outcome of the gpr reconstruction analysis after spatial
%smoothing and spatial normalization. N.B. the results of the GPR
%reconstruction analyses are not included in the data. To perform the
%analyses of the main results is necessary to run the GPR estimation pipeline 
% and the GPR reconstruction pipeline. 
output_measure_tag = 's6ws*4balanced*.img'
all_rois = dir(path_to_roi);
all_subjects = dir(strcat(results_dir,'s3*'));
current_roi = strcat(all_rois.folder,'/',all_rois.name);
fprintf('%s\n',current_roi)
roi_vol = spm_read_vols(spm_vol(current_roi));
for s_=1:length(all_subjects)
    current_sub = strcat(all_subjects(s_).folder,'/',all_subjects(s_).name);
    all_conditions = dir(strcat(current_sub,'/','*_*'));
    all_conditions(~[all_conditions.isdir])=[];
    for c_=1:length(all_conditions)
        current_condition = strcat(all_conditions(c_).folder,'/',all_conditions(c_).name);
        current_image = dir(strcat(current_condition,'/',output_measure_tag));
        data = spm_read_vols(spm_vol(strcat(current_image.folder,'/',current_image.name)));
        masked_data = data(logical(roi_vol));
        subject_condition_roi(s_,c_,:) = masked_data;
    end
end

% subject_condition_roi is a matrix of shape number of subjects x
% conditions x number of voxels in the roi
% the matrix contains the values of the measure to be considered for the
% analysis (i.e. reconstruction balanced accuracy above chance)
% condition 1 = high_report
% condition 2 = high_stimulus
% condition 3 = low_report
% condition 4 = low_stimulus
% condition 5 = mid_report
% condition 6 = mid_stimulus

%average across voxels of roi
mean_ROIs = mean(subject_condition_roi,3);
%rearrange the second dimensions from low_report to high_stimulus
mean_ROIs = [mean_ROIs(:,3),mean_ROIs(:,5),mean_ROIs(:,1),...
    mean_ROIs(:,4),mean_ROIs(:,6),mean_ROIs(:,2)];
