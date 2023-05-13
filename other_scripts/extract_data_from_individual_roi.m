% This template script can be used to extract information (i.e. balanced accuracy 
% values) from an individual ROI and store them in a matrix. This enables 
% the use of other software (i.e. JASP) to compute statistics.

% N.B. this script is meant as a tool to extract information from individual
% brain maps containing results of an analysis. We adopted this method to perform
% the statistical analyses illustrated in figure 6 of the manuscript (MT+). 
% and the exploratory analysis illustrated in figure 7. 
% More precisely:
% - For the analysis of results in MT+ we used a ROI defined with the
%   procedure described in the script:
%   "src/other_scripts/create_single_subjects_roi_mt.m". MT ROIs are
%   included in "data/fmri_data/ROIs/mt_individual_masks".
% - For the exploratory analysis illustrated in figure 7 of the manuscript
%   the criteria to define ROIs are described in the methods section of the
%   paper. Other details are illustrated in the script
%   "src/other_scripts/individual_roi_inverse_normalization.m". These ROIs
%   are included in "data/fmri_data/ROIs/xclass_individual_masks".  

clear all
SPM_path = ''
addpath(genpath(SPM_path))
results_dir = '../../data/fmri_data/gpr/';
%all_rois specifies the location of the individual rois from which to extract
%the data to analyze - in the example we specified the individual MT+ ROIs.
all_rois = dir('../../data/fmri_data/ROIs/mt_individual_masks/s*/bilateral_localizer_MT.nii');
all_subjects = dir('../../data/fmri_data/gpr/s3*');
%output_measure_tag refers to the name of the image from which we want to
%extract information. In this example we want to extract individual balanced accuracy
%from the outcome of the gpr reconstruction analysis after spatial
%smoothing. N.B. the results of the GPR reconstruction analyses are not
% included in the data. To perform the analyses of the main results is 
% necessary to run the GPR estimation pipeline and the GPR reconstruction pipeline. 
output_measure_tag = 's6s*4balanced*.img';

for s_=1:length(all_subjects)
    fprintf('%d\n',s_)
    current_roi = strcat(all_rois(s_).folder,'/',all_rois(s_).name);
    
    roi_vol = spm_read_vols(spm_vol(current_roi));
    
    current_sub = strcat(all_subjects(s_).folder,'/',all_subjects(s_).name);
    all_conditions = dir(strcat(current_sub,'/','*_*'));
    for c_=1:length(all_conditions)
        current_condition = strcat(all_conditions(c_).folder,'/',all_conditions(c_).name);
        current_image = dir(strcat(current_condition,'/',output_measure_tag));
        data = spm_read_vols(spm_vol(strcat(current_image.folder,'/',current_image.name)));
        masked_data = data(logical(roi_vol));
        subject_condition_roi(s_,c_,:) = mean(masked_data);
    end
end

% subject_condition_roi is a matrix of shape number of subjects x
% condition.
% the matrix contains the values of the measure to be considered for the
% analysis (i.e. reconstruction balanced accuracy above chance)
% condition 1 = high_report
% condition 2 = high_stimulus
% condition 3 = low_report
% condition 4 = low_stimulus
% condition 5 = mid_report
% condition 6 = mid_stimulus

%rearrange the second dimensions from low_report to high_stimulus
subject_condition_roi = [subject_condition_roi(:,3),subject_condition_roi(:,5),...
    subject_condition_roi(:,1),subject_condition_roi(:,4),subject_condition_roi(:,6),...
    subject_condition_roi(:,2)];