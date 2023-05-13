%script to create MT ROI on single subjects in the dataset
clear all
spm_path = %path to spm
addpath(genpath(spm_path));
addpath(genpath(strcat(pwd,'/functions'))); 
% this file contains coordinates of the clusters activated by coherent 
% VS incoherent motion (localizer), and located lateral to the 
% parietal-occipital sulcus. The areas can be viewed examining individual
% 'localizer_fwe.nii' files and going to the specified coordinates.
%
load('../../data/fmri_data/ROIs/mt_individual_masks/MT_coordinates.mat');
mask_files = dir('../../data/fmri_data/ROIs/mt_individual_masks/s*/mask.nii');
radius = 10;%radius of the sphere to center at the chosen coordinate

for i_ = 1:length(MT_coordinates)
    
    current_mask = fullfile(mask_files(i_).folder,mask_files(i_).name);
    current_localizer = fullfile(mask_files(i_).folder,'localizer_fwe.nii');
    
    vol_localizer=spm_read_vols(spm_vol(current_localizer));
    
    center_right = MT_coordinates(i_).right;
    center_left = MT_coordinates(i_).left;
    
    filename_combined = ('bilateral_localizer_MT.nii');
    % create sphere on right and left mt 
    [vol_r,~] = build_sphere(current_mask,radius,center_right);
    [vol_l,hdr] = build_sphere(current_mask,radius,center_left);
    
    hdr.fname = strcat(hdr.fname(1:end-8),filename_combined)
    % combine them in a single volume
    vol_bilateral = vol_r | vol_l;
    % intersection with localizer
    vol_bilateral_localizer = vol_bilateral & vol_localizer;
    
    spm_write_vol(hdr,vol_bilateral_localizer);
    
    vol_r=[];vol_l=[];vol_bilateral=[];vol_bilateral_localizer=[];
end
    
    
    