clear all
% the following script has been used to inverse normalize rois into individual
% subject space.
% all_fields are fieldmaps from structural scans and all_spaces are the
% individual sbjects structural scans. The output of this script will be in
% 0.8x0.8x0.8 mm and need to be resliced (see individual_roi_reslice.m).
% roi is the image resulting from the second level analysis on the GPR
% accuracy maps (whole-brain searchlight analysis). More specifically, it 
% is the intersection of voxels obtained by the contrast [1,1,1] defined 
% on two one-way factorial designs (one for stimulus and one for report) 
% with coherence level as a within-subject factor.
% The ROI was subsequently projected in the individual subject
% voxel space.

all_fields_files = '';
all_spaces_files = '';
roi_file = '';
all_fields = dir(all_fields_files);
all_spaces = dir(all_spaces_files);
roi = dir(roi_file)
target = '../data/fmri_data/xclass_individual_masks/';

for i_=1:length(all_fields)
    splt = strsplit(all_fields(i_).folder,'/');
    pref = strcat(splt{7},'_');
    fds = strcat(all_fields(i_).folder,'/',all_fields(i_).name);
    spc = strcat(all_spaces(i_).folder,'/',all_spaces(i_).name);
    
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fds};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {spc};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {roi};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {target};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = pref;
    
    spm_jobman('run', matlabbatch);
end