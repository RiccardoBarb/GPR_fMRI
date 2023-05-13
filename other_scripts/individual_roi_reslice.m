%this script must be used after inverse normalization individual_roi_inverse_normalization.m
%Once we have the mask in the individual anatomical space we reslice it into
%to the functoinal space.
% all_refs is the output of individual subjects accuracy map after
% smoothing.
% all_sources is the output of the script: individual_roi_inverse_normalization.m
%the results can be found in the folder data/fmri_data/xclass_individual_masks.
clear all
all_refs_files = ''
all_sources_files = ''
all_refs = dir(all_refs_dir);
all_sources = dir(all_sources_dir);

for i_=1:length(all_refs)
    
    ref = strcat(all_refs(i_).folder,'/',all_refs(i_).name);
    source = strcat(all_sources(i_).folder,'/',all_sources(i_).name);
    
    matlabbatch{1}.spm.spatial.coreg.write.ref = {ref};
    matlabbatch{1}.spm.spatial.coreg.write.source = {source};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch);
end
