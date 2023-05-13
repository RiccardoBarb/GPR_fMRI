function [center_xyz] = get_roi_center(roi)

% Function to obtain the center of a roi in 3d coordinates

% INPUT: roi = the path of the roi volume
% OUTPUT: center_xyz = coordinates of the center of the roi (in mm)
%
% Written by Riccardo Barbieri

vol_hdr = spm_vol(roi);
[vol_img,vol_xyz] = spm_read_vols(vol_hdr);

vox_idx = find(vol_img);
coordinates = vol_xyz(:,vox_idx);

max_xyz = max(coordinates,[],2);
min_xyz = min(coordinates,[],2);
center_xyz = (max_xyz-min_xyz)./2+min_xyz; 


end