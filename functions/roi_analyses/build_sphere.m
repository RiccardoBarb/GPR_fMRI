function [SL_vol,vol_hdr] = build_sphere(volume,radius,center)

% Function to produce a sphere of a certain radius in a given volume
% centered at a specific location. The sphere definition is based on the
% euclidean distance between the center and the voxels in the volume

% INPUT: volume = the path of the volume in which to write the sphere
%        radius = expressed in mm
%        center = coordinates of the center of the sphere (in mm)

% OUTPUT: SL_vol = matrix with size = volume, with ones in voxels where
%                  euclidean distance is <= than radius.
%         vol_hdr = the hdr file of the volume used to write the sphere. It
%         can be used for saving the image file with spm_write_vol
%         (remember to change vol_hdr.fname before saving)
%        
% Written by Riccardo Barbieri

if size(center,2)>1
    center = center';
end

vol_hdr = spm_vol(volume);
[vol_img,vol_xyz] = spm_read_vols(vol_hdr);

SL_vol = zeros(size(vol_img));
vox_idx = find(vol_img);
coordinates = vol_xyz(:,vox_idx);

distance_map  = sqrt(sum((center - coordinates) .^ 2));

SL_coordinates = coordinates(:,distance_map<=radius);
SL_idx = vox_idx(distance_map<=radius);

SL_vol(SL_idx) = 1;