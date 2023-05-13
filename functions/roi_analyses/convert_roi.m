function convert_roi(source,rois)

% Function to convert anatomical roi into the source image voxel size

% INPUT: source = the path of the source roi (having the voxel size to use
%                 for the conversion)
%        rois = structure with the folder content of the roi/s to be
%                transformed into source voxel size 
%                i.e. dir('pathoftherois/*.nii')
% 
%
% Written by Riccardo Barbieri

vol_hdr = spm_vol(source);
vol_img = spm_read_vols(vol_hdr);




for i = 1:length(rois)
    
   c_roi = strcat(rois(i).folder,'/',rois(i).name); 
   resampled_vol = zeros(size(vol_img));
   roi_hdr = spm_vol(c_roi);
   [roi_img,roi_xyz] = spm_read_vols(roi_hdr);


   new_coordinates = vol_hdr.mat \ [roi_xyz; ones(1, size(roi_xyz,2))];

       
   roi_vox = find(roi_img);
       
   roi_vox_new = round(new_coordinates(1:3,roi_vox));
       
   indexindex = sub2ind(size(vol_img),roi_vox_new(1,:),roi_vox_new(2,:),roi_vox_new(3,:));
       
       
   indexindex = unique(indexindex);
   resampled_vol(indexindex) = 1;
   vol_hdr.fname = [roi_hdr.fname(1:end-4) '_2mm.img'];
   spm_write_vol(vol_hdr,resampled_vol);
                        
end

