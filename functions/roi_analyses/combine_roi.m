function [new_roi,vol_hdr] = combine_roi(roi_a,roi_b,how,save)

% Function to combine a reference roi with another roi 

% INPUT: roi_a = the path of the "reference" roi - e.g. anatomical image
%        roi_b = the path of the  roi to be combined with the "reference" -
%        e.g. mask of roi obtained from a functional localizer.
%
% N.B:   the order of roi_a and roi_b matter only for determining the 
%        destination of the resulting vol_hdr in case of saving
%
%        how = string specifying the type of combination of the rois
%             'inter' compute the intersection of roi_a and roi_b, 'union'
%              computes the union; default = 'inter'.
%        save = if 1 save the vol_hdr in a new folder located in the path 
%               of roi_b with suffix specified by newf.
%
% OUT: new_roi = the binary 3D matrix of the new volume.
%      vol_hdr = the image to be written in an img file with newname
%
% Written by Riccardo Barbieri

if isempty (how)
    how = 'inter';
end
if isempty(save)
    save = 0;
end
if strcmp(how,'inter')
   namepart = '_inter.img';
elseif strcmp(how,'union')
    namepart = '_union.img';
else
    fprintf('%s method not recognized\nplease use inter or union',how)
    return
end

vol_hdr = spm_vol(roi_b);
roi_b_img = spm_read_vols(vol_hdr);
roi_a_img = spm_read_vols(spm_vol(roi_a));


new_roi = zeros (size(roi_a_img));

if strcmp(namepart, '_inter.img')
    idcombined = roi_a_img & roi_b_img;
elseif strcmp(namepart,'_union.img')
    idcombined = roi_a_img | roi_b_img;
end

new_roi(idcombined) = 1;

oldname = strsplit(vol_hdr.fname,'/');
newf = strcat(oldname{end-1},'_loc_combined');
newname = strcat(oldname{end}(1:end-4),namepart);
oldname{end-1} = newf;
oldname{end} = newname;
complete_newname = strcat(oldname,'/');

vol_hdr.fname = strcat(complete_newname{1:end});
vol_hdr.fname = vol_hdr.fname(1:end-1);

if save
    if ~exist(newf)
        mkdir(newf)
    end
    
    spm_write_vol(vol_hdr,new_roi);
end

end