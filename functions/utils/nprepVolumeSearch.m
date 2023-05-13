function [gis, cent_index]=nprepVolumeSearch(R, dims, MV)
% prepare a volume search-light index set
%   parameters:
%       R  - range (geodesic radius) of search-light in voxel.
%       dims - [X Y Z] dimensions of target space.
%       MV - 3D array of masking volume, dimensions has to agree with dims;
%   return:
%       gis - cell array of vertex indices, every element of the cell
%               is a linear index (in TARGET GRID) array.
    lin_index = 1:prod(dims);
    cent_index = lin_index( MV>0 );%
    nvox = numel(cent_index);
    ind = zeros(dims);  % pseudo-volume for re-indexing
    ind(cent_index) = 1:nvox;
    [MX,MY,MZ]=ndgrid(1:(2*R+1),1:(2*R+1),1:(2*R+1));   % meshgrid mixes up X and Y
    LG = [MX(:) MY(:) MZ(:)]';
    c0 = [R+1; R+1; R+1];
    center = repmat(c0, [1 size(LG,2)]);
    R2 = R * R; 
    rid = (sum(LG.*LG + center.*center) - 2*sum(LG.*center))<R2;
    MX = MX(rid); MY = MY(rid); MZ = MZ(rid); 
    
    gis = cell(nvox, 1);
    MV = MV(:); 

    for cnt=1:nvox % Count through all included voxels
        SV = false(numel(MV),1);
        [xc,yc,zc]=ind2sub(dims,cent_index(cnt));
        mcx = MX + xc - c0(1);
        mcy = MY + yc - c0(2);
        mcz = MZ + zc - c0(3);
        mind=(((mcx>0)&mcx<=dims(1)) & ((mcy>0)&mcy<=dims(2)) & ((mcz>0)&mcz<=dims(3)));
        SV(sub2ind(dims, mcx(mind), mcy(mind), mcz(mind))) = true;
        % linear array of squared distance
        gis{cnt} = ind(SV & (MV>0)); % re-indexing, MV>0 ensures nonzero index
    end;
end