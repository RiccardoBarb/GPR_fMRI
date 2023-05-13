function bins = bin_angles(angles, n_bins);


% Returns bin number for a given angle (in degree) and given number of bins around a circle
%
% IN:   angles  = column vector of angels in degree
%       n_bins  = scalar, number of bins
%
% OUT:  bins    = column vector containing the assinged bin for every angle
%               from the input variable 'angles'

% CMS 102019    : get bin of single angle
% FMT 230320    : extended to allow binning of vector of angles


range = 360 / n_bins;

bins=NaN(1,size(angles,2));

for ang_=1:size(angles,2)
    
    angle=angles(ang_);
    
    for b = 1:n_bins
        bin_center = b*range - range;
        if isnan(angles(ang_))
            bins(ang_) = NaN;
        else
            if (angle - bin_center < range/2 & angle >= bin_center - range/2)...
                    | angle-bin_center >= (bin_center + 360) - (range/2)
                bins(ang_) = b;
            end
        end
    end
end