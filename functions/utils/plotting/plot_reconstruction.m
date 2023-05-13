function [dx] = plot_reconstruction(predicted_directions,true_directions,bins)
% The fuction calculates and plots  an histogram of the absolute angular 
% deviation between the predicted directions obtained from the function
% ME_cmGPR_pred and "true directions", namely the test labels of each 
% cv-fold. The input variables are organized in cells to deal with unbalanced
% number of trials acros cv-folds
% IN: predicted_directions - cell of nSL X ncv-folds containing ntrials
%     predictions (in rad)
%     true_directions - cell of size ncv-folds X 1 containing the test 
%     labels of each cv-fold.
%     bins - optional argument. Integer number representing the bin-size of
%     the histogram (in degrees)
% OUT: histogram of dx in deg;
%      dx: absolute angular deviation in rad
%
% Author: Riccardo Barbieri

%-------------------Calculate Angular Deviation---------------------------%
if nargin < 3
    bins = 10;
end
true = vertcat(true_directions{1:end});
dx = zeros(size(predicted_directions,1),1);

for i_=1:size(predicted_directions,1)
    pred =  vertcat(predicted_directions{i_,1:end});
    dx = pred - true;
    dxi = abs(dx)>pi;
    dx(dxi) = -1*sign(dx(dxi)) * 2*pi + dx(dxi);
end    
%---------------------------Plot histogram--------------------------------%
figure('Name','Stimulus reconstruction across SL','NumberTitle','off');
histogram(rad2deg(dx),-180:bins:180,'Normalization','Probability');
xlabel('Reconstruction deviation from true direction')
ylabel('Proportion of SL')
end