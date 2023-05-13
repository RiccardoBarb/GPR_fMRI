function xr = ME_cmGPR_pred(Y, xp, mup, S2e)
% _
% Cyclic Multivariate Gaussian Process Regression: predict variables
% FORMAT mur = ME_cmGPR_pred(Y, mup, S2e)
% 
%     Y   - an n x v data matrix (n data points, v measurement channels)
%     xp  - an m x 1 vector of angles, the reconstruction grid
%     mup - an m x v matrix of predicted expectations
%     S2e - a  v x v estimated covariance matrix
% 
%     xr  - an n x 1 vector of reconstructed variables
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/05/2019, 16:10
%  Last edit: 20/08/2019, 17:00


% Get model dimensions
%-------------------------------------------------------------------------%
n = size(Y,1);                  % number of observations
v = size(Y,2);                  % number of measurements
m = numel(mup);                 % number of predictions

% Step 1: Calculate log-likelihood values
%-------------------------------------------------------------------------%
LL = logmvnpdf( Y', mup', S2e );

% Step 2: Perform maximum likelihood estimation
%-------------------------------------------------------------------------%
[LL_max, LL_ind] = max(LL,[],2);
 xr = xp(LL_ind);