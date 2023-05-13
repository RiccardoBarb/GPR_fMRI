function S2e = ME_cmGPR_cov(Y, mux, s2x, reg)
% _
% Cyclic Multivariate Gaussian Process Regression: ROI-based covariance
% FORMAT S2e = ME_cmGPR_cov(Y, mux, s2x, reg)
% 
%     Y   - an n x v data matrix (n data points, v measurement channels)
%     mux - an n x v matrix of fitted expectations
%     s2x - a  1 x v vector of channel variances
%     reg - a string specifying the method of regularization
% 
%     S2e - a  v x v estimated (and regularized) covariance matrix
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/03/2019, 17:40
%  Last edit: 26/08/2019, 16:30


% Set variances, if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(s2x)
    s2x = var(Y,[],1);
end;

% Set regularization, if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(reg)
    reg = 'logistic';           % options: none, step, arctan, logistic
end;

% Get model dimensions
%-------------------------------------------------------------------------%
xn = ~isnan(mux(:,1));          % indices of non-NaNs
n  = size(Y,1);                 % number of observations
v  = size(Y,2);                 % number of measurements
n  = n - sum(~xn);              % minus number of NaNs

% Step 1: Calculate mixing coefficient
%-------------------------------------------------------------------------%

switch reg
    case 'none'                 % no regularization
        r = 0;
    case 'step'                 % step-wise regularization
        r = 0; if v >= n, r = 1/2; end;
    case 'arctan'               % using arcustangens function
        r = 1/pi * (atan(log(v/n)) + pi/2);
    case 'logistic'             % using logistic function
        r = 1/(1 + exp(-log(v/n)));
end;

% Step 2: Estimate across-channel covariance
%-------------------------------------------------------------------------%
R   = Y(xn,:) - mux(xn,:);
S2e = (1/n) * (R'*R);
S2e = (1-r) * S2e + r * diag(s2x);