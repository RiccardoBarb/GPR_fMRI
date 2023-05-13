function [mup, mux, s2x, hyp] = ME_cmGPR_mean(Y, x, xp)
% _
% Cyclic Multivariate Gaussian Process Regression: voxel-wise means
% FORMAT [mup, mux, s2x] = ME_cmGPR_mean(Y, x, xp, tol)
% 
%     Y   - an n x v data matrix (n data points, v measurement channels)
%     x   - an n x 1 vector of angles, in radians, i.e. 0 <= x_i <= 2 pi
%     xp  - an m x 1 vector of angles, the reconstruction grid
% 
%     mup - an m x v matrix of predicted expectations
%     mux - an n x v matrix of fitted expectations
%     s2x - a  1 x v vector of channel variances
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/03/2019, 17:40
%  Last edit: 21/01/2020, 17:28


% Set reconstruction grid, if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(xp)
    xp = [1:1:360].*(pi/180);
end;

% Get model dimensions
%-------------------------------------------------------------------------%
xn = ~isnan(x);                 % indices of non-NaNs
n  = size(Y,1);                 % number of observations
v  = size(Y,2);                 % number of measurements
m  = numel(xp);                 % number of predictions

% Specify mean, cov, lik for GPML
%-------------------------------------------------------------------------%
objfunc  = @gp;                 % objective function (to be minimized)
meanfunc = [];                  % assumed mean function
covfunc  = @covPeriodic_RDM;    % assumed covariance function
likfunc  = @likGauss;           % assumed likelihood function
infmeth  = @infGaussLik;        % prefered inference method

% Step 1: Estimate hyperparameters of the Gaussian process
%-------------------------------------------------------------------------%
hp     = struct('mean', [], 'cov', [0], 'lik', -1);
hyp(v) = minimize(hp, objfunc, -100, infmeth, meanfunc, covfunc, likfunc, x(xn), Y(xn,v));
for i = 1:v
    hyp(i) = minimize(hp, objfunc, -100, infmeth, meanfunc, covfunc, likfunc, x(xn), Y(xn,i));
end;

% Step 2: Generate likelihood function on continuous variable
%-------------------------------------------------------------------------%
mup = zeros(m,v);
mux = NaN(n,v);
s2x = NaN(n,v);

for i = 1:v
    [mux(xn,i), s2x(xn,i)] = gp(hyp(i), infmeth, meanfunc, covfunc, likfunc, x(xn), Y(xn,i), x(xn));
    [mup(:,i)            ] = gp(hyp(i), infmeth, meanfunc, covfunc, likfunc, x(xn), Y(xn,i), xp);
end;

s2x = nanmean(s2x,1);