function [Prec, Precs] = avg_norm_circ_resp_dev_simulation(xr, x)
% _
% Average Normalized Circular Response Deviation
% FORMAT [Prec, Precs] = avg_norm_circ_resp_dev(xr, x)
% 
%     xr   - an n x v matrix of reconstructions (n observations, v channels)
%     x    - an n x 1 vector of angles, in radians, i.e. 0 <= x_i <= 2 pi
% 
%     Prec  - a  1 x v vector of avg. norm. circ. resp. dev. precision values
%     Precs - an n x v matrix of trial-wise (not averaged) precision values
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 26/08/2019, 16:45
%  Last edit: 04/03/2020, 12:22


% quantify precision
%-------------------------------------------------------------------------%
dx      = xr - repmat(x,[1 size(xr,2)]);
dxi     = abs(dx) > pi;
dx(dxi) = -1*sign(dx(dxi)) * 2*pi + dx(dxi);

% compute precision
%-------------------------------------------------------------------------%
Precs = (pi - abs(dx))./pi;
Prec  = nanmean(Precs, 1);