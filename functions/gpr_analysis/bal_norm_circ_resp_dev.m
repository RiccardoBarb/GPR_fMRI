function BAcc = bal_norm_circ_resp_dev(xr, x, meth)
% _
% Balanced Normalized Circular Response Deviation
% FORMAT BAcc = bal_norm_circ_resp_dev(xr, x)
% 
%     xr   - an n x v matrix of reconstructions (n observations, v channels)
%     x    - an n x 1 vector of angles, in radians, i.e. 0 <= x_i <= 2 pi
%     meth - method to use for numerical integration ('trapz' or 'rect')
% 
%     BAcc - a  1 x v vector of bal. norm. circ. resp. dev. accuracy values
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 21/01/2020, 17:48
%  Last edit: 22/01/2020, 14:38


% quantify accuracy
%-------------------------------------------------------------------------%
xn      = ~isnan(x);
dx      = xr(xn,:) - repmat(x(xn),[1 size(xr,2)]);
dxi     = abs(dx) > pi;
dx(dxi) = -1*sign(dx(dxi)) * 2*pi + dx(dxi);

% compute accuracies
%-------------------------------------------------------------------------%
Accs = (pi - abs(dx))./pi;

% sort accuracies
%-------------------------------------------------------------------------%
[xs, is] = sort(x(xn));
Accs     = Accs(is,:);

% integrate numerically
%-------------------------------------------------------------------------%
if strcmp(meth,'trapz')
    xs   = [xs; (2*pi)+min(xs)];
    Accs = [Accs; Accs(1,:)];
    BAcc = 1/(2*pi) * trapz(xs, Accs, 1);
end;
if strcmp(meth,'rect')
    xs   = [max(xs)-(2*pi); xs; (2*pi)+min(xs)];
    ds   = diff(xs);
    ws   = (ds(1:end-1)/2 + ds(2:end)/2)./(2*pi);
    BAcc = sum(repmat(ws,[1 size(xr,2)]).*Accs);
end;    