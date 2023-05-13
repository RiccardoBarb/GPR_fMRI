function L = logmvnpdf(X, Mu, Sigma)
% _
% Log Probability Density Function of the Multivariate Normal Distribution
% FORMAT L = logmvnpdf(X, Mu, Sigma)
% 
%     X     - an k x m matrix, the m values at which the PDF is evaluated
%     Mu    - an k x n matrix, the n means of the multivariate normal distribution
%     Sigma - an k x k matrix, the 1 covariance martrix of the multivariate normal
% 
%     L     - an m x n matrix, the log-PDF values of the m values at the n means
% 
% FORMAT L = logmvnpdf(x, mu, Sigma) calculates the logarithmized probability
% density function of the multivariate normal distribution with means Mu and
% covariance Sigma at the values X [1].
% 
% This function allows to calculate the log-PDF at *several* values x,
% using *several* means mu, but always the same covariance matrix Sigma.
% 
% References:
% [1] Wikipedia: "Multivariate normal distribution";
%     URL: https://en.wikipedia.org/wiki/Multivariate_normal_distribution.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 11/04/2019, 09:10
%  Last edit: 21/05/2019, 12:45


% Get matrix dimensions
%-------------------------------------------------------------------------%
m = size(X,2);
n = size(Mu,2);
k = size(Sigma,2);

% Compute constant part
%-------------------------------------------------------------------------%
S_ld  = logdet(Sigma);
S_inv = inv(Sigma);
L_con = - (k/2) * log(2*pi) - (1/2) * S_ld;
  
% Compute variable part
%-------------------------------------------------------------------------%
L_var = zeros(m,n);
for i = 1:m
    for j = 1:n
        L_var(i,j) = -(1/2) * (X(:,i)-Mu(:,j))' * S_inv * (X(:,i)-Mu(:,j));
    end;
end;

% Compute log-PDF
%-------------------------------------------------------------------------%
L = L_con + L_var;