% GPR MLE meta script
% _
% This script investigates the distribution of average and balanced
% precisions following reconstruction using Gaussian process regression
% (GPR) and maximum likelihood estimation (MLE).
% 
% The script can be used to reproduce the analyses described in the
% appendix 2 of the manuscript.
%
% If you want to plot results included in the data/simulation folder you
% can direclty run the code in section 3 and 4 after loading the matfile
% you want to plot.
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
%   Date: 04/03/2020, 15:30


clear
close all
GPR_ML_path = %path to GPR_ML_library http://gaussianprocess.org/gpml/code/matlab/doc/
addpath(genpath(GPR_ML_path))
addpath(genpath('../../functions'))
save_results = false;

%%% Step 0: Specify generative model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global settings
there_is_an_effect  = false;
angles_are_balanced = false;
rng(5);

if there_is_an_effect
    part1 = 'effect_'
else
    part1 = 'no_effect_'
end
if angles_are_balanced
    part2 = 'balanced'
else
    part2 = 'imbalanced'
end

% set parameters
v = 1000;                       % number of simulations (voxels)
r = 120;                        % searchlight radius (in voxels)
S = 10;                         % number of sessions
t = 16;                         % trials per session
n = S*t;                        % number of samples

% independent variables
K  = 1;                         % number of directions
mu = [0:((2/K)*pi):(2*pi-0.1)]; % prefered directions
ka = (2^K)*ones(size(mu));      % preference precisions
la = 1/K*ones(1,K);             % preference frequencies
if K == 2, la = [0.25, 0.75]; end;

% voxel tuning functions
L  = gamrnd(2,2,[1 v]);         % tuning smoothness
s2 = chi2rnd(4,[1 v]);          % noise variances

% voxel covariance matrix
ny = 0.25;                      % space constant
Sy = toeplitz(ny.^[0:1:(v-1)]); % spatial correlations

% specify reconstruction precision
prec = 1;                       % prediction precision


%%% Step 1: Simulate some data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample independent variables
if angles_are_balanced
    x = MD_unirnd(0, 2*pi, n);
else
    x = [];
    for k = 1:K
        x = [x; MD_vmrnd(mu(k), ka(k), la(k)*n)];
    end;
    x(x<0) = x(x<0) + 2*pi;
end;

% sample ground truth values
x0 = [(prec/180)*pi:(prec/180)*pi:2*pi]';
Yt = zeros(numel(x0),v);
for i = 1:v
    K0 = ME_GPR_kern(x0, x0, 'per', L(i));
    Yt(:,i) = mvnrnd(zeros(size(x0)), K0);
end;

% sample dependent variables
Y = zeros(n,v);
for i = 1:v
    for j = 1:n
        Y(j,i) = Yt(ceil((x(j)/(2*pi))*(360/prec)),i);
    end;
end;
V = eye(n);
E = matnrnd( zeros(n,v), V, diag(sqrt(s2))*Sy*diag(sqrt(s2)) );

% generate observed data
if there_is_an_effect
    Y = Y + E;
else
    Y = E;
end;

% generate searchlights
SL   = zeros(v,2*r+1);
VpSL = zeros(1,v);
for i = 1:v
    vox_ind = [max([i-r, 1]):1:min([v, i+r])];
    VpSL(i) = numel(vox_ind);
    SL(i,1:VpSL(i)) = vox_ind;
end;
clear vox_ind


%%% 2. Perform data analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get dimensions
t = size(Y,1);
v = size(Y,2);
l = size(SL,1);
tpS = t/S;                      % trials per data subset
t1  = (S-1)*tpS;                % number of in-sample trials

% reconstruction grid
xp = [prec:prec:360]'.*(pi/180);% build reconstruction grid
m  = numel(xp);                 % number of candidate angles

% preallocate results
mup = zeros(m,v,S);             % predicted expectations
mux = zeros(t1,v,S);            % fitted expectations
s2x = zeros(S,v);               % channel variances
xr  = zeros(t,l);               % reconstructed angles

% h  indexes CV fold
% i* indexes trials
% j* indexes voxels
% k  indexes searchlight

fprintf('\n-> Voxel-wise GPR:\n');

% Phase 1: estimation
for h = 1:S
    
    fprintf('   - CV fold %02d ... ', h);
    
    % get (out-of-sample) test data
    i2 = [(h-1)*tpS+1:h*tpS];
    Y2 = Y(i2,:);
    x2 = x(i2);
    
    % get (in-sample) training data
    i1 = setdiff([1:t],i2);
    Y1 = Y(i1,:);
    x1 = x(i1);
    
    %%% Step 1: voxel-wise GPR; all voxels are analyzed together %%%%%%%%%%
    [mup(:,:,h), mux(:,:,h), s2x(h,:)] = ME_cmGPR_mean(Y1, x1, xp);
    
    fprintf('successful!\n');
    
end;

fprintf('\n-> Searchlight-based MLE:\n');

% Phase 2: reconstruction
for h = 1:S
    
    fprintf('   - CV fold %02d ... ', h);
    
    % get (out-of-sample) test data
    i2 = [(h-1)*tpS+1:h*tpS];
    Y2 = Y(i2,:);
    x2 = x(i2);
    
    % get (in-sample) training data
    i1 = setdiff([1:t],i2);
    Y1 = Y(i1,:);
    x1 = x(i1);
    
    %%% Step 2: ROI-based MLE; loop over searchlights required %%%%%%%%%%%%
    for k = 1:l
        % get searchlight voxels
        jk  = SL(k,SL(k,:)~=0);
        % estimate covariance in searchlight
        S2e = ME_cmGPR_cov(Y1(:,jk), mux(:,jk,h), s2x(h,jk), 'logistic');
        % reconstruct angles from searchlight
        xr(i2,k) = ME_cmGPR_pred(Y2(:,jk), xp, mup(:,jk,h), S2e);
    end;
    
    fprintf('successful!\n');
    
end;

fprintf('\n');

% save results
if save_results
    result_name = ['../../../data/simulation/',...
        part1,part2,'.mat']
    save(result_name,'mup','mux','s2x','xr','x','L','s2')
end

%% %%% 3. Analyze results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% N.B. if you simply want to plot the results already included in the 
%% dataset you can load the mat files located in /data/simulation and start
%% from this point

% quantify precision
[Prec, Precs] = avg_norm_circ_resp_dev_simulation(xr, x);

[BPrec]   = bal_norm_circ_resp_dev_simulation(xr, x,'trapz');

%    xr is a t x l matrix of reconstructed angles,
% Precs is a t x l matrix of reconstruction precisions,
%  Prec is a 1 x l vector of reconstruction precisions,
% BPrec is a 1 x l vector of balanced precisions and

l= size(xr,2);

% calculate histograms (parameters)
dx = 0.5;
xs2= (dx/2):dx:max(s2);
xL = (dx/2):dx:max(L);
ns2= hist(s2, xs2);
nL = hist( L, xL);

% calculate histograms (angles)
dx = pi/12;
xh = (dx/2):dx:(2*pi-dx/2);
nt = hist(x, xh);
nr = hist(xr(:), xh);



% calculate histograms (precisions)
dx = 5;
xp = 0:dx:100;
np = hist(Precs(:)*100, xp);

% calculate histograms (overall precisions)
nap= hist(Prec*100, xp);
nbp= hist(BPrec*100, xp);

% process precisions
Precs_mean = mean(Precs*100,2);
Precs_SE   = std(Precs*100,[],2)./sqrt(l);

% calculate mean/SEM
Prec_mean  = mean(Prec*100);
BPrec_mean = mean(BPrec*100);
Prec_var   = var(Prec*100);
BPrec_var  = var(BPrec*100);
Prec_std   = std(Prec*100);
BPrec_std  = std(BPrec*100);
Prec_SE    = std(Prec*100)/sqrt(l);
BPrec_SE   = std(BPrec*100)/sqrt(l);


%% %%% 4. Plot results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot overall precisions
figure('Name', 'GPR-MLE-circ-per: precisions', 'Color', [1 1 1], 'Position', [1280 -160 1600 900]);

subplot(1,2,1); hold on;
plot(Prec*100, BPrec*100, 'ob', 'LineWidth', 2, 'MarkerSize', 10);
plot(Prec_mean, BPrec_mean, 'ok', 'MarkerFaceColor', [1 1 0], 'LineWidth', 2, 'MarkerSize', 10);
plot([0 100], [0 100], ':k', 'LineWidth', 2);
plot([50 50], [0 100], '--k', 'LineWidth', 2);
plot([0 100], [50 50], '--k', 'LineWidth', 2);
axis([0 100 0 100]);
axis square;
set(gca,'Box','On');
legend('precisions', 'mean', 'Location', 'SouthEast');
xlabel('average precision', 'FontSize', 12);
ylabel('balanced precision', 'FontSize', 12);
title('Precision of angle reconstructions', 'FontSize', 16);

subplot(2,2,2); hold on;
bar(xp, nap, 'FaceColor', [1 0 0]);
plot([Prec_mean, Prec_mean], [0, (11/10)*max(nap)], '-k', 'LineWidth', 2);
axis([0, 100, 0, (11/10)*max(nap)]);
set(gca,'Box','On');
xlabel('average precision', 'FontSize', 12);
ylabel('number of searchlights', 'FontSize', 12);
title(sprintf('Distribution of average precision (mean = %1.2f, std = %1.4f)', Prec_mean, Prec_std), 'FontSize', 16);

subplot(2,2,4); hold on;
bar(xp, nbp, 'FaceColor', [0 1 0]);
plot([BPrec_mean, BPrec_mean], [0, (11/10)*max(nbp)], '-k', 'LineWidth', 2);
axis([0, 100, 0, (11/10)*max(nbp)]);
set(gca,'Box','On');
xlabel('balanced precision', 'FontSize', 12);
ylabel('number of searchlights', 'FontSize', 12);
title(sprintf('Distribution of balanced precision (mean = %1.2f, std = %1.4f)', BPrec_mean, BPrec_std), 'FontSize', 16);