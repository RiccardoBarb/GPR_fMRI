function predict_gprML(scan_dir,subjectscan,target_dir,sl_radius,ffs,c1,c2)
% this function start the GPR prediction on preprocessed beta-images of
% fMRI data. The prediction is performed by applying multivariate MLE on
% 1) the continuous activity pattern (mup) and
% 2) the estimated spatial covariance (S2_est) in a given searchlight.

% The voxel-wise estimations of mup mux s2 hyp, need to be performed with the
% function estimate_gprML

% The analysis is written to run on a server with 24 workers. We
% parallelize the cross-validation for loop to speed up the analysis.


%--------------------% set-up parallel pool %-----------------------------%
Nworkers = 24;

p=(gcp('nocreate'));

if isempty(p)
    parpool(Nworkers);
    maxNumCompThreads(Nworkers);
else
    fprintf('Already connected to %d Workers \n close and restart them \n',p.NumWorkers);
    delete(p);                                                             %felix 210819 delete pools to try if they then release their momory
    parpool(Nworkers);                                                     %felix 210819 restart pools
    maxNumCompThreads(Nworkers);
end

%turn warnings off
pctRunOnAll warning off
%-------------------------------------------------------------------------%
%% Part 1: SL creation/loading
fprintf('Setting up GPR estimation for %s, on %s coherence for %s direction \n',subjectscan,c1,c2)

% get or create/save search-light indices.
% NB: during the GPR estimation phase we will use every voxel in the volume

id_sl_dir = strfind(target_dir,subjectscan);
sl_dir = strcat(target_dir(1:[id_sl_dir-1]),'SL');

if ~exist(sl_dir)
    mkdir(sl_dir)
end

maskmat = fullfile(sl_dir, [subjectscan '_volume_search_' num2str(sl_radius) '.mat']);
maskmap = fullfile(scan_dir, sprintf('mask.nii'));
betamap = fullfile(scan_dir, sprintf('beta_0001.nii'));

resultsvol_hdr=spm_vol(betamap);%Get headers of betas to write results

maskvol_hdr=spm_vol(maskmap);%Get headers of mask
mask_vol=spm_read_vols(maskvol_hdr);%%Get actual mask image
sz=size(mask_vol);
% if the searchlight indices don't exist we create them, otherwise we load
% them from the existing matfile
if ~exist(maskmat)
    
    fprintf('constructing volume searchlight indices...');
    
    % gid will contain all the indices for interested voxels
    [gis, gid] = nprepVolumeSearch(sl_radius, size(mask_vol), mask_vol);
    
    fprintf('\nvolume searchlight indices constructed\n');
    
    save(maskmat, 'gis', 'gid');
else
    load(maskmat);
end

clear mask_vol maskvol_hdr

nvox = numel(gid);

fprintf('Subject %s Total voxels in ROI: %d\n',subjectscan, nvox);


%% Part 2: Get Beta images estimated with 1st level-GLM
%%          and the voxel-wise estimated mup mux s2 hyp
gid = gid';

%get the number of the first and final scans
s1=ffs(1);
s2=ffs(2);

%total number of runs
nruns=s2-s1+1;

% Assign each coherence-condition to the relevant betas to be taken from
% the first level-GLM

if strcmpi(c1,'low'), coherence = 0;
elseif strcmpi(c1,'mid'), coherence = 16;
elseif strcmpi(c1,'high'), coherence = 32;end
if strcmpi(c2,'stimulus'), condition = 3;
elseif strcmpi(c2,'report'), condition = 4;end

betanumbers = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]+coherence;
ntrials = numel(betanumbers);
%--------------------% Get Betas and lables %-----------------------------%
fprintf('Extracting subject %s        \n\n',subjectscan);

[extracted_labels,vectors] = ...
    extract_glm_design(scan_dir,condition,ntrials,betanumbers,gid,nruns);

fprintf('\nDone!');
%---------------------% Get mup mux s2 hyp %------------------------------%
fprintf('\nloading estimated GPR hyperparameters and predictions...')

load (fullfile(target_dir,'prep_mu_periodic2b.mat'));

disp('\nDone!')
%------------------------% Analysis setup %-------------------------------%
%index for run-wise cv
ind = logical(repmat([0 ones(1, nruns-1)], [1, 1]))';

% reconstruction grid
prec = 1;                         % estimation precision in degrees
xp = [prec:prec:360]'.*(pi/180);  % build reconstruction grid
m  = numel(xp);                   % number of candidate angles

%preallocate outcome variables

acc_vol = zeros(nvox, 1);                 % SL cross-validated accuracy
bal_acc_vol = zeros(nvox, 1);             % SL balanced cross-validated accuracy
predicted_directions = cell(nvox,nruns);  % SL trial-wise predicted direction
true_directions = cell(nruns,1);          % tesing lables of each cv-fold

%if starting from day 2 itrain and itest don't work properly. Here we
%adjust s1 and s1 to match first (=1) run
s1=1;

%% Part 3: Perform searchlight-based cyclic GPR prediction

fprintf(1,'\n %s %s %s %s %s %s %s \n','GPR Prediction for',subjectscan ,c1,'coherence','for',c2,'direction');

for ii=1:nvox
    current_index = ii;
    fprintf('\n %0.1f%%\n',(ii/nvox).*100); % \b is backspace
    
    % Bring only data from the current searchlight into the parallel
    % for loop
    vv_sl = vectors(:,gis{ii},:);  % trial-wise betas
    hyp_sl = hyp(:,gis{ii});       % voxel-wise kernel hyperparameters
    mup_sl = mup(:,gis{ii},:);     % predicted continuous voxel response
    mux_sl = mux(:,gis{ii},:);     % fitted trial-wise voxel response
    s2x_sl = s2x(:,gis{ii});       % voxels variance
    n_sl_vox = numel(gis{ii});     % number of voxels in the sl
    parfor r=1:nruns
        
        %training and testing index for current cv fold
        itrain = find( [1:nruns]~=r-s1+1);
        itest  = r-s1+1;
        
        %vectorize training data for current cv fold
        vectors_train = vv_sl(:,:,itrain);
        train_samples = permute(vectors_train,[3,1,2]);
        train_samples = reshape(train_samples,[ntrials*(nruns-1),n_sl_vox]);
        
        labels_train=extracted_labels(:,itrain)';
        labels_train=labels_train(:)./360*2*pi;     %Labels are converted in rad
        
        %vectorize testing data for current cv fold
        
        test_samples = vv_sl(:,:,itest);
        labels_test=extracted_labels(:,itest);
        labels_test=labels_test./360*2*pi;     %Labels are converted in rad
        
        % cyclic multivariate Gaussian process regression
        
        % estimate spatial covariance in sl
        S2_est = ME_cmGPR_cov(train_samples, mux_sl(:,:,r), s2x_sl(r,:), 'logistic');
        
        % reconstruct direction from sl activity
        predicted_directions{ii,r} = ME_cmGPR_pred(test_samples, xp, mup_sl(:,:,r), S2_est);
        true_directions{r}=labels_test;
        
        
        % empty variables to free workspace
        itrain=[]; itest=[]; vectors_train=[]; labels_train=[]; S2est=[];
        train_samples=[]; train_labels=[]; test_samples=[];labels_test=[];
        
    end
    
    % quantify crossvalidated accuracy and balanced crossvalidated
    % accuracy
    current_predictions = [predicted_directions{ii,1:end}]';
    current_true_dir =[true_directions{1:end}]';
    acc_vol(ii) = avg_norm_circ_resp_dev(current_predictions(:),current_true_dir(:));
    bal_acc_vol(ii)= bal_norm_circ_resp_dev(current_predictions(:),current_true_dir(:),'trapz');
    
    acc_vol(ii) = acc_vol(ii)*100;
    bal_acc_vol(ii) =  bal_acc_vol(ii)*100;
    
    
    hyp_sl=[]; mup_sl=[]; mux_sl=[]; s2x_sl=[]; vv_sl=[]; n_sl_vox=[];
    current_predictions=[]; current_true_dir=[];
end
cd(target_dir)

chance = 50;
%write volume with Accuracy above chance
resultsvol_m=zeros(sz)+chance;
resultsvol_m(gid) = acc_vol;
str=fullfile(target_dir, [subjectscan, '_sl_', num2str(sl_radius)]);
resultsvol_hdr.fname=[str 'accuracy_minus_chance.img'];
spm_write_vol(resultsvol_hdr,resultsvol_m - chance); % Results are percent above chance
%write volume with balanced Accuracy above chance

resultsvol_m=zeros(sz)+chance;
resultsvol_m(gid) = bal_acc_vol;
resultsvol_hdr.fname=[str 'balanced_accuracy_minus_chance.img'];
spm_write_vol(resultsvol_hdr,resultsvol_m - chance); % Results are percent above chance

disp('Saving prediction results')

save('predictions.mat','predicted_directions','true_directions','-v7.3');
end
