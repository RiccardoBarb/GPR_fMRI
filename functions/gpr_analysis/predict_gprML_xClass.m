function predict_gprML_xClass(scan_dir,subjectscan,sl_radius,ffs,train_on,test_on,roi_file,source_SL)
% this function start the GPR prediction on preprocessed beta-images of
% fMRI data. The prediction is performed by applying multivariate MLE on
% 1) the continuous activity pattern (mup) and
% 2) the estimated spatial covariance (S2_est) in a given searchlight.

% The voxel-wise estimations of mup mux s2 hyp, need to be performed with the
% function estimate_gprML.
%
% This particular function is used to constrain the reconstructoin analysis on
% searchlights obtained from a roi. We used this function to perform the
% cross-prediction analyses described in the paper materials and methods
% "additional exploratory analyses: fMRI data" subsection "model
% generalization".

%% Part 1: SL loading/creation and roi referencing into the model voxel space

maskmat = fullfile(source_SL, [subjectscan '_volume_search_' num2str(sl_radius) '.mat']);
maskmap = fullfile(scan_dir, sprintf('mask.nii'));
maskvol_hdr=spm_vol(maskmap);%Get headers of mask
mask_vol=spm_read_vols(maskvol_hdr);%%Get actual mask image
clear maskvol_hdr
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

nvox = numel(gid);
%get the roi mask and combine the two masks
roi = spm_read_vols(spm_vol(roi_file{:}));
ids = mask_vol & roi;
nvox_roi = numel(find(ids));
%we reference the combination of the roi and the mask onto the
%linear index. Now the voxels of the roi mask can be referenced to the
%model (mup) voxel space.
roi_voxels = ids(gid);
v_ids = find(roi_voxels);
fprintf('Subject %s Total voxels in ROI: %d\n',subjectscan, nvox_roi);

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

%% get training data and labels

part_train = strsplit(train_on,'_');

if strcmpi(part_train{1},'low'), coherence = 0;
elseif strcmpi(part_train{1},'mid'), coherence = 16;
elseif strcmpi(part_train{1},'high'), coherence = 32;end
if strcmpi(part_train{2},'stimulus'), condition = 3;
elseif strcmpi(part_train{2},'report'), condition = 4;end

betanumbers = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]+coherence;
ntrials = numel(betanumbers);
%--------------------% Get Betas and lables %-----------------------------%
fprintf('Train data for subject %s          \n\n',subjectscan);

[extracted_labels_train,betas_train] = ...
    extract_glm_design(scan_dir,condition,ntrials,betanumbers,gid,nruns);

fprintf('\nDone!');

%---------------------% Get mup mux s2 hyp %------------------------------%
fprintf('\nloading estimated GPR hyperparameters and predictions...')

source_files = fullfile(['../data/fmri_data/out/gpr/',subjectscan,'/',train_on]);

load (fullfile(source_files,'prep_mu_periodic2b.mat'));

disp('\nDone!')

%% get testing data and labels
for t_ = 1:length(test_on)
    current_test_on = test_on{t_};
    fprintf('Setting up GPR for %s, Xclass Train on %s Test on %s \n',...
        subjectscan,train_on, current_test_on)
    target_dir = fullfile(['../data/fmri_data/out/gpr/Xclass'...
        ,subjectscan,'/','train_',train_on,'_test_',current_test_on]);
    
    part_test = strsplit(current_test_on,'_');
    
    if strcmpi(part_test{1},'low'), coherence = 0;
    elseif strcmpi(part_test{1},'mid'), coherence = 16;
    elseif strcmpi(part_test{1},'high'), coherence = 32;end
    if strcmpi(part_test{2},'stimulus'), condition = 3;
    elseif strcmpi(part_test{2},'report'), condition = 4;end
    
    betanumbers = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]+coherence;
    ntrials = numel(betanumbers);
    
    %--------------------% Get Betas and lables %-----------------------------%
    fprintf('Test data for subject %s         \n\n',subjectscan);
    
    [extracted_labels_test,betas_test] = ...
        extract_glm_design(scan_dir,condition,ntrials,betanumbers,gid,nruns);
    
    fprintf('\nDone!');
    
    %------------------------% Analysis setup %-------------------------------%    
    % reconstruction grid
    prec = 1;                         % estimation precision in degrees
    xp = [prec:prec:360]'.*(pi/180);  % build reconstruction grid
    
    %preallocate outcome variables
    
    bal_acc_vol = nan(nvox, 1);             % SL balanced cross-validated accuracy
    predicted_directions = cell(nvox,nruns);  % SL trial-wise predicted direction
    true_directions = cell(nruns,1);          % tesing lables of each cv-fold
    
    %if starting from day 2 itrain and itest don't work properly. Here we
    %adjust s1 and s1 to match first (=1) run
    s1=1;
    
    %% Part 3: Perform searchlight-based cyclic GPR prediction
    
    fprintf(1,'\n %s %s %s %s %s %s %s \n','GPR Prediction for',subjectscan ,...
        'xClass train on',train_on,'test on',current_test_on);
    
    for ii=1:nvox_roi
        current_index = v_ids(ii) ;
        fprintf('\n %0.1f%%\n',(ii/nvox_roi).*100); % \b is backspace
        
        % Bring only data from the current searchlight into the parallel
        % for loop
        vv_sl_tr = betas_train(:,gis{current_index},:);  % trial-wise betas training set
        vv_sl_te = betas_test(:,gis{current_index},:); % trial-wise betas testing set
        
        hyp_sl = hyp(:,gis{current_index});       % voxel-wise kernel hyperparameters
        mup_sl = mup(:,gis{current_index},:);     % predicted continuous voxel response
        mux_sl = mux(:,gis{current_index},:);     % fitted trial-wise voxel response
        s2x_sl = s2x(:,gis{current_index});       % voxels variance
        
        n_sl_vox = numel(gis{current_index});     % number of voxels in the sl

        for r=1:nruns
            
            %training and testing index for current cv fold
            itrain = find( [1:nruns]~=r-s1+1);
            itest  = r-s1+1;
            
            %vectorize training data for current cv fold
            vectors_train = vv_sl_tr(:,:,itrain);
            train_samples = permute(vectors_train,[3,1,2]);
            train_samples = reshape(train_samples,[ntrials*(nruns-1),n_sl_vox]);
            
            labels_train=extracted_labels_train(:,itrain)';
            labels_train=labels_train(:)./360*2*pi;     %Labels are converted in rad
            
            %vectorize testing data for current cv fold
            
            test_samples = vv_sl_te(:,:,itest);
            labels_test=extracted_labels_test(:,itest);
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
        
        % quantify balanced crossvalidated accuracy
        current_predictions = [predicted_directions{ii,1:end}]';
        current_true_dir = [true_directions{1:end}]';
        bal_acc_vol(current_index) = bal_norm_circ_resp_dev(current_predictions(:),current_true_dir(:),'trapz');
        bal_acc_vol(current_index) =  bal_acc_vol(current_index)*100;
        
        
        hyp_sl=[]; mup_sl=[]; mux_sl=[]; s2x_sl=[]; vv_sl=[]; n_sl_vox=[];
        current_predictions=[]; current_true_dir=[];
    end
    
    mkdir(target_dir)
    cd(target_dir)
    
    
    chance = 50;
    bal_acc = bal_acc_vol(~isnan(bal_acc_vol))-chance;
        
    disp('Saving prediction results')
    save('bal_acc.mat','bal_acc');
    save('predictions.mat','predicted_directions','true_directions','-v7.3');
end
