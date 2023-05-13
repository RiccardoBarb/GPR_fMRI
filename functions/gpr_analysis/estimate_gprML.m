function estimate_gprML(scan_dir,subjectscan,target_dir,sl_radius,ffs,c1,c2)

% this function start the GPR estimation on preprocessed beta-images of
% fMRI data.

% The analysis is written to run on a server with 24 workers. We
% parallelize the voxels for loop to speed up the analysis.

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

maskvol_hdr=spm_vol(maskmap);%Get headers of mask
mask_vol=spm_read_vols(maskvol_hdr);%%Get actual mask image

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

fprintf('Subject %s Total voxels in volume: %d\n',subjectscan, nvox);

%% Part 2: Get Beta images estimated with 1st level-GLM and setup analysis
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

%------------------------% Analysis setup %-------------------------------%
%index for run-wise cv
ind = logical(repmat([0 ones(1, nruns-1)], [1, 1]))'; 

% reconstruction grid
prec=1;                         % estimation precision in degrees
xp = [prec:prec:360]'.*(pi/180);% build reconstruction grid
m  = numel(xp);                 % number of candidate angles

%preallocate outcome variables
mup = zeros(m,nvox,nruns);                 % predicted continuous voxel response
mux = zeros(ntrials*(nruns-1),nvox,nruns); % fitted trial-wise voxel response
s2x = zeros(nruns,nvox);                   % voxels variance
hyp(nruns,nvox) = struct('mean',...
    [], 'cov', [0], 'lik', -1);            % voxels Kernel hyperparameters
%if starting from day 2 itrain and itest don't work properly. Here we
%force s1 to match first (=1) run
s1=1;
%% Part 3: Perform voxel-wise cyclic GPR estimation

fprintf(1,'\n %s %s %s %s %s %s %s \n','GPR estimation for',subjectscan ,c1,'coherence','for',c2,'direction');

parfor ii=1:nvox
    
    if mod(ii,round(nvox/100))==0
        fprintf('\n%d\n',ii); % multiple voxels are estimated simultaneusly
    end
    
    vv = squeeze(vectors(:,ii,:));
    
    for r=1:nruns

        %training index for current cv fold
        itrain =  [1:nruns]~=r-s1+1;

        %vectorize data for current cv fold
        vectors_train=vv(:,itrain)';
        Samples=vectors_train(:);
        
        labels_train=extracted_labels(:,itrain)';
        Labels=labels_train(:)./360*2*pi;     %Labels are converted in rad
        
        % cyclic multivariate Gaussian process regression
        
        [mup(:,ii,r), mux(:,ii,r), s2x(r,ii), hyp(r,ii)] = ME_cmGPR_mean(Samples, Labels,xp);
        
        % empty variables to free workspace
        
        itrain=[]; vectors_train=[]; labels_train=[]; Samples=[]; Labels=[];
        
    end
    vv=[];
end

%create target directory and save estimated results
if ~exist(target_dir)
    mkdir(target_dir)
end

disp('Saving estimation results')

save(fullfile(target_dir,'prep_mu_periodic2b.mat'),'mup','mux','s2x','hyp','-v7.3');
disp('done')
