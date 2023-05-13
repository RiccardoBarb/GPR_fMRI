% This is the main script for performing the GPR estimation pipeline. For
% each coherence level we perform an estimation of the voxelsresponse 
% profiles across stimulus or report direction. The method is described in
% details on the manuscript. We have also provided a tutorial script which
% should help in understanding the basic functioning of the estimation
% procedure (See part 2 of the tutorial_gprML.m).
% The analysis will automatically run for each specified subject and
% condition, and will go through each subject voxels and estimate the
% response profiles with a 10 fold cross validation procedure. N.B. the
% analysis takes a long time to run for all of the subject and conditions:
% for each subject a total of 6 analyses are performed (2 conditions -
% stimulus and report - X 3 coherence level) and the average estimation
% time for a single analysis is about 10h when running on a cluster
% computer with 24 workers and 250GB of RAM. The script is developed to run
% on such cluster computer and parallelizes the computations performed for each
% voxel. You can however deactivate the parallel computing pool by
% commenting the first 20 lines of the function "estimate_gprML". The
% outputs of the analyses are very large .mat files (between 4 and 6 GB) named
% 'prep_mu_periodic2b.mat' which will be saved in the path specified by the 
% variable "estimation_dir" under the appropriate subfolder associated with
% the subject id and the condition id (combination of coherence_reconstructed
% label). ES. The output of the analysis on high coherence - stimulus
% direction for subject 301 will be located in
% estimation_dir/s301/high_stimulus/prep_mu_periodic2b.mat.
% N.B. running the estimation procedure is necessary to compute the
% subsequent stimulus-report reconstruction pipeline.

% add paths
SPM_path =  %path to SPM
GPR_ML_path = %path to GPR_ML_library http://gaussianprocess.org/gpml/code/matlab/doc/
addpath(genpath(strcat(pwd,'/functions')));  %path to functions
addpath(genpath(SPM_path));
addpath(genpath(GPR_ML_path));
scan_dir = '../data/fmri_data/trialwise_glm/';
estimation_dir = '../data/fmri_data/gpr/';

%% setup estimation design

ffscan = [1,10];%first and final scans
sl_radius = 1; %searchlight radius (in voxels) the estimation is done
%on single voxels
subject_pool = dir(strcat(scan_dir,'s*')); 
subjects = {}; 
for s = 1:length(subject_pool) 
    sub = strsplit(subject_pool(s).name,'_'); 
    subjects{s} = sub{1}; 
end
coherences = {'high','mid','low'};
conditions = {'stimulus','report'};

%% start GPRML estimation for every voxel

for s_  = 1:length(subjects)
    beta_dir = strcat(scan_dir,sprintf('%s',subjects{s_}));
    for cond_ = 1:length(conditions)
        for ch_ = 1:length(coherence)
            target_dir = fullfile([estimation_dir...
                ,subjects{s_},'/',coherences{ch_},'_',conditions{cond_}]);
            
            estimate_gprML(beta_dir,subjects{s_},target_dir,...
                sl_radius,ffscan,coherences{ch_},conditions{cond_});
            
        end
    end
end
