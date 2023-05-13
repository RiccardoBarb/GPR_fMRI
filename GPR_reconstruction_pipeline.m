% This is the main script for performing the GPR reconstruction pipeline. 
% The method is described in details on the manuscript. 
% We have also provided a tutorial script which should help in 
% understanding the basic functioning of the reconstruction procedure 
% (See part 3 of the tutorial_gprML.m).
% The analysis will automatically run for each specified subject and
% condition, and will first compute a searchlight of radius 4 centered around 
% each voxel. Then, it will go through each subject Searchlight index,
% extract the corresponding response profiles (estimated with the
% script GPR_estimation_pipeline) for the training set and perform stimulus
% or report reconstruction via MLE on the test set. This will result in a
% set of predicted directions for each searchlight index which will be 
% compared with the actual condition directional label 
% (either stimulus or report direction).
% This will result in a measure of performance called BFCA 
% (see manuscript eq. 15) for each searchlight index, which we will write 
% into an image that can be interpreted as an "accuracy map". For the sake 
% of completeness our script will produce 2 separate outputs for
% each analysis: accuracy_minus_chance and balanced_accuracy_minus_chance.
% The first is computed with a standard measure of continuous accuracy (FCA
% described in eq. 2 of the manuscript), the second with the balananced
% measure of accuracy (BFCA-eq 15) the the comparison of the two outputs is
% described in Appendix 2 of the manuscript.
% N.B. the analysis takes a long time to run for all of the subject and conditions:
% for each subject a total of 6 analyses are performed (2 conditions -
% stimulus and report - X 3 coherence level) and the average estimation
% time for a single analysis is about 10h when running on a cluster
% computer with 24 workers and 250GB of RAM. The script is developed to run
% on such cluster computer and parallelizes the computations performed for each
% voxel. You can however deactivate the parallel computing pool by
% commenting the first 20 lines of the function "predict_gprML".
% The analyse outputs 2 separate img files "accuracy_minus_chance" and 
% "balanced_accuracy_minus_chance" as well as a relatively large mat file 
% whichs stores the predicted voxel-wise predicred directions for each CV 
% fold. Such files will be saved in the path specified by the variable
% "recon_dir" under the appropriate subfolder associated with
% the subject id and the condition id (combination of coherence_reconstructed
% label). ES. The output of the analysis on high coherence - stimulus
% direction for subject 301 will be located in
% recon_dir/s301/high_stimulus/.
% N.B. To run the reconstruction procedure is necessary to compute first the
% GPR estimation pipeline.
% add paths
SPM_path =  %path to SPM
GPR_ML_path = %path to GPR_ML_library http://gaussianprocess.org/gpml/code/matlab/doc/
addpath(genpath(strcat(pwd,'/functions'))); %path to functions
addpath(genpath(SPM_path));
addpath(genpath(GPR_ML_path));
scan_dir = '../data/fmri_data/trialwise_glm/';
recon_dir = '../data/fmri_data/gpr/';

%% setup reconstruction design

ffscan = [1,10];%first and final scans
sl_radius = 4; %searchlight radius (in voxels) the reconstruction is multivariate
subject_pool = dir(strcat(scan_dir,'s*')); 
subjects = {}; 
for s=  1:length(subject_pool) 
    sub = strsplit(subject_pool(s).name,'_'); 
    subjects{s} = sub{1}; 
end
coherences = {'high','mid','low'};
conditions = {'stimulus','report'};

%% start GPR - based reconstruction

for s_=1:length(subjects)
    beta_dir = strcat(scan_dir,sprintf('%s',subjects{s_}));
    for cond_=1:length(conditions)
        for ch_=1:length(coherence)
            target_dir = fullfile([recon_dir...
                ,subjects{s_},'/',coherences{ch_},'_',conditions{cond_}]);
            
            predict_gprML(beta_dir,subjects{s_},target_dir,...
                sl_radius,ffscan,coherences{ch_},conditions{cond_});
            
        end
    end
end