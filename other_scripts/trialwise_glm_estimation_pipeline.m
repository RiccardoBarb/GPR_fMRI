%% design matrix editor v3.0
% Riccardo: updated on August 2019

clear all
close all

SPM_path =  %path to SPM
addpath(genpath(strcat(pwd,'/functions')));  %path to functions
addpath(genpath(SPM_path));


%% General info about dataset

%logsfile
logs_dir = '../data/logfiles_experiment/behavior_and_eyetracking/';%location of the logfiles
%functional scans
scan_dir = '';%location of functional brain scans
all_s = dir(strcat(scan_dir,'s*'));
subjectscan = 1:length(all_s);% vector of subject order in the dataset
%folder eg. [1,2,3...]
subj_id = {all_s(subjectscan).name};% cell with subjects id eg.{'s303'},{'s318'}
%number of runs
runs_experimental = [1:5,8:12];%run number for the trialwise design
total_runs_day1 = 7;%total number of runs for day 1 (including the localizer runs)
left_out_tr = 4;

%% specify single-trial GLM
for s_ = 1:length(subjectscan)
    target_singletrial = strcat('../data/fmri_data/trialwise_glm/',...
        sprintf('%s',subj_id{s_}));
    [time, task] = extract_design_info(logs_dir, subj_id{s_});
    movfull = rdk_getscans(scan_dir,subjectscan(s_),runs_experimental,...
        total_runs_day1,'r','txt','fmri',0);%head-motion parameters
    datafull = rdk_getscans(scan_dir,subjectscan(s_),runs_experimental,...
        total_runs_day1,'rf','img','fmri',0);
    %setup design matrix and estimate GLM
    [matlabbatch] = rdk_create_design_singletrial(time,task,datafull.scans,...
        movfull.scans,subjectscan(s_),runs_experimental,...
        target_singletrial,left_out_tr);
end