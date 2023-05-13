% Preprocessing pipeline

clear all
close all

SPM_path =  %path to SPM
addpath(genpath(strcat(pwd,'/functions')));  %path to functions
addpath(genpath(SPM_path));
n_subjects = 23;
%% data conversion
scan_dir = %path to raw data folder;
target_dir = '../data/fmri_data/preprocessed';

% the subjects you want to convert, based on the amount of folders located
% in the directory

subjects = [1:n_subjects*2];% this converts data for all subjects

% the specific runs you want to convert, in this case 5 is the structural T1
% 10 is the first experimental run and 30 is the first localizer run. 
%The number is based on the amount of subfolder located in the directory.
nruns = [5,10,14,18,22,26,30,34];
[data] = rdk_convertrawdata(scan_dir,subjects,nruns,'dcm',target_dir);

%% realignment
% here the realignment happens differently for day 1 and day 2, in
% particular we reslice the day 1 images but only estimate the day 2

clear all;
scan_dir = '../data/fmri_data/preprocessed/';
n_subjects = 23;
subject = [1:n_subjects];% this performs realignment on all subjects
runs = 1:14; %this includes 10 experimental runs (5 each day) and 4 localizer
           %(2 each day)
total_runs_day1 = 7;
left_out_tr = 4; %we exclude the first 4 volumes
for s_=1:length(subject)
     sub = subject(s_);
     data = rdk_getscans(scan_dir,sub,runs,total_runs_day1...
         ,'f','img','fmri',left_out_tr);
     data = rdk_scan_preprocessing(data,'realign_custom');
end
%% coregistration
% in the coregistration phase we perform different operations for day1 and
% day2: we coregister the realigned and resliced images of day 1 with the
% t1 image of day1. Then we coregister the realigned images of day 2
% (only estimate) with the realigned and resliced images of day 1 and
% reslice them.
clear data data1 data2 structural
for s_ = 1:length(subject)
    sub=subject(s_);
    %here we load the structural and source images for the coregistration
    data1 = rdk_getscans(scan_dir,sub,1:7,total_runs_day1,'r','img','fmri',0);%all the other images day1
    data2 = rdk_getscans(scan_dir,sub,8:14,total_runs_day1,'f','img','fmri',...
        left_out_tr);%all the other images day2
    structural = rdk_getscans(scan_dir,sub,1,total_runs_day1,'s','img','t1',0);% this will be the reference for both day 1 and day 2
    data1.structural = structural.scans;
    source1 = rdk_getscans(scan_dir,sub,1,total_runs_day1,'m','img','fmri',0);% this is the source of day 1
    source2 = rdk_getscans(scan_dir,sub,8,total_runs_day1,'f','img','fmri',0);% source of day 2
    data1.source = source1.scans;
    data2.source{1} = [source2.scans{1,8,400}];%take only the first image of day2
    data2.structural = data1.source;%structural.scans;
    rdk_scan_preprocessing(data1,'coreg_custom');
    rdk_scan_preprocessing(data2,'coreg_custom');
end
