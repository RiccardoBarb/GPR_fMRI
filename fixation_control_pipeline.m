%% Analyze fixation 


% FMT 230320: This script is written to test if a subject did contain
% fixation during a certain trial or not.The script contains parts from 
% open accessible code from:

% Anne Urai:
% - Urai, A.E., Braun, A. & Donner, T.H. Pupil-linked arousal is driven by decision uncertainty and alters serial choice bias. Nature Communications 8,14637 (2017). DOI: 10.1038/ncomms14637
% - https://github.com/anne-urai/pupil_preprocessing_tutorial
%
% It relays on 
%
% Fieldtrip-Toolbox:
% - Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.
% - ftp://ftp.fieldtriptoolbox.org/pub

% Please cite these people if you use the code


% felix.toepfer@bccn-berlin.de

% Berlin, 28.02.2019
%% Eyelink PREPROCESSING
%
% 1. convert edf file to asc
% 2. create a FieldTrip-style data structure
% 3. interpolate Eyelink-defined and additionally detected blinks
%
% Anne Urai, 2016
% modified Felix T�pfer, 2019

%% setup path
clear; clc; close all;

% add the current path, because this is your working directory
thispath = pwd;
path_to_toolbox = %path to Fieldtrip-Toolbox;

addpath(genpath(strcat(thispath,'/functions')));

% set path to FieldTrip - get this from http://www.fieldtriptoolbox.org/download
% it is important to just 'addpath' the fieldtrip folder. Do not 'genpath'
% the fieldtrip folder. Because this adds also the subfolders and leads to
% function referencing problems

%addpath(fullfile(thispath, 'fieldtrip-20190205'))
addpath(genpath(path_to_toolbox))

% run default setting of fieldtrip
ft_defaults;

% set path to the folder where your data is located and a string that
% specifies the experimental folders in that path that you want to analyse

directory = '../data/logfiles_experiment/behavior_and_eyetracking/';
search_index = 's3';

% get the folders that contain the data that you want to work on an check
% if they contain asc files
folders = get_files(directory, search_index);


%% Get eye data and preprocess it

for z=1:length(folders)
    
    if ~folders(z).is_Tasc
        fprintf('Error: %s does not contain an eyetracking file\n',folders(z).subject)
        continue
        
    else
        
        %% preprocess eyetracking data
        
        
        ascFile=folders(z).Tasc_dir;
            
        [dat, event] = prepro_eye(ascFile);
        
        
        %% get the timing from the behavior *.mat files

        
        behavFile = folders(z).mat_dir;
        
        [behavior, display]= extract_behavior(behavFile);

        
        %% get interesting periods
        

        [timing]= get_traces(dat, event, behavior);
        

        %% fixation control
        
        radius = 2; % in degress of visual angle (dva)

        duration = 200; % in ms
        
        [timing] = control_fixation (dat, display, timing, radius, duration);
        

        %% save results

        savename=[ascFile(1:end-3) 'mat'];


        save(savename,'dat','behavior','timing','display')
        
        clear dat
        clear timing
    end
end