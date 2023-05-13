%% Full pipeline to estimate gpr based respone profiles of x and y 
%% coordinates of gaze data of one (dominant) eye seperate for each time point of interest.
%% Further this pipeline incorporats the reconstruction of physical 
%% stimulus direction and choice by means of maximum likelyhood estimation

clear; clc; close all;

% shuffle randomization seed

analysis.rng = rng('shuffle');

% add paths
path_to_toolbox = %path to Fieldtrip-Toolbox;
GPR_ML_path = %path to GPR_ML_library http://gaussianprocess.org/gpml/code/matlab/doc/
% add the current path, because this is your working directory
thispath = pwd;
addpath(genpath(thispath));
addpath(genpath('../functions'))
addpath(genpath(GPR_ML_path))
% add fieldtrip toolbox (https://www.fieldtriptoolbox.org/) this is necessary
% for the processing of the gaze data
addpath(path_to_toolbox)
% run default setting of fieldtrip
ft_defaults;

% find find the folders where the files of interest are located

data_dir = '../../data/eye_tracking_for_GPR';
search_index = 's3';
folder_dir = get_data_folder_directory(data_dir,search_index);

% find the files of interest

search_str = '_eye';
file_extension = '.mat';
[files_dir, subj_name]  = get_files_directory(folder_dir,search_str,file_extension);

%% calculation

% save the estimated results (1)
% do not save the estimated results (0)

save_flag=0;

% set the direction labels to be used for the estimation and reconstruction 
% 'stimulus'= physical stimulus labels, 'report'=reported direction of the
% subjects (=choices)

labels={'stimulus','report'};

% start parallel processing

start_parpool(24)   


% run GPR based analysis with the use of the distinct labels

for lab_=1:length(labels)

    bal_acc = eye_gprML_efficient(files_dir,labels{lab_}, save_flag);
  
end

% end parpool

end_parpool

%save successfull analysis

analysis.labels = labels;
analysis.date = date;
analysis.time = fix(clock);
analysis.files_dir = files_dir;

save([pwd,'/analysisinfo',sprintf('%d',analysis.time(1:3)) '.mat'], 'analysis')


