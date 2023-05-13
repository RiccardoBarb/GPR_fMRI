%% Create_meta(data_dir)

% this script is writte to create a meta file, that holds the information
% in which trial a certain subject did exceed fixation. Further it contains
% the trials in which the data was to noisy. Trials where a certain subject
% exceed fixatoin of the data was to noisy will not be analysed further.


% FMT 230320


% do you want to save the results in the data directory? 0=no / 1=yes

saveflag=1;

% set path to the folder where your data is located and a string that
% specifies the experimental folders in that path that you want to analyse

data_dir = '../data/logfiles_experiment/behavior_and_eyetracking/';
search_index = 's3';
target_dir = '../data/logfiles_experiment/behavior_and_eyetracking/';

% get the folders that contain the data that you want to work on an check
% if they contain asc files
folders = get_files(data_dir, search_index);


%% determine the mean noise level
std=[];

for z=1:length(folders)
    
    if ~folders(z).is_Tasc
        fprintf('Error: %s does not contain an eyetracking file\n',folders(z).subject)
        continue
        
    else
        folders(z).subject(1:end-2)
        file=[folders(z).subject(1:end-2) 'd' folders(z).subject(end) 'T1.mat'];
        load(fullfile(data_dir, folders(z).subject, file));
       
        std=[std [timing(:).std]];
        
    end
    
end

save(strcat(target_dir,'/std.mat'), 'std');

%% Set threshold of gaze noise by standard deviation

% because there are some traces where the subjects eye leshes where hiding
% the eyes or the subject closed there eye to much, we obtained some traces
% that are very noisy. To prevent loosing to much data we want to exclude
% those runs from the exclusion analysis.
% the std data of all subjects does not follow any fit ( it looks not
% normal, not log normal, ...). Therefore we fit a smooth density parameter over
% the histogram and get threshold as the border where the area under the
% curve (complete area = 1) is 0.9. that means we get 90%  of the data. and
% define the rest as noisy

% define threshold


load(strcat(target_dir,'/std.mat'));



[noise_thresh, xi]=ksdensity(std,0.9,'Function','icdf');




for z=1:length(folders)
    
    if ~folders(z).is_Tasc
        fprintf('Error: %s does not contain an eyetracking file\n',folders(z).subject)
        continue
        
    else
        folders(z).subject
        file=[folders(z).subject(1:end-2) 'd' folders(z).subject(end) 'T1.mat'];
        load(fullfile(data_dir, folders(z).subject, file));

        meta(z).subject=folders(z).subject;       
        meta(z).noise_thresh=noise_thresh;

        
        for r_=1:length(timing)
            
            meta(z).exclud_trial_id_per_run(r_).noise_level=timing(r_).std;
            
            if timing(r_).std<noise_thresh
                meta(z).exclud_trial_id_per_run(r_).trial=timing(r_).stim_exceeded_trial_id;
            else
                meta(z).exclud_trial_id_per_run(r_).trial=[];
            end
        end
        
        if isempty(cell2mat({meta(z).exclud_trial_id_per_run(:).trial}))
            meta(z).is_exclude=0;
        else
            meta(z).is_exclude=length(cell2mat({meta(z).exclud_trial_id_per_run(:).trial}));
        end
        
    end
    
end

% exclude training run ( run 1 )
 meta_fixation = transform_meta(meta);

if saveflag==1

savename = [strcat(target_dir,'meta_fixation.mat')];
save(savename,'meta_fixation');

end
