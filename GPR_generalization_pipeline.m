%% predict with xClass design - model generalization and consistency.
% in this script we use the GPR reconstruction pipeline to predict
% the reported and the stimulus directions based on models estimated in a
% different coherence-condition; e.g. train on the GPR estimated from 
% condition low-report and reconstruct condition high-stimulus.
% The procedure requires users to have already estimated the GPR tuning
% profiles (see GPR_estimation_pipeline.m) for all of the conditions and 
% subjects of interest.
% The analysis is constrained to a specific ROI obtained by specifying
% two one-way factorial designs (one for stimulus and one for report) 
% with coherence level as a within-subject factor. The ROI is the 
% intersection of the voxels resulting from the contrast [1,1,1] defined on
% both models. The ROI was subsequently projected in the individual subject
% voxel space (see other_scripts/individual_roi_inverse_normalization.m and
% other_scripts/individual_roi_reslice.m). The above-mentioned rois are 
% included in the dataset and can be found in 
% /data/fmri_data/xclass_individual_masks.
% The first part of the script is used to perform the cross-reconstruction
% analysis by specifying a 'train_on' string variable and a 'test_on' cell
% containing one or multiple strings. This procedure should be performed by
% specifying both training directions for each condition: e.g. to evaluate 
% the model generalization performance between low-report and high-stimulus
% one should run two separate analyses, one with train_on = 'low_report', 
% test_on = {'high_stimulus'} and the other with train_on = 'high_stimulus', 
% test_on = {'low_report'}.
% After this, the second part of the script is used to average the
% reconstruction performances resulting from the two analyses.
% The third and final part of the script is used to produce a summary
% plot displaying the model generalization and model consistency analysis
% described in Figure 8 of the manuscript.
%% Part 1: cross-reconstruction analysis
clear all
SPM_path =  %path to SPM
GPR_ML_path = %path to GPR_ML_library http://gaussianprocess.org/gpml/code/matlab/doc/
addpath(genpath(strcat(pwd,'/functions'))); %path to functions
addpath(genpath(SPM_path));
addpath(genpath(GPR_ML_path));
source_SL = %path to the folder containing all searchlight volumes computed with a radius of 4 voxels
scan_dir = '../data/fmri_data/trialwise_glm/';
roi_dir = '../data/fmri_data/xclass_individual_masks/';
subject_pool = dir(strcat(scan_dir,'s*')); 
subjects = {}; 
for s=  1:length(subject_pool) 
    sub = strsplit(subject_pool(s).name,'_'); 
    subjects{s} = sub{1}; 
end
% train_on should be a string with coherencelabel_conditionlabel
% possible options: 'low_report','mid_report','high_report','high_stimulus',
% 'low_stimulus','mid_stimulus'
train_on = ''
% test_on should be a cell with one or multiple strings describing
% coherencelabel_conditionlabel: example {'high_stimulus','mid_stimulus',
% 'low_stimulus','high_report','mid_report','low_report'}
test_on = {}
for s_=1:length(subjects)
    beta_dir = strcat(scan_dir,sprintf('%s',subjects{s_}));
    ffscan = [1,10];%first and final scans
    sl_radius = 4; %searchlight radius (in voxels)
    roi_file = strcat(roi_dir,'/','r',subjects{s_},...
        '_average_eff_stim_rep_inter.img');% the ROI MASK
    predict_gprML_xClass(beta_dir,subjects{s_},sl_radius,ffscan,train_on,test_on,roi_file,source_SL)
end

%% Part 2 - average the accuracy values from both training directions
% This part of the script can only be run if part 1 was completed on all 
% of the possible pairs of coherence-condition. Once completed, the outputs
% of the analyses should be located in /data/fmri_data/gpr/Xclass/s*/
% train_*_test_*/bal_acc.mat - where s* is the subject id train_* is the
% pair of coherence-condition used for training and _test_* is the 
% pair of coherence-conditions used for the reconstruction. N.B. if train_*
% and _test_* have the same coherence-condition (e.g.
% train_high_report_test_high_report it's not necessary to average
% reconstruction performance values from both training directions, since
% it doesn't represent a model generalization/consistency design.

clear all
close all
all_s = %subject list;
base_folder=%output folder of the model generalization analyses;
conditions = {'low_report','high_stimulus','low_stimulus','mid_stimulus',...
    'mid_report','high_report'};
for c_= 1:length(conditions)
    for s_ = 1:length(all_s)
        subject = all_s(s_).name;
        for i_=1:length(conditions)
            if c_==i_
                filez = dir(strcat(base_folder,conditions{c_},'/',subject,...
                    '/train_',conditions{c_},'_test_',conditions{i_}...
                    ,'/bal_acc.mat'));
                bal= load(strcat(filez.folder,'/',filez.name));
                mean_bal_acc= bal.bal_acc;
                target_dir = strcat(base_folder,'means/',subject,'/',...
                    conditions{c_},'_',conditions{i_});
                if ~exist(target_dir)
                    mkdir(target_dir);
                end
                filename = strcat(target_dir,'/','mean_bal_acc.mat');
                save(filename,'mean_bal_acc')
            else
                try
                    filea = dir(strcat(base_folder,conditions{c_},'/',subject,...
                        '/train_',conditions{c_},'_test_',conditions{i_}...
                        ,'/bal_acc.mat'));
                    fileb = dir(strcat(base_folder,conditions{c_},'/',subject,...
                        '/train_',conditions{i_},'_test_',conditions{c_}...
                        ,'/bal_acc.mat'));
                    bal1= load(strcat(filea.folder,'/',filea.name));
                    bal2= load(strcat(fileb.folder,'/',fileb.name));
                    mean_bal_acc = mean([bal1.bal_acc,bal2.bal_acc],2);
                    target_dir = strcat(base_folder,'means/',subject,'/',...
                        conditions{c_},'_',conditions{i_});
                    if ~exist(target_dir)
                        mkdir(target_dir);
                    end
                    filename = strcat(target_dir,'/','mean_bal_acc.mat');
                    save(filename,'mean_bal_acc')
                catch ME
                    if (strcmp(ME.identifier,'MATLAB:load:pathIsDirectory'))
                        fprintf('value_already computed, skipping pair \n')
                    end
                end
            end
        end
    end
end
%% Part 3 - plot  summary 
%this procedure generates the plot used for figure 8. Uncomment the fisrt
%part to generate a matrix with the average balanced accuracy values for
%each subject and condition pair. N.B. this only works if you successfully 
%completed part 1 and pert 2 of the script.
%Otherwise load the provided results of the analyses and plot the figures.
clear all
close all
% all_s = %subject list;
% base_folder= %folder with the mean results of the previous analysis;
% 
% conditions = {'low_report','mid_report','high_report','low_stimulus',...
%     'mid_stimulus','high_stimulus'};
% % this matrix represents the order of the conditions to be averaged
% idsx = [3,5,1,4,6,2;
%         2,5,4,3,6,1;
%         3,5,1,4,6,2;
%         2,5,3,4,6,1;
%         2,5,4,3,6,1;
%         6,4,1,3,5,2];
% 
% 
% for s_ = 1:length(all_s)
%     subject = all_s(s_).name;
%     for i_=1:length(conditions)
%         ff = dir(strcat(base_folder,'/',subject,'/*',conditions{i_},'*'));
%         for ii_ = 1:length(ff)
%             load(strcat(ff(idsx(i_,ii_)).folder,'/',ff(idsx(i_,ii_)).name,'/','mean_bal_acc.mat'));
%             mean_bal_acc_c(s_,i_,ii_) = mean(mean_bal_acc);
%             clear mean_bal_acc
%         end
%     end
% end
%% results of the analyses already calculated
conditions = {'low_report','mid_report','high_report','low_stimulus',...
    'mid_stimulus','high_stimulus'};
load('../data/fmri_data/xclass_results_figure_8.mat')
n_subj = 23;
figure();
for p_ = 1:size(mean_bal_acc_c,2)
    current = squeeze(mean_bal_acc_c(:,p_,:));
    stm = std(current)./sqrt(n_subj);
    subplot(1,size(mean_bal_acc_c,2),p_)
    bar(mean(current),'Facecolor',[0.5,0.5,0.5])
    hold on
    errorbar(mean(current),stm,'color',[0,0,0],'Linestyle','none')
    ylim([-0.5,5])
    set(gca,'xticklabel',{'LR','MR','HR',...
    'LS','MS','HS'});
    title(conditions{p_})
end