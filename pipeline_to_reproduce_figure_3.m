% This script can be used to reproduce the figure 3 of the manuscript
% Barbieri, Toepfer et al., 2023. The plot displays GPR-based tuning 
% response profiles of 16 voxels of subject 318. The following script is
% divided in 3 parts. The first shows the procedure for selecting the 
% searchlight to from which to extract the voxels response profiles.
% The second part extracts the beta estimates and the tuning response
% profiles estimated via GPR from the chosen subject and condition
% (s318-high-stimulus).
%The third part is for plotting purposes. N.B. to reproduce the end-to-end
%pipeline is necessary to perform the whole GPR estimation and
%reconstruction pipeline. Alternatively you can reproduce the figure by
%using the data included in the dataset
% '../data/fmri_data/s318_tuning_curves_figure_3.mat'
%% First part
% the first part of the script identifies the best searchlight as the one
% with the highest mean balanced accuracy value across conditions and 
% coherence. The searchlight is extracted from a ROI defined for 
% exploratory purposes resulting from the second level analysis on the GPR
% accuracy maps (whole-brain searchlight analysis). More specifically, it 
% represents an area coding for both the stimulus and the report (more
% details provided in the method section of the manuscript).
% The ROI was subsequently projected in the individual subject
% voxel space.

clear all
SPM_path = ''
addpath(genpath(strcat(pwd,'/functions')))
addpath(genpath(SPM_path))
coherence = {'high','mid','low'};
condition = {'stimulus','report'};
subject = '';%name of subject id
subject_folder = strcat('../data/fmri_data/trialwise_glm','/',subject,'/');
source_mask_file = strcat(subject_folder,'mask.nii');%this is not the ROI MASK
mask = spm_read_vols(spm_vol(source_mask_file));
roi_file = strcat('../data/fmri_data/xclass_individual_masks','/','r',subject,...
    '_average_eff_stim_rep_inter.nii');% this IS the ROI MASK
for cond_ = 1:numel(condition)
        for coh_= 1:numel(coherence)
            accuracy_file = strcat('../data/fmri_data/gpr/',subject,...
                '/',coherence{coh_},'_',condition{cond_},...
                '/',subject,'_sl_4balanced_accuracy_minus_chance.img')
            accuracy = spm_read_vols(spm_vol(accuracy_file));
            %we get the roi mask and combine the two
            roi = spm_read_vols(spm_vol(roi_file));
            ids = mask & roi;
            accuracy_in_ROI(cond_,coh_,:) = accuracy(ids); 
        end
end
mean_accuracy_in_ROi = mean(squeeze(mean(accuracy_in_ROI,1)),1);
[~,best_SL_id] = max(mean_accuracy_in_ROi);

%% Second part
% The second part of the script extracts beta estimates from the trial-wise
% first level GLM from the chosen condition and the corresponding tuning 
% response profiles estimated via gpr
do_predictions=0;
do_betaestimate=1;
mapping_file = '' %this is the path of the searchlight radius 4 computed 
% during the gpr reconstruction phase.
coherence = {'high'};
condition = {'stimulus'};
if strcmpi(coherence(1),'low'), coh_ref = 0;
elseif strcmpi(coherence(1),'mid'), coh_ref = 16;
elseif strcmpi(coherence(1),'high'), coh_ref = 32;end
if strcmpi(condition(1),'stimulus'), cond_ref = 3;
elseif strcmpi(condition(1),'report'), cond_ref = 4;end

betanumbers = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]+coh_ref;
condition_folder = strcat('../data/fmri_data/gpr/',subject,...
    '/',coherence{coh_},'_',condition{cond_},'/');
%the tuning_file is the result of running the GPR estimation pipeline.
%There should be one per subject for each condition and coherence level
%(total of 6 per subject)
tuning_file = strcat(condition_folder,'/','prep_mu_periodic2b.mat');
[roi_tuning,~,~,betaestimate,extracted_labels] = get_tuning_from_roi...
    (tuning_file,roi_file,source_mask_file,mapping_file,...
    do_predictions,do_betaestimate,betanumbers,cond_ref,best_SL_id);
%% Third part
% the third part is for plotting purposes. In order to reproduce the figure
% we already computed the necessary steps and provided them with a mat file
% containing the results of part 1 and 2 of this script on subject318.
% to reproduce the Figure 3 of the manuscript run the following lines.
load('../data/fmri_data/s318_tuning_curves_figure_3.mat')
rng (69)
nvoxels = 16;
s2 = zeros(nvoxels);
voxel_toplot = zeros(nvoxels,360);
ssqr_toplot = zeros(nvoxels);
vox= randi(size(roi_tuning,2),1,nvoxels) ; %assign 16 random voxels 

tuning_condition = roi_tuning-mean(roi_tuning,1);%mean center
betas = permute(betaestimate,[3,1,2]);
betas_condition = reshape(betas,[160,size(tuning_condition,2)]);
betas_condition = betas_condition-mean(betas_condition,1);
betas_toplot(:,:) = betas_condition(:,vox);
labels=extracted_labels(:)';
labels_toplot(:)=labels;
voxel_toplot(:,:) = tuning_condition(:,vox)';
lab_id = round(labels);
lab_id(lab_id==0)=360;

s2 = std(betas_condition(:,vox));
sqres = (betas_condition(~isnan(lab_id),vox)-tuning_condition(lab_id(~isnan(lab_id)),vox)).^2;
totalres =  (betas_condition(~isnan(lab_id),vox) - mean(tuning_condition(lab_id(~isnan(lab_id)),vox))).^2;
ssqr_toplot = 1-(sum(sqres)./sum(totalres));

clear roi_tuning betaestimate betas_condition tuning_condition betas extracted_labels labels

f= figure('rend','painters','pos',[1280 -160 1600 900],'visible','on');
x = 1:360;%which directions
for v_=1:nvoxels

    vxs = squeeze(voxel_toplot(v_,:));
    subplot(4,4,v_)
    plot(x,vxs,'k','Linewidth',2);
    hold on
    plot(labels_toplot(:),squeeze(betas_toplot(:,v_)),'k.');
    hold on
    patch([x fliplr(x)], [vxs-s2(v_)  fliplr(vxs+s2(v_))], [0,0,0],'facealpha',0.1,'edgealpha',0);
    xlim([0,360])
    ylim([-10,10])
    axis square
end


