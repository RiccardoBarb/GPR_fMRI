% The following script can be used to perform the behavioral analyses 
% described in the manuscript. More precisely, the script will go through
% the experimental data and extract the trial-by-trial stimulus motion
% direction and the corresponding reported direction of motion for each 
% coherence level and each subject.
% Once the data are extracted, the script produces the plots displayed in
% figure 2 of the manuscript. The first plots are response deviations from
% the true stimulus for each coherence level and the boxplot of the FCA
% pooled across all subjects (see eq.2 in the methods section). Then we
% produce descriptive scatterplots of the pooled response distributions
% against the true stimulus directions for all of the coherence levels (NB.
% in the 0 coherence condition there is no true stimulus motion direction).
% We also provide code to perform the same scatterplots for single subject
% responses.
% The final part of the script includes the code to run checks on the amount
% of excluded trials due to eye movements.
clear all
%% Behavioral data extraction
addpath(genpath(strcat(pwd,'/functions')));
coherence = [0,16,32];
condition = [3,4];
slist = dir('../data/fmri_data/trialwise_glm/s*');
n_subj = length(slist);
acc = nan(n_subj,length(coherence));
all_stim = cell(n_subj,length(coherence));
all_rep = cell(n_subj,length(coherence));
for i_=1:length(slist)
    cs =  slist(i_).name;
    name = strsplit(cs,'_');
    subject =   name{1};
    scan_dir = sprintf('../data/fmri_data/trialwise_glm/%s',subject);
    %--------------------------% Get lables %-----------------------------%
    fprintf('Extracting subject %s        \n\n',subject);
    for c_ = 1:length(coherence)
        betanumbers = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]+coherence(c_);
        ntrials = numel(betanumbers);
        [stimulus] = ...
            extract_behavior_design(scan_dir,condition(1),ntrials,betanumbers,10);
        [report] = ...
            extract_behavior_design(scan_dir,condition(2),ntrials,betanumbers,10);
        stimulus = stimulus';
        report = report';
        s_rad = stimulus(:)./360*2*pi;
        r_rad = report(:)./360*2*pi;
        acc(i_,c_) = avg_norm_circ_resp_dev(s_rad,r_rad);
        all_stim{i_,c_} = stimulus;
        all_rep{i_,c_} = report;
    end
    fprintf('\n');
end
fprintf('\nDone!');
%% the following code produces Figure 2 B) and C) from the Paper
acc = [acc(:,1),acc(:,2),acc(:,3)];
sub_var =std(acc,[],1);
sub_std = sub_var./sqrt(23);
meanss = mean(acc);
figure('Renderer', 'painters', 'Position', [10 10 1024 760]),
boxplot(acc)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),[0.5,0.5,0.5],'FaceAlpha',.5);
end
%plot means as black asterisks.
hold on
plot(1,meanss(1), 'k*')
plot(2,meanss(2), 'k*')
plot(3,meanss(3), 'k*')

ylim([0.4,1])

figure('Renderer', 'painters', 'Position', [10 10 1024 760]),
t = {'zero','mid','high'}
for c_ = [3,2,1]
    subplot(1,3,c_)
    r = vertcat(all_rep{1:end,c_});
    s = vertcat(all_stim{1:end,c_});
    dev = angdiff2(s(:),r(:));
    h = histogram(dev,[-180:10:180],'Normalization','Probability','Facecolor',[0.5,0.5,0.5]) 
    ylim([0,0.4])
    xlim([-180,180])
    title(t{c_})
    
end
%% the following code produces Figure 2 A) from the Paper
t = {'zero','mid','high'};
save_path = pwd;
for c_ = [3,2,1]
    figure('Renderer', 'painters', 'Position', [10 10 1024 760]),
    r = vertcat(all_rep{:,c_});
    s = vertcat(all_stim{:,c_});
    dev = angdiff2(s(:),r(:));
    
    s=scatterhist(s(:),r(:),'direction','out','Nbins',[18,18],'color',[0.5,0.5,0.5]);
    axis square;
    set(gca,'Box','On');
    ylabel ('Reported direction (\circ)','FontWeight','Bold','FontSize', 24)
    xlabel ('Stimulus direction (\circ)','FontWeight','Bold','FontSize', 24)
    
    title(t{c_})
    %% to save pic uncomment the following
    %a= sprintf('%s_',t{c_})
    %saveas(gcf,strcat(save_path,a),'png')
end
%% the following code produces the same figures for each individual participant
%% i.e. those used for appendix2-Figure1 in the paper
t = {'zero','mid','high'}
for i_=1:23
for c_ = [3,2,1]
    figure('Renderer', 'painters', 'Position', [10 10 1024 760]),
    r = vertcat(all_rep{i_,c_});
    s = vertcat(all_stim{i_,c_});
    dev = angdiff2(s(:),r(:));
    s=scatterhist(s(:),r(:),'direction','out','Nbins',[18,18],'color',[0.5,0.5,0.5]);
    axis square;
    set(gca,'Box','On');
    ylabel ('Reported direction (\circ)','FontWeight','Bold','FontSize', 24)
    xlabel ('Stimulus direction (\circ)','FontWeight','Bold','FontSize', 24)
    title(t{c_})
    a= sprintf('%s_%i',t{c_},i_)
    %% to save pic uncomment the following
    %saveas(gcf,strcat(save_path,a),'png')
end
end
%% Analysis to check excluded trials per subject/coherence
close all
for i_=1:size(all_rep,1)
    for ii_=1:size(all_rep,2)
        nd(i_,ii_) = sum(sum(isnan(all_rep{i_,ii_})));
    end
end