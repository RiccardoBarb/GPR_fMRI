%% This function is written to plot the result of the gpr based 
%% reconstruction analysis of the eye tracking data obtained in the realted 
%% fMRI project.

% This script is a merge between a script by Gi-Yeul Bae
% to assess significant clusters in time continuous data
% Bae, Gi-Yeul, and Steven J. Luck. "Decoding motion direction using the topography of sustained ERPs and alpha oscillations." NeuroImage 184 (2019): 242-255.
% And a Charlie M Sexton 
% All credits go to these authors

% Written Felix M Töpfer 15112020
% last edit 26012021
% N.B this script must be run after running
% pipeline_for_estimation_and_reconstruction_from_eyetracking.m and after
% running permutation_analysis_gpr_eyetracking.m


%% set paths


clear; clc; close all;

% add the current path, because this is your working directory
thispath=pwd;
boundedline_path =%necessary to plot lines with errobars https://www.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m
addpath(genpath(thispath));
addpath(genpath('../functions'));
addpath(genpath(boundedline_path));


%% define variables 

rng('shuffle')

label='stimulus';                 % define labels 'stimulus' or 'report'
coh = {'low','mid','high'};     % define coherence
tp = 101;                       % number of timepoints of interest


save_folder='../../data/eye_tracking_for_GPR/results';
save_figures = 1; % save outputs to folder
save_files =1;
plot_individual = 0; % plot individual accuracies

%% get directory of interesting files

  % get paths
    data_dir='/analysis/condec/data/222_FMRI/eye_tracker';
    search_index='s3';
   
    folder_dir=get_data_folder_directory(data_dir,search_index);


    for f_=1:numel(folder_dir)
        data_direct=folder_dir{f_};
        search_index='gpr';
        folder_dir(f_)=get_data_folder_directory(data_direct,search_index);
    end
    
    search_str=label;
    file_extension='.mat';
    [files_dir, subj_name]=get_files_directory(folder_dir,search_str,file_extension);
    
%% get data

    % preallocate
        subj_accu=nan(numel(files_dir), numel(coh), tp);
        
    % load data and smooth
        for subj_=1:numel(files_dir)
            fprintf('\nload subject %i\n',subj_)
            load(files_dir{subj_}) 
            subj_accu(subj_,:,:)=bal_acc/100;
        end
    
    % estimate separately for each coherence level
    
    for coh_=1:numel(coh)

        % smooth subject data
            subj_smo=nan(numel(files_dir),tp);
            for subj_=1:numel(files_dir)
                subj_smo(subj_,:)=subj_accu(subj_,coh_,:);      
            end

        % estimate group mean and error
            group_mean = squeeze(nanmean(subj_smo,1));
            group_se = squeeze(nanstd(subj_smo,1))/sqrt(numel(files_dir));

        %% mass T cluster analysis

            releventTime = 1:tp; % all timepoints

            Ps = nan(2,length(releventTime));
            for i = 1:length(releventTime) % make sure this time range is correct
                tp = releventTime(i);

                [H,P,CI,STATS] =  ttest(subj_smo(:,tp), 0.5,'tail','right'); % Test Against chance (=50%)

                Ps(1,i) = STATS.tstat;
                Ps(2,i) = P;
            end
            % find significant points
            candid = Ps(2,:) <= .05;
            if any(candid)
                disp(['Significant at: ' coh{coh_}]);
            end

            candid_marked = zeros(1,length(candid));
            candid_marked(1,1) = candid(1,1);
            candid_marked(1,length(candid)) = candid(1,length(candid));
            %remove orphan time points
            for i = 2:length(releventTime)-1

                if candid(1,i-1) == 0 && candid(1,i) ==1 && candid(1,i+1) ==0
                    candid_marked(1,i) = 0;
                else
                    candid_marked(1,i) = candid(1,i);
                end

            end

            % combine whole time range with relevent time & significant information
            clusters = zeros(length(tp),1);
            clusterT = zeros(length(tp),1);
            clusters(releventTime,1) = candid_marked;
            clusterT(releventTime,1) = Ps(1,:);
            clusterTsum = sum(Ps(1,logical(candid_marked)));

            %%find how many clusters are there, and compute summed T of each cluster
            tmp = zeros(10,25);
            cl = 0;
            member = 0;
            for i = 2:length(clusters)-1


                if clusters(i-1) ==0 && clusters(i) == 1 && clusters(i+1) == 1
                    cl = cl+1;
                    member = member +1;
                    tmp(cl,member) = i;

                elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 0
                    member = member +1;
                    if cl == 0
                        cl = cl+1; % in case first index in clusters equals 1
                    end
                    tmp(cl,member) = i;
                    member = 0;
                elseif clusters(i-1) ==1 && clusters(i) == 1 && clusters(i+1) == 1
                    member = member +1;
                    if cl == 0
                        cl = cl+1; % in case first index in clusters equals 1
                    end
                    tmp(cl,member) = i;

                else


                end
            end


            HowManyClusters = cl;
            a = tmp(1:cl,:);
            eachCluster = a(:,logical(sum(a,1)~=0));

            %now, compute summed T of each cluster
            dat_clusterSumT = zeros(HowManyClusters,1);
            for c = 1:HowManyClusters
                dat_clusterSumT(c,1) = sum(clusterT(eachCluster(c,eachCluster(c,:) ~=0)));
            end
    
        %% find .05 point - need to repeat with new permutation files for decisions
            if strcmp(label,'report')
               % load([data_folder filesep 'simulation_100_interation_press' filesep 'simulation_' data_to_plot '_' period '_' coh '_Decisions.mat']) % load different null distribution of t-mass.
                load([data_dir filesep 'permutation_test' filesep 'GPR_simulation' '_' coh{coh_} '_report.mat']) % load different null distribution of t-mass.
            elseif strcmp(label,'stimulus')
                load([data_dir filesep 'permutation_test' filesep 'GPR_simulation' '_' coh{coh_} '_stimulus.mat']) % load different null distribution of t-mass.
            end
            iteration = 1000;
            cutOff = iteration - iteration * 0.05; %one tailed
            % sortedTvlaues = sort(EmpclusterTvalue,2);
            critT = simulationT(cutOff); % 2 tailed iteration * 0.025
            sigCluster = dat_clusterSumT > critT;
            draw = eachCluster(sigCluster,:); % CMS - this is used to shade areas that are significant - not currently using it though
            draw = sort(reshape(draw,1,size(draw,1)*size(draw,2)));
            draw = draw(draw>0);

            for si = 1 :size(dat_clusterSumT,1)

                [~,where] = min(abs(simulationT - dat_clusterSumT(si)));
                pv = 1 - where/iteration;
            end

            % cS: save info relevant to plotting all coherences together in
            % separate structure
            allStats(coh_).Coh=coh;
            allStats(coh_).critT = critT;
            allStats(coh_).Ps = Ps;
            allStats(coh_).UsedSubj=numel(files_dir);
            allStats(coh_).excludedSubj={};
            allStats(coh_).draw = draw;
            allStats(coh_).subAverage = group_mean;
            allStats(coh_).seAverage = group_se;
            allStats(coh_).eachCluster = eachCluster(sigCluster,:); % makes draw sig lines easier    


    end
    
%% plot

%% plot significant clusters

%     cl=colormap(parula(50));
cl = get(gca,'ColorOrder');

figure

% % add clusters
hold on
for c = 1:length(coh)
    subAverage = allStats(c).subAverage;
    seAverage = allStats(c).seAverage;
    draw = allStats(c).draw;
    eachCluster = allStats(c).eachCluster;
    
    accEst = squeeze(subAverage);
    w = zeros(tp,1)';
    w(draw)=1;
    
    %         a = area(1:length(tm), accEst.*w,'HandleVisibility','off');
    %         a.EdgeColor = 'none';
    %         a.FaceColor = [0.8,0.8,0.8];
    %         child = get(a,'Children');
    %         set(child,'FaceAlpha',0.9)
    
    % this draw lines over sections above chance instead of shaded area
    % as above. will not work with multiple sig areas atm
    %         if any(draw)
    %             line([draw(1) draw(end)],[.195 .195],'LineWidth',4,'Color',cl(c,:));
    %         end
    
    for cc = 1:size(eachCluster,1)
        cluster = eachCluster(cc,:); cluster = cluster(cluster~=0);
        l(cc) = line([cluster(1) cluster(end)],[.58 .58]-.003*c,'LineWidth',4);
        l(cc).Color = [cl(c,:) 0.5];
    end
    
    hold on
    %         [lines(c) patches(c)] = boundedline(1:length(tm),subAverage,seAverage, 'cmap',cl(50-16*c,:),'alpha','transparency',0.3);
    %         plot(1:length(tm),subAverage,'LineWidth',1.5,'Color',cl(50-16*c,:)); % make lines bolder
    [lines(c) patches(c)] = boundedline(1:tp,subAverage,seAverage,'cmap',cl(c,:),'alpha','transparency',0.2);
    plot(1:tp,subAverage,'LineWidth',1.5,'Color',cl(c,:)); % make lines bolder
end
xlabel('Time (ms)', 'Fontsize', 20);
ylabel('Decoding Accuracy', 'Fontsize', 20)
ax = gca;
ax.YLim = [.45, .6];

    ax.XTick = [1 26 51 76 100]; % set xTick locaiton. Colum Location in the image matrix

    ax.XTickLabel = {'0','500','1000','1500','2000'};
    ax.XLim = [1 125];
    ax.FontSize=16;

h = line(1:tp,0.5* ones(1,tp),'HandleVisibility','off');
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
[~,hObj] = legend(lines, coh,'Location','NorthEast');
hL=findobj(hObj,'type','line');
set(hL,'linewidth',4)      
legend boxoff;
title([' gpr reconstruction of ' label ])
set(gcf,'color','w');
set(gca,'FontSize',14);
 
if save_figures
    fname = [save_folder filesep 'GPR_' label '_reconstruction_accuracy.png'];
    saveas(gcf,fname)
end
if save_files

    name=['GPR_stats_'  label '.mat'];
    fname= fullfile(save_folder,name);
    save(fname,'allStats');
end
