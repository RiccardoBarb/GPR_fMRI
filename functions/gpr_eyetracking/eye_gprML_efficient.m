function bal_acc=eye_gprML_efficient(files_dir, label, save_flag)

% This code runs the gaze position-based motion direction reconstruction
% analysis. It's divided in 2 parts: chanel-wise GPR estimation & stimulus reconstruction.
% For the reconstruction both channels (x and y) are used. Every time point
% of interest is estimated seperately

% Written by Felix M. Toepfer on the basis of a GPR tutorial script for fMRI
% written by Joram Soch & Riccardo Barbieri with contributions from Felix M. Toepfer & Carsten Bogler.

% plot intermediate steps for visualization (1)
% do not plot intermediate steps for visualization (0)
plotit = 0;

%-------------------------------------------------------------------------%
% %% Part 1: Data preparation
% 
% % For the eye-based analysis the independent variable xm is the motion
% % direction shown to the participants or the direction of motion that
% % participants reported to perceived. The dependent variable is 
% % the distance from the gaze to the presented fixation point on the
% % screen. This distance is resolved by the x and y ordinates of the gaze.
% % Each ordinate has its own response profile dependent on the direction of
% % motion xm. The response profiles might change for the time t. Therefore
% % the response profiles are seperately estimated for each time point of
% % interest toi.

% 
% %------------------------% Data properties %------------------------------%


% specify data dimensions

t = 16;                         % number of trials per condition (16 trials per coherence)
v = 2;                          % number of channels (=ordinates)
S = 10;                         % number of runs
VpSL = 2;                       % number of channels per searchlight

coh=3;                          % number of coherences


toi=[501:20:2500, 2500];                    % define time points of interest. traces consist of 500ms prestimulus baseline + 2000 ms stimulus 


prec = 1;                                   %precision of x0 samplig (in degrees)

x0 = [(prec/180)*pi:(prec/180)*pi:2*pi]';   %x0 is the space of independent
                                            %variable i.e. all of the possible
                                            %motion directions


%% Load data

% loop trough all subjects of interest ( all files of interest)

for subj_=1:length(files_dir)

    % preallocate
    bal_acc=nan(coh,numel(toi));
    pred_direct=nan(coh,numel(toi),S,t);
    true_direct=nan(coh,numel(toi),S,t);
    
    % load file
    load(files_dir{subj_})
 
    % print to commandline
    fprintf('subject %i \n',subj_)
%% baseline correction

    eye.EYE(1,:,:) = ft_preproc_baselinecorrect(squeeze(eye.EYE(1,:,:))', 1, 500)';
    eye.EYE(2,:,:) = ft_preproc_baselinecorrect(squeeze(eye.EYE(2,:,:))', 1, 500)';

%% obtain data structure seperately for each coherence level

for coh_=1:coh

    % get data and labels  
    
    % initialize data matrix
    data=nan(size(eye.EYE(:,:,eye.coherence==eye.coh_list(coh_))));
    % get data
    data(:,:,:)=eye.EYE(:,:,eye.coherence==eye.coh_list(coh_)); % channel:tp:all_trials

    % labels
    if strcmp(label,'stimulus')
        labels=eye.stim(eye.coherence==eye.coh_list(coh_));

    elseif strcmp(label,'report')
        labels=eye.response(eye.coherence==eye.coh_list(coh_));

    else
        error('Wrong label string. Use ''stimulus'' or ''report''.')
    end
    
    % exlcude trials where subject did not fixate. This meta information
    % is obtained from the fixation controll pipeline. A trial was not
    % analyzed in the fMRI data analysis and will not be incorporated here
    % if the gaze of a subject was 2 dva away from fixation for more than
    % 200 ms.
    % For excluded trials data and labels are exchanged to NaNs    
    nfix=eye.not_fixate(eye.coherence==eye.coh_list(coh_));
    
    data(:,:,logical(nfix))=nan;
    labels(:,logical(nfix))=nan;

    % plot x channel gaze position over label for every timepoint of interest
    if plotit
        for t_=numel(toi)
            figure(t_)
            plot(labels,squeeze(data(1,toi(t_),:)),'*')
        end
    end
    %reshape to separate runs for later crossvalidation
    data= reshape(data(:,:,:,:),[2,2500,16,10]);
    labels=reshape(labels, [16, 10]);
    
    % get timepoints of interest
    data=(data(:,toi,:,:));

    % permute to fit structure of analysis
    yt=permute(data,[3,1,2,4]); % trials x channel x timepoints x runs(aka crossvalidationfolds)

    if plotit
        %plot at the last timepoint
        plot_raw(squeeze(yt(:,:,1,:)),labels)
    end
    
    
%% Part 2: Perform GPR estimation

    % The goal of this part is to obtain a continuous estimate of the 
    % response profile (Y) as a function of motion direction (x).

    % The model follows the assumption that the gaze position in a certain trial is informative
    % about the motion direction, displaying a smooth response profile.
    % However, our experimental design will only include a
    % limited number of sampled directions and gaze ordinates.
    % GPR fills the gaps by exploiting the covariance matrix (kernel) of
    % each ordinate position graph across trials (see covPeriodic_RDM for details on the covariance function).
    % This allow to obtain an estimate of the ordinates's mean response profile (mup),
    % its variance s2x, and the kernel's hyperparameters (hyp), which will be
    % then used to perform the motion direction or choice reconstruction.

    %------------------------% Analysis setup %-------------------------------%
    
    % pre-allocate

    m  = numel(x0);                     % number of possible motion directions

    %preallocate outcome variables
    mup = zeros(m,v,numel(toi),S);                 % predicted continuous ordinate response profile // 360 (=circle with precision 1 deg):channels:runs
    mux = zeros(t*(S-1),v,numel(toi),S);           % fitted trial-wise voxel response // trials*runs-1 : voxel : runs
    s2x = zeros(S,v,numel(toi));                   % voxels variance
    hyp(S,v,numel(toi)) = struct('mean',...
    [], 'cov', [0], 'lik', -1);     % voxels Kernel hyperparameters
    
    
    for chan=1:v

        fprintf('\n channel %i \n',chan)

        %------------------------% Channel-wise GPR %-------------------------------%
       
        % estimate separately for each timepoint of interest
        parfor ii = 1:numel(toi)

            vv = squeeze(yt(:,chan,ii,:)); % data= all trials: of 1 channel : all runs : ONLY 1 tp

            for r = 1:S

                %training index for current cv fold
                itrain =  [1:S] ~= r;
                
                %vectorize data for current cv fold
                vectors_train = vv(:,itrain)';
                Samples = vectors_train(:);

                labels_train = labels(:,itrain)';
                Labels = labels_train(:)./360*2*pi;     %Labels are converted in rad

                % cyclic multivariate Gaussian process regression

                [mup(:,chan,ii,r), mux(:,chan,ii,r), s2x(r,chan,ii), hyp(r,chan,ii)] = ME_cmGPR_mean(Samples, Labels,x0);

            end    
        end
    end
    %% Part 3: Perform searchlight based stimulus reconstruction

    % Now we can use the response profiles (mup) to perform direction
    % reconstruction via MLE. We first combine the response profiles of
    % each ordinate (x and y) within a searchlight (mup_SL).
    % We then estimate the searchlight spatial
    % covariance (S2_est - see ME_cmGPR_cov).The important difference between
    % the fMRI analysis is here, that we do not need to apply regularisation 
    % when estimating the spatial covariance, because we have less response
    % profiles than single trials. Then we compute the Log PDF of the
    % multivariate normal distribution with means = mup_SL and covariance = S2_est
    % (see logmvnpdf). The predicted motion direction will be the one maximizing
    % the log-liklehood given the current data and response profiles in a SL.
    % This procedure is performed in a run-wise cv procedure. Finally we calculate
    % the absolute angular deviation between the true and the predicted direction
    % to compute a crossvalidated continuous accuracy measure (see
    % avg_norm_circ_resp_dev and bal_norm_circ_resp_dev).

    %------------------------% Analysis setup %-------------------------------%
    % define searchlights

    % eye tracking: 1 searchlight with 2 voxels
    SL = zeros(v-VpSL+1,VpSL);  
    for k = 1:size(SL,1)
        SL(k,:) = [k:(k+VpSL-1)];
    end
    
        fprintf('\n\n-> SL-based reconstruction:\n');
        
    %preallocate outcome variables
    for tp=1:numel(toi)

        acc_vol = zeros(size(SL,1), 1);           % SL cross-validated accuracy
        bal_acc_vol = zeros(size(SL,1), 1);       % SL balanced cross-validated accuracy
        predicted_directions = cell(size(SL,1),S);% SL trial-wise predicted direction
        true_directions = cell(S,1);              % tesing lables of each cv-fold

        %-----------------% SL-based stimulus reconstruction %--------------------%


            for ii = 1:size(SL,1)

                % Bring only data from the current searchlight into the loop

                vv_sl = squeeze(yt(:,SL(ii,:),tp,:));       % trial-wise ordinates
                hyp_sl = squeeze(hyp(:,SL(ii,:),tp));       % voxel-wise kernel hyperparameters
                mup_sl = squeeze(mup(:,SL(ii,:),tp,:));     % predicted continuous ordinate position
                mux_sl = squeeze(mux(:,SL(ii,:),tp,:));     % fitted trial-wise ordinate position
                s2x_sl = squeeze(s2x(:,SL(ii,:),tp));       % ordinate variance
                n_sl_vox = numel(SL(ii,:));     % number of channels in the sl
                
        %force real in case of complex values
            if ~isreal(mup_sl) || ~isreal(mux_sl) || ~isreal(s2x_sl) 
                 mup_sl = real(squeeze(mup(:,SL(ii,:),tp,:)));     
                 mux_sl = real(squeeze(mux(:,SL(ii,:),tp,:)));     
                 s2x_sl = real(squeeze(s2x(:,SL(ii,:),tp)));       
            end               
                
                
                for r = 1:S % crossvalidation

                    %training and testing index for current cv fold
                    itrain = find( [1:S] ~= r);
                    itest  = r;

                    %vectorize training data for current cv fold
                    vectors_train = vv_sl(:,:,itrain);
                    train_samples = permute(vectors_train,[3,1,2]);
                    train_samples = reshape(train_samples,[t*(S-1),n_sl_vox]);

                    labels_train = labels(:,itrain)';
                    labels_train = labels_train(:)./360*2*pi;     %Labels are converted in rad

                    %vectorize testing data for current cv fold

                    test_samples = vv_sl(:,:,itest);
                    labels_test = labels(:,itest);
                    labels_test = labels_test./360*2*pi;     %Labels are converted in rad

                    % cyclic multivariate Gaussian process regression

                    % estimate spatial covariance in sl
                    % important difference to the fMRI analysis pipeline is
                    % that for the eye tracking data analysis we have less
                    % channels than trials, therefore we do not need to
                    % regularize die covariance matrix
                    
                    S2_est = ME_cmGPR_cov2(train_samples, mux_sl(:,:,r), s2x_sl(r,:), 'steplog');

                    % reconstruct direction from sl activity
                    predicted_directions{ii,r} = ME_cmGPR_pred(test_samples, x0, mup_sl(:,:,r), S2_est);
                    true_directions{r} = labels_test;

                end

                % quantify crossvalidated accuracy and balanced crossvalidated
                % accuracy
                current_predictions = [predicted_directions{ii,1:end}]';
                
                pred_direct(coh_,tp,:,:)=current_predictions;

                
                current_true_dir =[true_directions{1:end}]';
                true_direct(coh_,tp,:,:)= current_true_dir;
                
                acc_vol(ii) = avg_norm_circ_resp_dev(current_predictions(:),current_true_dir(:));
                bal_acc_vol(ii)= bal_norm_circ_resp_dev(current_predictions(:),current_true_dir(:),'trapz');

                acc_vol(ii) = acc_vol(ii)*100;
                bal_acc_vol(ii) =  bal_acc_vol(ii)*100;

                fprintf('successful!\n');

            end

            if plotit
                plot_fit(yt,labels,mup,1,[1:v],'2d')
                plot_fit(yt,labels,mup,1,[1:v],'3d')
                dx = plot_reconstruction(predicted_directions,true_directions,10);
            end

        bal_acc(coh_,tp)=bal_acc_vol;
    end
end

if save_flag
    savefolder=[files_dir{subj_}(1:end-12) 'gpr'];
    if ~exist(savefolder,'dir')
        mkdir(savefolder)
    end
    
    savename=[savefolder filesep 'gpr_results_' label '.mat'];
    save(savename,'bal_acc','true_direct','pred_direct')
end

end
end
