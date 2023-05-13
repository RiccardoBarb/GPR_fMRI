function [time, task] = extract_design_info(logs_dir, subj_id)
% _
% Extract design information from RDM experiment in fMRI
% 
% Version History:
% - RB,     2017: written
% - RB, 26/10/18: converted into function
% - RB, 17/12/18: added catch trial info
% - RB, 20/02/19: modified for two days
% - JS, 16/04/19: rewritten
% - JS, 14/05/19: added button response info
% - RB, 12/08/19: added trial exclusion based on eyetracking
% 
% Function Authors:
% - RB = riccardo.barbieri@bccn-berlin.de
% - JS = joram.soch@bccn-berlin.de


% experiment dimensions
D = 2;                          % number of days
S = D*5;                        % number of sessions
T = 5;                          % revolution time [s]

% load logfiles
load(strcat(logs_dir,subj_id,'_1','/',subj_id,'_mri_day1.mat'));
data1 = data;                   % data from day 1
load(strcat(logs_dir,subj_id,'_2','/',subj_id,'_mri_day2.mat'));
data2 = data;                   % data from day 2
%load eyetracking data from metafile
load('../data/logfiles_experiment/behavior_and_eyetracking/meta_fixation.mat');
fix1=strcmp({meta_fix.subject},strcat(subj_id,'_1'));
fix2=strcmp({meta_fix.subject},strcat(subj_id,'_2'));
%check if we have eyetracking data or if something is to exclude
if sum(fix1)>0 & meta_fix(fix1).is_exclude >0
    exclude1=1;
else 
    exclude1=0;
end
if sum(fix2)>0 & meta_fix(fix2).is_exclude >0
    exclude2=1;
else 
    exclude2=0;
end
% for each session
for j = 1:S
    
    % if day 1
    if j <= S/D
        t = data1.task.nTraining;
        % stimulation and response phase
        time(j).stim.onsets   = data1.task.Timing.StimTiming{t+j}';
        time(j).stim.duration = data1.task.StimDur;
        time(j).resp.onsets   = data1.task.Timing.ResponseTiming{t+j}';
        time(j).resp.duration = data1.task.Timing.ITITiming{t+j}(2)-data1.task.Timing.ResponseTiming{t+j}(1);
        % button press for response
        time(j).button.resp.onsets   = time(j).resp.onsets+data1.response.RTLine{t+j}';
        time(j).button.resp.duration = 0;
        task(j).button.resp.RT       = data1.response.RTLine{t+j}';
        % button press for catch response
        l = 0;
        for k = find(data1.task.Catch{t+j})
            % get response and catch response angles
            resp_start_angle  = data1.response.LineStartAngle{t+j}(k);
            resp_RT           = data1.response.RTLine{t+j}(k);
            catch_start_angle = data1.response.DirectionCatchStart{t+j}(k);
            catch_end_angle   = data1.response.DirectionCatch{t+j}(k);
            catch_RT          = data1.response.RTCatch{t+j}(k);
            % determine direction of the response bar
            catch_end_if_cw   = mod(catch_start_angle + (catch_RT/T)*360, 360); % clock-wise
            catch_end_if_cc   = mod(catch_start_angle - (catch_RT/T)*360, 360); % counter-clock-wise
            if abs(catch_end_angle - catch_end_if_cw) < abs(catch_end_angle - catch_end_if_cc),
                resp_bar_dir = +1; else, resp_bar_dir = -1; end;
            % determine time from response to catch onset
            resp_end_angle  = mod(resp_start_angle + resp_bar_dir * (resp_RT/T)*360, 360);
            resp_catch_time = mod(resp_bar_dir * (catch_start_angle - resp_end_angle), 360)/360 * T;
            % determine catch response onset as total time
            if ~isnan(catch_RT) & (catch_RT>0)
                l = l + 1;
                time(j).button.catch.onsets(l,1) = time(j).resp.onsets(k) + data1.response.RTLine{t+j}(k) + resp_catch_time + catch_RT;
                time(j).button.catch.duration    = 0;
                task(j).button.catch.RT(l,1)     = data1.response.RTCatch{t+j}(k);
            end;
        end;
        % coherence level, actual direction and reported direction
        task(j).stim.catch     = data1.task.Catch{t+j}';
        task(j).stim.coherence = data1.task.coherence{t+j}';
        task(j).stim.direction = data1.response.ShownAngles{t+j}';
        task(j).resp.direction = data1.response.RecordedAngles{t+j}';
        task(j).resp.nofixation=zeros(size(task(j).resp.direction));
        if exclude1>0
        task(j).resp.nofixation(meta_fix(fix1).exclud_trial_id_per_run(j).trial)=1; 
        end
    % if day 2
    else
        t = data2.task.nTraining;
        % stimulation and response phase
        time(j).stim.onsets   = data2.task.Timing.StimTiming{t+j-S/D}';
        time(j).stim.duration = data2.task.StimDur;
        time(j).resp.onsets   = data2.task.Timing.ResponseTiming{t+j-S/D}';
        time(j).resp.duration = data2.task.Timing.ITITiming{t+j-S/D}(2)-data2.task.Timing.ResponseTiming{t+j-S/D}(1);
        % button press for response
        time(j).button.resp.onsets   = time(j).resp.onsets+data2.response.RTLine{t+j-S/D}';
        time(j).button.resp.duration = 0;
        task(j).button.resp.RT       = data2.response.RTLine{t+j-S/D}';
        
        % button press for catch response
        l = 0;
        for k = find(data2.task.Catch{t+j-S/D})
            % get response and catch response angles
            resp_start_angle  = data2.response.LineStartAngle{t+j-S/D}(k);
            resp_RT           = data2.response.RTLine{t+j-S/D}(k);
            catch_start_angle = data2.response.DirectionCatchStart{t+j-S/D}(k);
            catch_end_angle   = data2.response.DirectionCatch{t+j-S/D}(k);
            catch_RT          = data2.response.RTCatch{t+j-S/D}(k);
            % determine direction of the response bar
            catch_end_if_cw   = mod(catch_start_angle + (catch_RT/T)*360, 360); % clock-wise
            catch_end_if_cc   = mod(catch_start_angle - (catch_RT/T)*360, 360); % counter-clock-wise
            if abs(catch_end_angle - catch_end_if_cw) < abs(catch_end_angle - catch_end_if_cc),
                resp_bar_dir = +1; else, resp_bar_dir = -1; end;
            % determine time from response to catch onset
            resp_end_angle  = mod(resp_start_angle + resp_bar_dir * (resp_RT/T)*360, 360);
            resp_catch_time = mod(resp_bar_dir * (catch_start_angle - resp_end_angle), 360)/360 * T;
            % determine catch response onset as total time
            if ~isnan(catch_RT) & (catch_RT>0)
                l = l + 1;
                time(j).button.catch.onsets(l,1) = time(j).resp.onsets(k) + resp_RT + resp_catch_time + catch_RT;
                time(j).button.catch.duration    = 0;
                task(j).button.catch.RT(l,1)     = data2.response.RTCatch{t+j-S/D}(k);
            end;
        end;
        % coherence level, actual direction and reported direction
        task(j).stim.catch     = data2.task.Catch{t+j-S/D}';
        task(j).stim.coherence = data2.task.coherence{t+j-S/D}';
        task(j).stim.direction = data2.response.ShownAngles{t+j-S/D}';
        task(j).resp.direction = data2.response.RecordedAngles{t+j-S/D}';
        task(j).resp.nofixation=zeros(size(task(j).resp.direction));
        if exclude2>0
        task(j).resp.nofixation(meta_fix(fix2).exclud_trial_id_per_run(j-S/D).trial)=1; 
        end
        
    end;
    
end;