function [timing] = control_fixation (dat, display, timing, radius, duration)

        fprintf('fixation control\n')
        
%% get fixation radius

rad=angle2pix(display, radius);

%% get fixation information for every whole trial


% NaN out the interpolated data, because this obscures the fixation
% detection

    x=dat.gazex;
    y=dat.gazey;

    for n=1:length(dat.interpol)
        x(dat.interpol(n,1):dat.interpol(n,2))=NaN;
        y(dat.interpol(n,1):dat.interpol(n,2))=NaN;
    end
    
    for n=1:length(dat.saccsmp)
        x(dat.saccsmp(n,1):dat.saccsmp(n,2))=NaN;
        y(dat.saccsmp(n,1):dat.saccsmp(n,2))=NaN;
    end
    

        for r_=1:length(timing)
            


%            % now analyse fixation but on the traces without interpolation
%             
%             for t_=1:length(timing(r_).trial)
%                 gazex_trial=x(timing(r_).trial(t_).onset:timing(r_).trial(t_).offset);
%                 gazey_trial=y(timing(r_).trial(t_).onset:timing(r_).trial(t_).offset);
% 
%                 figure(r_);hold on; plot(gazex_trial,gazey_trial,'k')
%                 figure(r_);hold on; plot_circle(timing(r_).center(1),timing(r_).center(2), 134);
% 
%                 diff_fix=pdist2([gazex_trial',gazey_trial'],[timing(r_).center(1),timing(r_).center(2)]);
% 
%                 exceed=find(diff_fix>rad);
%                 figure(r_);hold on; plot(gazex_trial(exceed),gazey_trial(exceed),'c')
%                 exceed_dif=[1 (find(diff(exceed)>1))'];
%                 dur=diff(exceed_dif);
%                 id_diff=find(dur>duration);
%                 timing(r_).trial(t_).exceed=length(id_diff);
%             end
% 
% 
%             timing(r_).exceed=sum([timing(r_).trial(:).exceed]);
%             timing(r_).exceeded_trial_id=find([timing(r_).trial(:).exceed]>0);
        end

% analyse fixation ONLY during STIMULUS 

% get the number of trials per timing


        for r_=1:length(timing)
            
            
            % check the noisyness of the data, we do that with observing the velocity
            % values on the interpolated data. the interpolation smooths the data in a
            % way because the high velocities were already excluded. If the data is
            % very noise there are still huge varities in the velocity and we want to
            % exclude these timings
% 
%                 timing(r_).noise_level(1)=std(diff(dat.gazex(timing(r_).onset:timing(r_).offset)));
%                 timing(r_).noise_level(2)=std(diff(dat.gazey(timing(r_).onset:timing(r_).offset)));

%                 x_noise=x(timing(r_).onset:timing(r_).offset); x_noise(isnan(x_noise))=[];
%                 y_noise=y(timing(r_).onset:timing(r_).offset); y_noise(isnan(y_noise))=[];
%                 timing(r_).noise_level(1)=std(diff(x_noise))/sqrt(length(x_noise));
%                 timing(r_).noise_level(2)=std(diff(y_noise))/sqrt(length(y_noise));

            x_timing=x(timing(r_).onset:timing(r_).offset);
            y_timing=y(timing(r_).onset:timing(r_).offset);
            
            % define standart deviation of real world eye trace
            diff_fix=pdist2([x_timing',y_timing'],[timing(r_).center(1),timing(r_).center(2)]);
            timing(r_).std=nanstd(diff_fix);
            
            
            for t_=1:length(timing(r_).trial)
                gazex_stim=x(timing(r_).trial(t_).stim_onset:timing(r_).trial(t_).stim_offset);
                gazey_stim=y(timing(r_).trial(t_).stim_onset:timing(r_).trial(t_).stim_offset);

        %        figure(r_);hold on; plot(gazex_stim,gazey_stim,'b')
                % get euclidian distance from fixation to gaze position
                diff_fix=pdist2([gazex_stim',gazey_stim'],[timing(r_).center(1),timing(r_).center(2)]);
                % find gazeposition that exceed the fixation radius, one
                % need to ad a value at the end otherwise it misses the
                % values at the end
                exceed=[find(diff_fix>rad); length(gazex_stim)*2];
                % define the periods where the gaze exceeds the fixations
                exceed_dif=[1 (find(diff(exceed)>1))' ];
                dur=diff(exceed_dif);
                
                gazex_stim(isnan(gazex_stim))=[];
                
                %of there is a fixation that is larger than 200 ms over
                %2dva and the trace during the stimulus presentation is
                %longer than 1000ms (we want to analyse at least half of the trace) the trial needs to
                %be excluded
               if find(dur>200) & (length(gazex_stim)>1000)
                    timing(r_).trial(t_).stim_exceed=1;
               else
                    timing(r_).trial(t_).stim_exceed=0;
               end
            end
            timing(r_).stim_exceeded=sum([timing(r_).trial(:).stim_exceed]);
            timing(r_).stim_exceeded_trial_id=find([timing(r_).trial(:).stim_exceed]>0);
        end
        
%% get fixation information for every whole timing

        % for r_=1:length(timing)
        % %diff_fix=abs(timing(n).gazey-timing(n).center(2));
        %         gazex_timing=dat.gazex(timing(r_).onset:timing(r_).offset);
        %         gazey_timing=dat.gazey(timing(r_).onset:timing(r_).offset);
        %         
        % %         figure(r_);hold on; plot(gazex_timing,gazey_timing,'k')
        % %         figure(r_);hold on; plot(timing(r_).center(1),timing(r_).center(2),'k+')
        % %         figure(r_);hold on; plot_circle(timing(r_).center(1),timing(r_).center(2), 134);
        % %         
        %         
        %         diff_fix=pdist2([gazex_timing',gazey_timing'],[timing(r_).center(1),timing(r_).center(2)]);
        % 
        %         exceed=find(diff_fix>134);
        %         
        % %         figure(r_);hold on; plot(gazex_timing(exceed),gazey_timing(exceed),'c')
        %         
        %         exceed_dif=[1 (find(diff(exceed)>1))'];
        %         dur=diff(exceed_dif);
        %         id_diff=find(dur>200);
        % 
        %         timing(r_).exceed=length(id_diff);
        % 
        % end