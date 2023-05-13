function [timing]= get_traces(dat, event, behavior)


        fprintf('separate interesting periods\n')

        type={event.type};
        run_mess=find(strcmp(type,'Run:')==1);


        value={event.value};
        a=strfind(value,'Trail:');
        trial_message=find(cellfun(@isempty,a)==0);

        % devide into runs and trials

        trial_mess_last=find(diff(trial_message)>1);
        trial_mess_first= [1 trial_mess_last+1];
        trial_mess=[trial_mess_first',[trial_mess_last length(trial_message)]'];


        trial_count=diff(trial_mess')+1;

        m=0;

        for n=1:length(trial_mess)

        % timing(n).sample(1)=event(trial_message(trial_mess(n,1))).sample;
        % timing(n).sample(2)=event(trial_message(trial_mess(n,2))+1).sample;

        timing(n).onset=event(trial_message(trial_mess(n,1))).sample;
        timing(n).offset=event(trial_message(trial_mess(n,2))+1).sample;

        % timing(n).gazex=dat.gazex(timing(n).sample(1):timing(n).sample(2));
        % timing(n).gazey=dat.gazey(timing(n).sample(1):timing(n).sample(2));

        %timing(n).gazex=dat.gazex(timing(n).onset:timing(n).offset);
        %timing(n).gazey=dat.gazey(timing(n).onset:timing(n).offset);

        timing(n).center(1)=mean(dat.gazex(timing(n).onset:timing(n).offset));
        timing(n).center(2)=mean(dat.gazey(timing(n).onset:timing(n).offset));


        for t=1:trial_count(n)      % get sample start and end for every trial
            m=m+1;
            if t==trial_count(n) % last trial needs and other referenzing

        %         timing(n).trial(t).sample(1)=event(trial_message(m)).sample;
        %         timing(n).trial(t).sample(2)=event(trial_message(m)+1).sample;
        %         timing(n).trial(t).gazex=dat.gazex(timing(n).trial(t).sample(1):timing(n).trial(t).sample(2));
        %         timing(n).trial(t).gazey=dat.gazey(timing(n).trial(t).sample(1):timing(n).trial(t).sample(2));

                timing(n).trial(t).onset=event(trial_message(m)).sample;
                timing(n).trial(t).offset=event(trial_message(m)+1).sample;
        %        timing(n).trial(t).gazex=dat.gazex(timing(n).trial(t).onset:timing(n).trial(t).offset);
        %        timing(n).trial(t).gazey=dat.gazey(timing(n).trial(t).onset:timing(n).trial(t).offset);


            else    
        %         timing(n).trial(t).sample(1)=event(trial_message(m)).sample;
        %         timing(n).trial(t).sample(2)=event(trial_message(m+1)).sample;
        %         timing(n).trial(t).gazex=dat.gazex(timing(n).trial(t).sample(1):timing(n).trial(t).sample(2));
        %         timing(n).trial(t).gazey=dat.gazey(timing(n).trial(t).sample(1):timing(n).trial(t).sample(2));

                timing(n).trial(t).onset=event(trial_message(m)).sample;
                timing(n).trial(t).offset=event(trial_message(m+1)).sample;
        %        timing(n).trial(t).gazex=dat.gazex(timing(n).trial(t).onset:timing(n).trial(t).offset);
        %        timing(n).trial(t).gazey=dat.gazey(timing(n).trial(t).onset:timing(n).trial(t).offset);

            end
        end

        end


        % seperate in periodes per trial
        % we can not simply use the Run: x, Trial: y information in the event file,
        % because each trial starts with the intertrial interval followed by the
        % stimulation and so on. the event information therefor co-incides with the
        % appearance of the intertrialinterval (that has a differnet duration every
        % trial).

        % iti
        for n=1:length(timing)
        a=num2cell([timing(n).trial(:).onset])';
        [timing(n).trial(:).iti_onset]=deal(a{:});

        b=num2cell([timing(n).trial(:).onset]+(round(behavior(n).iti_dur*1000)))';
        [timing(n).trial(:).iti_offset]=deal(b{:});

        % stim
        c=num2cell([timing(n).trial(:).iti_offset])';
        [timing(n).trial(:).stim_onset]=deal(c{:});

        d=num2cell([timing(n).trial(:).stim_onset]+(round(behavior(n).stim_dur*1000)))';
        [timing(n).trial(:).stim_offset]=deal(d{:});

        % resp

        e=num2cell([timing(n).trial(:).stim_offset])';
        [timing(n).trial(:).resp_onset]=deal(e{:});

        f=num2cell([timing(n).trial(:).offset])';
        [timing(n).trial(:).resp_offset]=deal(f{:});

        end