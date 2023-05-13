function [behavior, display]= extract_behavior(behavFile)
%
%   The function extract relevant information from a given behavior file

% in:   - behavFile:    directory of the behavior file that needs to be
%                       analysed
% out:  - behavior:     a structur that contains the relevant information

% tells u what happens    
    fprintf('get behavior data\n')

% load the file
    load(behavFile);

%get number of runs
    nruns=data.task.run;
    
% get number of bins
    nbins = data.task.NumbOfDirection;

% extract data of intrest sorted by runs
    for r_=1:nruns
        behavior(r_).coh=data.task.Coherences;
        behavior(r_).direction=data.task.direction{r_};
        behavior(r_).direction_bin=bin_angles(data.task.direction{r_}, nbins);
        behavior(r_).response=data.response.RecordedAngles{r_};
        behavior(r_).response_bin=bin_angles(data.response.RecordedAngles{r_}, nbins);
        behavior(r_).coherence=data.task.coherence{r_};
        behavior(r_).iti_dur=data.task.Timing.StimOnset{r_}-data.task.Timing.ITIOnset{r_}(1:end-1);
        behavior(r_).stim_dur=data.task.Timing.RespOnset{r_}-data.task.Timing.StimOnset{r_};
        behavior(r_).resp_dur=data.task.Timing.ITIOnset_real{r_}-data.task.Timing.RespOnset{r_};
    end
        
% get the values of the display the stimulation was presented to make sure
% you later on calculate the fixation radius correctly
    display=data.display;

end  
