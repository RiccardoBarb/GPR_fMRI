function [dat, event] = prepro_eye(ascFile)

% this function is written to preprocess eye link specific eye tracking
% data. It is mainly written by Anne Urai and only adapted by Felix Töpfer.
% The following functions are all created by Anne Urai '18

% in:   - ascFile   : directory to the *.asc file that should be analysed

% out:  - dat       : preprocessed and blink interpolated data
%       - event     : events that were written in the eyelink recorded data

% Felix 170719

        % ============================================== %
        % 2. create a FieldTrip-style data structure
        % ============================================== %
        % read in the asc EyeLink file
        
        asc = read_eyelink_ascNK_AU(ascFile);

        
        % create events and data structure, parse asc
        
        [data, event, blinksmp, saccsmp] = asc2dat(asc);
        
        
        % ============================================== %
        % 3. interpolate Eyelink-defined and additionally detected blinks
        % ============================================== %
       
        plotMe = false; % possible to plot the outcome if set to: true
              
        [newpupil, newblinksmp, nanIdx, dat] = blink_interpolate(data, blinksmp, plotMe);
        
        
        dat.saccsmp=saccsmp;
        
        