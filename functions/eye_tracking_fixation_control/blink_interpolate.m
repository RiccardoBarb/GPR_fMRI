function [newpupil, newblinksmp, nanIdx, dat] = blink_interpolate(data, blinksmp, plotme)
% interpolates blinks and missing data
% Anne Urai, 2016

% get the stuff we need
dat.time        = data.time{1};
dat.pupil       = data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyepupil')),:);
dat.gazex       = data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyeh')),:);
dat.gazey       = data.trial{1}(~cellfun(@isempty, strfind(lower(data.label), 'eyev')),:);
padding         = 0.2; % how long before and after do we want to pad?

dat.blinksmp    = blinksmp; % get sample idx from asc

% initialize settings
if ~exist('plotme', 'var'); plotme = true; end % plot all this stuff

if ~isempty(blinksmp),
    disp('interpolating EL-defined blinks');
    
    % ====================================================== %
    % STEP 1: INTERPOLATE EL-DEFINED BLINKS
    % ====================================================== %
    
    if plotme,
        clf;  sp1 = subplot(411); plot(dat.time,dat.pupil, 'color', [0.5 0.5 0.5]);
        axis tight; box off; ylabel('Raw');
        set(gca, 'xtick', []);
    end
    
    % merge 2 blinks into 1 if they are < 250 ms together (coalesce)
    coalesce = 0.250;
    for b = 1:size(blinksmp, 1)-1,
        if blinksmp(b+1, 1) - blinksmp(b, 2) < coalesce * data.fsample,
            blinksmp(b, 2) = blinksmp(b+1, 2);
            blinksmp(b+1, :) = nan;
        end
    end
    % remove those duplicates
    blinksmp(isnan(nanmean(blinksmp, 2)), :) = [];
    
    % pad the blinks
    padblinksmp(:,1) = round(blinksmp(:,1) - padding * data.fsample);
    padblinksmp(:,2) = round(blinksmp(:,2) + (padding+0.2) * data.fsample);
    
    % avoid idx outside range
    if any(padblinksmp(:) < 1), padblinksmp(find(padblinksmp < 1)) = 1; end
    if any(padblinksmp(:) > length(dat.pupil)), padblinksmp(find(padblinksmp > length(dat.pupil))) = length(dat.pupil); end
    
    % make the pupil NaN at those points
    for b = 1:size(padblinksmp,1),
        dat.pupil(padblinksmp(b,1):padblinksmp(b,2)) = NaN;
        dat.gazex(padblinksmp(b,1):padblinksmp(b,2)) = NaN;
    end
    
        % felix020819: we do not use this part of the analysis. because it
        % actually does not define missing data, but data with peak in
        % difference over 4*std. This does work for Anne, because in an
        % pupil experiment you take a lot of care about the luminance
        % stability during the stimulus presentation. We havnt done this.
        % this leads to stronger changes in the pupil diameter on each
        % stimulus onset ( as here is higher luminance ). This in turn
        % leads to the fact that gaze data is erased even though just hte
        % pupil changed more strongly ( please remember that here the
        % threshold definded blinks are definded using the pupil not the
        % gaze!) Anyways, comparing the gaze and pupil defined additional
        % blinks the pupil seems to work more consistent without detecting
        % larger saccades as blinks. So we use the code from Anne as is,
        % but only dont to this additional thresholding
        %     % also set the pupil to zero when there were missing data
        %     dat.pupil(dat.pupil < nanmedian(dat.pupil)-4*nanstd(dat.pupil)) = nan;
        %     dat.pupil(dat.pupil > nanmedian(dat.pupil)+4*nanstd(dat.pupil)) = nan;
        %     
        %     % store these points
        %     a=dat.pupil < nanmedian(dat.pupil)-4*nanstd(dat.pupil);a=find(a==1);
        %     b=dat.pupil > nanmedian(dat.pupil)+4*nanstd(dat.pupil);b=find(b==1);
        %     c=[a b]; c=sort(c);
        %     d=[0 find(diff([c max(c)*2])>1)];
        %     
        %     for n=1:(length(d)-1)
        %        nandat(n,1)= c(d(n)+1);
        %        nandat(n,2)= c(d(n+1));
        %     end

        %    dat.pupil(dat.pupil < 10) = nan;
    
    % save NaN idx
    nanIdx = find(isnan(dat.pupil));
    
    % interpolate linearly
    dat.gazex(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.gazex(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'linear');
    dat.gazey(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.gazey(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'linear');

    % now the pupil itself
    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'linear');

    % to avoid edge artefacts at the beginning and end of file, pad in seconds
    edgepad = 1;
    dat.pupil(1:edgepad*data.fsample)           = NaN;
    dat.pupil(end-edgepad*data.fsample : end)   = NaN;
    
    dat.gazey(isnan(dat.pupil)) = NaN;
    dat.gazex(isnan(dat.pupil)) = NaN;
    
%    also extrapolate ends
    dat.gazex(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.gazex(~isnan(dat.pupil)), find(isnan(dat.pupil)), ...
        'nearest', 'extrap');
    dat.gazey(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.gazey(~isnan(dat.pupil)), find(isnan(dat.pupil)), ...
        'nearest', 'extrap');
    
    % now the pupil itself
    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), ...
        'nearest', 'extrap');
      
      if plotme, sp2 = subplot(411); hold on;
        % show how well this worked
        plot(dat.time, dat.pupil, 'b');
        axis tight; box off; ylabel('Interp');
        set(gca, 'xtick', []);
    end
else
    nanIdx = [];
end

% ====================================================== %
% STEP 2: INTERPOLATE PEAK-DETECTED BLINKS
% ====================================================== %

assert(~any(isnan(dat.pupil)));
win             = hanning(11);
pupildatsmooth  = filter2(win.', dat.pupil, 'same');

dat.pupildiff   = diff(pupildatsmooth) - mean(diff(pupildatsmooth)) / std(diff(pupildatsmooth));
[peaks, loc]    = findpeaks(abs(dat.pupildiff), 'minpeakheight', 3*std(dat.pupildiff), 'minpeakdistance', 0.5*data.fsample);

if plotme, sp2 = subplot(412);
    plot(dat.time(2:end), dat.pupildiff);
    hold on; plot(dat.time(loc), peaks, '.');
    box off; ylabel('Peaks');
    set(gca, 'xtick', []); ylim([-500 500]);
    sp3 = subplot(413); hold on;
end

if ~isempty(peaks),
    
    if plotme,
        plot(dat.time, dat.pupil, 'color', [0.5 0.5 0.5]);
        hold on;
    end
    
    % convert peaks into blinksmp
    newblinksmp = nan(length(peaks), 2);
    for p = 1:length(peaks),
        newblinksmp(p, 1) = loc(p) - padding * data.fsample; % peak detected will be eye-opening again
        newblinksmp(p, 2) = loc(p) + (padding+0.2) * data.fsample;
    end
    
    % merge 2 blinks into 1 if they are < 250 ms together (coalesce)
    coalesce = 0.250;
    for b = 1:size(newblinksmp, 1)-1,
        if newblinksmp(b+1, 1) - newblinksmp(b, 2) < coalesce * data.fsample,
            newblinksmp(b, 2) = newblinksmp(b+1, 2);
            newblinksmp(b+1, :) = nan;
        end
    end
    % remove those duplicates
    newblinksmp(isnan(nanmean(newblinksmp, 2)), :) = [];
    
    % make sure none are outside of the data range
    newblinksmp(newblinksmp < 1) = 1;
    newblinksmp(newblinksmp > length(dat.pupil)) = length(dat.pupil) -1;
    
    % make the pupil NaN at those points
    for b = 1:size(newblinksmp,1),
        dat.pupil(newblinksmp(b,1):newblinksmp(b,2)) = NaN;
    end
    
    % save NaN idx
    nanIdx = [nanIdx find(isnan(dat.pupil))];
    
    % interpolate linearly
    dat.gazex(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.gazex(~isnan(dat.pupil)), find(isnan(dat.pupil)), ...
        'nearest', 'extrap');
    dat.gazey(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.gazey(~isnan(dat.pupil)), find(isnan(dat.pupil)), ...
        'nearest', 'extrap');
    
    % now the pupil itself
    dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
        dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'linear');
else
    newblinksmp = [];
end

% remove remaining nans (probably at the end)
dat.pupil(isnan(dat.pupil)) = interp1(find(~isnan(dat.pupil)), ...
    dat.pupil(~isnan(dat.pupil)), find(isnan(dat.pupil)), 'nearest', 'extrap');
newpupil = dat.pupil;

% link axes
if plotme,
    plot(dat.time, dat.pupil, 'b');
    ylim([min(dat.pupil)*0.9 max(dat.pupil)*1.1]);
    box off; ylabel('Clean');
    try
        linkaxes([sp1 sp2 sp3], 'x');
        set([sp1 sp2 sp3], 'tickdir', 'out');
    end
    xlim([dat.time(1)-10 dat.time(end)+10]);
end
% newblinksmp = [];
% newpupil = [];

% combine interpolated data matrix

%interpol=[padblinksmp;newblinksmp; nandat];
interpol=[padblinksmp;newblinksmp];
interpol=sortrows(interpol);
dat.interpol=interpol;

end
