function [data] = rdk_getscans(scan_dir,nsubj,nruns,totald1,tag,type,folder,correction)

%function [data] = rdk_getscans(scan_dir,nsubj,nruns,tag,type,folder,correction)
%
%  The function get the full path of the fmri data for further processing
%
%
%  The data are assumed to be organized with the following structure:
%
%  -scan_folder
%    -subject_folder
%      -run_folder
%        -scan_file.hdr
%        -scan_file.img
%
%   IN: scan_dir is the data directory containing separate folders for each
%       subject;
%       nsubj is the vector of subject numbers of the dataset
%       nruns is a vector containing the run number in ascendent order
%       tag is the character corresponding to the first letter of the image
%           file (eg. 'r','a','f','s')
%       type is the extension of the image we want to analyze (eg. 'img')
%       correction is an integer representing the number of TR we
%       discard to allow for magnetic saturation effects (only make sense
%       if we get functional scans)
%
%   OUT: data is a structure arranged in length(nsubj) X length(nruns)
%             .subject: subject name
%             .run: run number
%             .scans: full path of the scans
%
% Riccardo 16/01/2019: written

form=[tag,'*.',type];


subj=[dir(fullfile(scan_dir,'/s*'))];
fprintf('\n loading....');

for subj_=1:length(nsubj)
    sname=subj(nsubj(subj_)).name;
    fprintf('\n subject %s which is the number %02i in the dataset:', sname, nsubj(subj_));
    if strcmp(folder,'t1')
        runs=[dir(fullfile(scan_dir,sname,'/','T1w*'))];
    else
        runs=[dir(fullfile(scan_dir,sname,'/','rfMRI_*'))];
    end
    for r_=nruns
        if r_>=8 %this number depend on the amount of runs in the day1 folder
            d2=[dir(fullfile(scan_dir,sname,'/','s*day2'))];
            runs=[dir(fullfile(scan_dir,sname,'/','s*day2','/','rfMRI_*'))];
            r2_=r_-totald1;% amount of runs for day1
            rname=strcat(d2.name,'/',runs(r2_).name);
            fprintf('\n RUN %i (%s)',r2_,rname);
            data.subject{nsubj(subj_)}=sname;
            data.run{subj_,r_}=rname;
            scans=[dir(fullfile(scan_dir,sname,rname,form))];
            for sc_=1:length(scans)-correction
                data.scans{subj_,r_,sc_}=fullfile(scan_dir,sname,rname,scans(sc_+correction).name);
            end
        else
            rname=runs(r_).name;
            fprintf('\n RUN %i (%s)',r_,rname);
            data.subject{nsubj(subj_)}=sname;
            data.run{subj_,r_}=rname;
            scans=[dir(fullfile(scan_dir,sname,rname,form))];
            for sc_=1:length(scans)-correction
                data.scans{subj_,r_,sc_}=fullfile(scan_dir,sname,rname,scans(sc_+correction).name);
            end
        end
    end
end
end

