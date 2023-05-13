function [data] = rdk_convertrawdata(scan_dir,nsubj,nruns,type,target_dir1)

%function [data] = rdk_getrawdata(scan_dir,nsubj,nruns,type)
%
%  The function get the full path of the fmri data for further processing
%
%
%  The data are assumed to be organized with the following structure:
%
%  -scan_folder
%    -subject_folder
%      -scans folder (4scout,1mprage,2field,1epi+2field+1epi... and so on
%        -scan_file.dcm
%
%   IN: scan_dir is the data directory containing separate folders for each
%       subject;
%       nsubj is the vector of subject numbers of the dataset
%       nruns is a vector containing the run number in ascendent order
%       type is the extension of the image we want to analyze (eg. 'img')
%       darget_dir is the target directory of the converted data
%       
%   OUT: data is a structure arranged in length(nsubj) X length(nruns)
%             .subject: subject name
%             .run: run number
%             .scans: full path of the scans
%
% Riccardo 18/01/2019: written

form = ['*.',type];

cd(scan_dir);
subj = dir('s*');
fprintf('\n loading....');

for subj_ = nsubj

    target_dir = target_dir1;%reset target_dir name

    cd(scan_dir);
    sname = subj(subj_).name;
    fprintf('\n converting data for subject %s which is the number %02i in this analysis:', sname, subj_);
    target_dir = strcat(target_dir,'/',subj(subj_).name(1:end-3));
    if exist(target_dir)==7
      fprintf('\n this folder already exist! \n assuming the data belong to the same subject they will be saved in day2 subfolder')
       target_dir=strcat(target_dir,'/',subj(subj_).name(1:end-3),'_day2');            
    else
        fprintf('\n saving data in %s',target_dir);
    end
    cd(sname);
    runs = dir('0*');
    
    for r_ = nruns
        
        rname = runs(r_).name;
        fprintf('\n RUN %i (%s)',r_,rname);
        data.subject{subj_} = sname;
        data.run{subj_,r_} = rname;
        cd(scan_dir);
        cd(sname);
        cd(rname);
        scans = dir(form);
        clear files
        for sc_=1:length(scans)
        data.scans{subj_,r_,sc_}=fullfile(scan_dir,sname,rname,scans(sc_).name);   
        end
        hdr = spm_dicom_headers(strvcat(data.scans{subj_,r_,:}), true);
        out = spm_dicom_convert(hdr,'all','series','img',target_dir);
    end

end
end