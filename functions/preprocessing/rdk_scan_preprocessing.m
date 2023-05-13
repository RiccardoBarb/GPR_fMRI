function [data]=rdk_scan_preprocessing(data,type)

%% realign
if strcmp(type,'realign_reslice')
    
    for s_=1
        for r_=1:size(data.run,2)
            
            
            
            clear P;
            
            clear files;
            
            for i=1:size(data.scans,3)
                files{i}=data.scans{s_,r_,i};
            end
            
            P  = strvcat(files{~isempty(files),:});
            
            
            flags.quality = 0.9;
            flags.fwhm    = 5;
            flags.sep     = 4;
            flags.rtm     = 0;
            flags.PW      = '';
            flags.interp  = 2;
            flags.wrap    = [0 0 0];
            fprintf('\n REALIGNING RUN %s SUBJECT %s ',data.run{r_},data.subject{s_});
            spm_realign(P,flags);
            P            = char(P);
            flags.mask   = 1;
            flags.interp = 4;
            flags.which  = [2 1];
            flags.wrap   = [0 0 0];
            flags.prefix = 'r';
            fprintf('\n RESLICING RUN %s SUBJECT %s ',data.run{r_},data.subject{s_});
            spm_reslice(P,flags);
        end;
        
    end
end
if strcmp(type,'realign_est')
    
    for s_=1
        for r_=1:size(data.run,2)
            
            
            
            clear P;
            
            clear files;
            
            for i=1:size(data.scans,3)
                files{i}=data.scans{s_,r_,i};
            end
            
            P  = strvcat(files{~isempty(files),:});
            
            
            flags.quality = 0.9;
            flags.fwhm    = 5;
            flags.sep     = 4;
            flags.rtm     = 0;
            flags.PW      = '';
            flags.interp  = 2;
            flags.wrap    = [0 0 0];
            fprintf('\n REALIGNING RUN %s SUBJECT %s ',data.run{r_},data.subject{s_});
            spm_realign(P,flags);
        end
        
    end
end

if strcmp(type,'realign_custom')
    
    %this setup allow for realignment of day 1 scans with reslicing.
    %We only estimate the realignment for day2 without reslicing, as this
    %will happen in the next iteration of the code.
    
    res=1;
    for s_=1
        for r_=1:size(data.run,2)
            if r_>=8%again this number depends on the amount of scans that belong to day1
                res=0;
            end
            clear P;
            
            clear files;
            
            for i=1:size(data.scans,3)
                files{i}=data.scans{s_,r_,i};
            end
            
            P  = strvcat(files{~isempty(files),:});
            
            flags.quality = 0.9;
            flags.fwhm    = 5;
            flags.sep     = 4;
            flags.rtm     = 0;
            flags.PW      = '';
            flags.interp  = 2;
            flags.wrap    = [0 0 0];
            fprintf('\n REALIGNING RUN %s SUBJECT %s ',data.run{r_},data.subject{s_});
            spm_realign(P,flags);
            if res
                P            = char(P);
                flags.mask   = 1;
                flags.interp = 4;
                flags.which  = [2 1];
                flags.wrap   = [0 0 0];
                flags.prefix = 'r';
                fprintf('\n RESLICING RUN %s SUBJECT %s ',data.run{r_},data.subject{s_});
                spm_reslice(P,flags);
            end
        end
        
    end
end
%% coregistration
if strcmp(type,'coreg_custom')
    clear P
    %here we coregister the realigned and resliced images of day 1 with the
    %t1 image of day1. Then we coregister the realigned images of day 2
    %(only estimate) with the realigned and resliced images of day 1 and
    %reslice them.
    for s_=1
        clear files
        fprintf('\n Coregistering Day 1 SUBJECT %s ',data.subject{s_});
        %day 1
        if length(data.run)<8
            matlabbatch{1}.spm.spatial.coreg.estimate.ref =data.structural;
            matlabbatch{1}.spm.spatial.coreg.estimate.source = data.source;
            notempty = ~cellfun('isempty',data.scans);
            files=data.scans(notempty);
            matlabbatch{1}.spm.spatial.coreg.estimate.other = files;
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
           %save('coreg.mat','matlabbatch')
            spm_jobman('run',matlabbatch);
        else
            % day 2
            clear matlabbatch
            clear files
            fprintf('\n Coregistering Day 2 SUBJECT %s ',data.subject{s_});
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = data.structural;
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = data.source;
            notempty = ~cellfun('isempty',data.scans);
            files=data.scans(notempty);
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = files;
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
            %save('coreg2.mat','matlabbatch')
            spm_jobman('run',matlabbatch);
        end
    end
    
    
end



