function [extracted_labels,vectors] = extract_glm_design(scan_dir,condition,nconditions,betanumbers,gid,nruns)
n_multi_reg = 6; %6 head motion params as regressors

load(fullfile(scan_dir,'trialwise_design.mat'));

betas_per_run = length(matlabbatch{1,1}.spm.stats.fmri_spec.sess(1).cond)+n_multi_reg;

nvox =  numel(gid);

vectors = zeros(nconditions, nvox, nruns);

extracted_labels = zeros(nconditions,nruns);

for cond=1:nconditions
    
    fprintf('\b\b\b\b\b\b %0.1f%%',(cond/nconditions).*100); % \b is backspace
    for r=1:nruns
        betanumber = betanumbers(cond); % specify the beta number for a decoding condition
        betanumber = betanumber+(r-1)*betas_per_run;
        beta_image = sprintf('%s/beta_%04d.nii',scan_dir,betanumber);
        vol_hdr = spm_vol(beta_image);
        vol_vol = spm_read_vols(vol_hdr);
        vectors(cond,:,r)=vol_vol(gid);
        a = strsplit(matlabbatch{1,1}.spm.stats.fmri_spec.sess(r).cond(betanumbers(cond)).name);
        if size(a,2)>4 
            if a{5}=='0'
                extracted_labels(cond,r)=str2num(a{condition});
            elseif a{5}=='1'
                extracted_labels(cond,r)=NaN;% if the trial needs to be excluded
            end
        else
            extracted_labels(cond,r)=str2num(a{condition});
        end
    end
end
