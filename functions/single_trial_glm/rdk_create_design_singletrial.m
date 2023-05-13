function [matlabbatch]=rdk_create_design_singletrial(timing,id,scans,mov,subjects,runs,target_dir,correction)
%this script has been updated on august 2019 to run data analyses on the
%complete set of subjects of the CONDEC experiment.
%implemented the extraction of timings and condition informations based on
%the function "extract_design_info.m", it now include the trials to exclude
%based on the eyetracking data

clear matlabbatch
clear SPM
nruns = length(runs);
nsubj = length(subjects)
correction = correction*0.8;
if ~exist(target_dir)
    mkdir (target_dir);
    fprintf('creating new folder %s',target_dir);
    fprintf('/nsaving the current design matrix into %s',target_dir);
    
else
    
    fprintf('saving the current design matrix into %s',target_dir);
end


matlabbatch{1,1}.spm.stats.fmri_spec.dir = {target_dir};
matlabbatch{1,1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1,1}.spm.stats.fmri_spec.timing.RT = 0.8;
matlabbatch{1,1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1,1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

matlabbatch{1,1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1,1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1,1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1,1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1,1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1,1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1,1}.spm.stats.fmri_spec.cvi = 'AR(1)';
%% model specification

for r_=1:nruns
    tt_=0;
    coherences=unique(id(r_).stim.coherence);
    [s_coherences,cids]=sort(id(r_).stim.coherence);
    cscans={scans{1,runs(r_),1:end}}';
    idempty=cellfun(@isempty,cscans);
    matlabbatch{1}.spm.stats.fmri_spec.sess(r_).scans ={cscans{~idempty}}';
    
    for t_=1:numel(cids)
        %this is the name of each regressor: c (coherence level) (true
        %direction) (reported direction) (is to exclude?)
        matlabbatch{1}.spm.stats.fmri_spec.sess(r_).cond(1,t_).name=...
            ['c ' num2str(s_coherences(t_)) ' ' num2str(id(r_).stim.direction(cids(t_)))...
            ' ' num2str(id(r_).resp.direction(cids(t_))) ' ' num2str(id(r_).resp.nofixation(cids(t_)))];
        matlabbatch{1}.spm.stats.fmri_spec.sess(r_).cond(1,t_).onset=timing(r_).stim.onsets(cids(t_))-correction;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r_).cond(1,t_).duration=2;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r_).cond(1,t_).tmod=0;
        matlabbatch{1}.spm.stats.fmri_spec.sess(r_).cond(1,t_).pmod=struct([]);
        matlabbatch{1}.spm.stats.fmri_spec.sess(r_).cond(1,t_).orth=1;
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess(r_).multi_reg ={mov{runs(r_)}}';
end

cd(target_dir);
%save the matlabbatch file
save('trialwise_design.mat','matlabbatch');
%specify the model
spm_jobman('run', matlabbatch);
%     %% estimate 1st level model
load SPM;
spm_spm(SPM);
end
