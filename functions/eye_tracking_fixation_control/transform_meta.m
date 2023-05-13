function [meta_fix] = transform_meta(meta)

% this function is written to exclude the training run (= the first run)
% from the meta data structure
% FMT 270320

for n=1:length(meta)
   
    meta_fix(n).subject=meta(n).subject;
    meta_fix(n).noise_thresh=meta(n).noise_thresh;
    
    if isstruct(meta(n).exclud_trial_id_per_run)
    meta_fix(n).exclud_trial_id_per_run=meta(n).exclud_trial_id_per_run(2:6);
    meta_fix(n).is_exclude=length([meta_fix(n).exclud_trial_id_per_run(:).trial]);
    
    else
        meta_fix(n).exclud_trial_id_per_run=[];
    end
    
end

end