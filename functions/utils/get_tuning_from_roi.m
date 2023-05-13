function [mean_roi_tuning,roi_precision,roi_tuning,beta_estimate,extracted_labels] = get_tuning_from_roi(tuning_file,roi_file,source_mask_file,mapping_file,...
    do_predictions,do_betaestimate,betanumbers,condition,argmax)

if isempty(do_predictions)
    do_predictions = 0;
end
fprintf('\nreading mask file\n')
%first we get the subject mask from the 1st level glm
mask = spm_read_vols(spm_vol(source_mask_file));
%then we get the roi mask and combine the two
roi = spm_read_vols(spm_vol(roi_file));
ids = mask & roi;
%here we obtain the linear index of the individual mask (same as doing
%find(mask))
load(mapping_file,'gid','gis');
%finally we reference the combination of the roi and the mask onto the
%linear index. Now the voxels of the roi mask can be referenced to the
%model (mup) voxel space.
roi_voxels = ids(gid);
%the following part is only to get the voxels in specific SL, for only
%center of SL this should not be done
v_ids = find(roi_voxels);
v_ids = v_ids(argmax);%comment out if not necessary to find best sl
v = zeros(size(v_ids,2),260);
for i_=1:size(v,1)
    v(i_,:)=[gis{v_ids(i_)}',repmat(0,1,size(v,2)-numel(gis{v_ids(i_)}))];
end
uniq_vox = unique(v);
roi_voxels = uniq_vox(2:end);
%until here

fprintf('loading voxel tuning of current subject...\n');
load(tuning_file,'mup');
fprintf('done!')
roi_tuning = zeros(size(mup,1),length(find(roi_voxels)),size(mup,3));
roi_tuning = mup(:,roi_voxels,:);
clear mup 
%commented part takes thee mean tuning instead of each cv folds
mean_roi_tuning = roi_tuning;%mean_roi_tuning = mean(roi_tuning,3);   

if do_predictions
    fprintf('\ncalculating balanced precision in current roi...')
    part = strsplit(tuning_file,'/');
    part{end} = 'predictions.mat';
    unionz = strcat(part,'/');
    prediction_file = strcat(unionz{1:end});
    load(prediction_file(1:end-1));
    predicted_roi = predicted_directions(roi_voxels,:);
    clear predicted_directions
    roi_precision = zeros(size(predicted_roi,1),1);
    for v_=1:size(predicted_roi,1)
        p = vertcat(predicted_roi{v_,:});
        d = vertcat(true_directions{:});
        roi_precision(v_) = bal_norm_circ_resp_dev(p,d,'trapz')-0.5;
    end
        fprintf('done!')
else
    roi_precision = [];
end
if do_betaestimate
    %--------------------% Get Betas and lables %-----------------------------%
scan_dir = source_mask_file(1:end-9); 
    fprintf('Extracting subject betas');
ntrials = numel(betanumbers);
[extracted_labels,vectors] = ...
    extract_glm_design(scan_dir,condition,ntrials,betanumbers,gid,10);
beta_estimate = vectors(:,roi_voxels,:);
fprintf('\nDone!');
    
else
    beta_estimate =[];
end