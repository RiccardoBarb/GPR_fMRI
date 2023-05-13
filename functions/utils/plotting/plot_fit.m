function plot_fit(yt,labels,mup,cv,vox,opt)
% The fuction plots the predicted continuous activity for the selected cv-fold
% in a certain amount of selected voxels. There are 2 plotting modes:
% 2d - very similar to plot_raw with an overlay of the predicted continuous 
% activity for the current voxels.
% 3d - a 3d view of the multivariate continuous activity pattern for the 
% selected voxels
% IN: yt - matrix of the sampled dependent variable of size (ntrials X 
%     nvoxels X ncv-folds) 
%     labels - matrix of the sampled independent variable of size (ntrials 
%     X ncv-folds
%     mup - continuous voxel-wise predicted activity. Output of the
%     function ME_cmGPR_mean.
%     cv - which cv fold you want to plot
%     vox - Integer vector representing the list of voxels to plot. If the 
%     list contains more than 10 elements we just plot the first 10 
%     opt - optional argument with plotting mode. Default = '2d'
% Author: Riccardo Barbieri

%set default plotting mode
if nargin < 6
    opt = '2d';
end
%set default if voxels to plot are missing
if nargin < 6 && size(yt,2)<=10
    sz =  size(yt,2);
    vox = 1:sz;
elseif nargin < 6 && size(yt,2)>10
    sz = 10;
    vox = 1:sz;
else
    if length(vox)>10
        fprintf('Plotting only the first 10 voxels in the list')
        vox =vox(1:10);
    end
    sz = length(vox);
end

%----------------------------------2d mode--------------------------------%
if lower(opt) == '2d'
    
    yl = max(max(max(yt)));
    
    figure('Name','Voxel activity vs predicted activity','NumberTitle','off');
    for j = 1:sz
        subplot(ceil(sz/2),2,j)
        hold on
        plot (labels(:),reshape(yt(:,vox(j),:),size(yt,1)*size(yt,3),1),'.b');
        hold on
        plot(squeeze(mup(:,vox(j),cv)),'r-')
        title(sprintf('voxel number %d',j))
        xlim ([0,360]);
    end
%----------------------------------3d mode--------------------------------%
elseif lower(opt) == '3d'
    
    for j=1:sz
        m_mup(:,j) = squeeze(mup(:,vox(j),1)-mean(mup(:,vox(j),cv)));
    end
    
    figure('Name','Continuous activity pattern in current slice','NumberTitle','off',...
        'rend','painters','pos',[500 500 2000 300],'visible','on');
    
    subplot(1,2,1)
    surf(m_mup(:,vox,cv),'MeshStyle','Column','EdgeColor','k','FaceAlpha',0.5),
    xlabel('Voxels in current slice')
    ylabel('Motion directions')
    zlabel('Mean centered predicted activity')
    view(40,40)
    subplot(1,2,2)
    imagesc(m_mup')
    ylabel('Voxels in current slice')
    xlabel('Motion directions')
    c = colorbar
    c.Label.String= 'Mean centered predicted activity'
end
