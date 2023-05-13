function plot_raw(yt,labels,vox)
% Simple utility fuction for plotting the relationship between the sampled
% dependent and independent variable.
% Single voxels are plotted in different subplot. If there are more than 10
% voxels, we only plot 10.
% IN: yt - matrix of the sampled dependent variable of size (ntrials X 
%     nvoxels X ncv-folds) 
%     labels - matrix of the sampled independent variable of size (ntrials 
%     X ncv-folds
%     vox - optional argument. Integer vector representing the list of 
%     voxels to plot: e.g. plot_raw(yt,labels,[1,2,3]). If the list 
%     contains more than 10 elements we just plot the first 10 
% Author: Riccardo Barbieri

yl = max(max(max(yt)));

%set default if voxels to plot are missing
if nargin < 3 && size(yt,2)<=10
    sz =  size(yt,2);
    vox = 1:sz;
elseif nargin < 3 && size(yt,2)>10
    sz = 10;
    vox = 1:sz;
else
    if length(vox)>10
        fprintf('Plotting only the first 10 voxels in the list')
        vox =vox(1:10);
    end
    sz = length(vox);
end


figure()

for j = 1:sz
    subplot(ceil(sz/2),2,j)
    plot (labels(:),reshape(yt(:,vox(j),:),size(yt,1)*size(yt,3),1),'.b');
    hold on
    line ([0,360],[0,0],'color','r','linestyle','--');
    title(sprintf('voxel number %d',vox(j)))
    ylim([-yl-(2/3).*yl,yl+(2/3).*yl])
    xlim([0,360])
end

end