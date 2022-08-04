%%%% output images of selected ROIs with object number
% argin is output from get_tuning_curves.m

function make_analyzed_rois_images(results)


fhandle = figure('units','normalized','outerposition',[0 0 1 1]);

% different-colored masks of each roi
allcells = zeros([size(results.meanImage),3]);
mm = results.meanImage;
mm = mm/(max(max(mm)));
zz = false([size(results.meanImage),3]);
for i = 1:height(results.sftf_tuning)
    randcolor = rand(3,1);
    inds1 = zz; inds1(:,:,1) = results.sftf_tuning.cellimage{i};
    inds2 = zz; inds2(:,:,2) = results.sftf_tuning.cellimage{i};
    inds3 = zz; inds3(:,:,3) = results.sftf_tuning.cellimage{i};
%     allcells(inds1) = randcolor(1)*mm(results.sftf_tuning.cellimage{i}); % scale to intensity from meanImage   
%     allcells(inds2) = randcolor(2)*mm(results.sftf_tuning.cellimage{i});
%     allcells(inds3) = randcolor(3)*mm(results.sftf_tuning.cellimage{i});
    allcells(inds1) = randcolor(1); % don't scale to intensity from meanImage
    allcells(inds2) = randcolor(2);
    allcells(inds3) = randcolor(3);
end
image(allcells)
xtickvals = get(gca,'Xtick');
ytickvals = get(gca,'Ytick');
xlimvals = get(gca,'XLim');
ylimvals = get(gca,'YLim');
saveas(fhandle,'roi_masks','tif')

% map of numbered indices of each roi center corresponding to rows in
% results.sftf_tuning
scatter(results.sftf_tuning.centeryx(:,2),results.sftf_tuning.centeryx(:,1),'.')
text(results.sftf_tuning.centeryx(:,2),results.sftf_tuning.centeryx(:,1),cellstr(num2str([1:height(results.sftf_tuning)]')))
set(gca,'Ydir','reverse')
set(gca,'Xtick',xtickvals);
set(gca,'Ytick',ytickvals);
set(gca,'XLim',xlimvals);
set(gca,'YLim',ylimvals)
saveas(fhandle,'roi_inds_map','tif')

% time-averaged image of F values
imagesc(results.meanImage);
colormap(gray);
set(gca,'Xtick',xtickvals);
set(gca,'Ytick',ytickvals);
set(gca,'XLim',xlimvals);
set(gca,'YLim',ylimvals)
saveas(fhandle,'mean_image','tif')

close(fhandle)