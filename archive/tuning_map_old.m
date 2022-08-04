%%%% overlay heatmap of a variable on the time-averaged image of cells


param_to_map = 'sf_pref';
% param_to_map = 'sf_anovap';
logval = 1; % take log of specified parameter value
colormap_for_vals = parula;
patchoutlines_color = 'w';

% imagesc(t.reg_struct.movingimage_warped);
% hold on

imagesc(t.patchdata.im)
% imagesc(t.patchdata.imroiblurred) % blurred patches
cc = repmat(linspace(0,1,64)',1,3); % bw
cc = cc.^2;
set(gca,'Colormap',cc)
[ii jj] = find(t.patchdata.bnds_large_image);
hold on
scatter(jj,ii,'.',patchoutlines_color)



allvals = s{:,param_to_map};
if logval
    allvals = log(allvals);
end
minval = min(allvals);
maxval = max(allvals);
normedvals = [allvals - minval] ./ [maxval-minval];

nrois = height(s);
for iroi = 1:nrois
    clear ii jj
    [ii jj] = find(s.roi_image_reg{iroi});
    val_in_cmap = colormap_for_vals( round(normedvals(iroi)*63)+1,: );
    scatter(jj,ii,[],repmat(val_in_cmap,length(ii),1),'.')
end
% set(gca,'XTick',[])
% set(gca,'YTick',[])
hold off
    