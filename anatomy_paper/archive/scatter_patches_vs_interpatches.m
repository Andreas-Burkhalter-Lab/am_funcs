%%% scatter_patches_vs_interpatches
% allpaths must be open
%%% updated 2019/69/23

tablerows_to_plot = 2:5;
single_marker_color = 0;
    marker_color = [0 0 1];
plot_log_intensities = 0; 
use_filtered_proj_image = 0; % use cirnormed proj image rather than raw proj image
    
cc = 0; % counter
for i = tablerows_to_plot
    cc = cc+1; 
    if use_filtered_proj_image
        projimage_to_use = allpaths.projfilt{i};
    elseif ~use_filtered_proj_image
        projimage_to_use = allpaths.proj_orig{i};
    end        
    p(cc) = analyze_patchiness(allpaths.m2filt{i},projimage_to_use,allpaths.roi{i},allpaths.baseline{i},allpaths.zoom(i),allpaths.scope{i});
end
figure
for i = 1:length(tablerows_to_plot)
    hold on
    if plot_log_intensities
        ss = scatter(log(p(i).perimeter_test.patchTable.patchIntens),log(p(i).perimeter_test.patchTable.interpatchIntens),'o','filled');
    elseif ~plot_log_intensities
        ss = scatter(p(i).perimeter_test.patchTable.patchIntens,p(i).perimeter_test.patchTable.interpatchIntens,'o','filled');        
    end
    if single_marker_color
        ss.MarkerFaceColor = marker_color;
        ss.MarkerEdgeColor = marker_color;
    end
    hold off
end

xymax = max([ylim xlim]); 
% xlim([0 xymax]); 
% ylim([0 xymax]);
xx = .001:.001:max([xlim ylim]);
hold on
plot(xx,xx,'k--')
legend(cellstr(num2str(cell2mat(allpaths.sub(tablerows_to_plot))))')
% legend(cellstr(num2str(allpaths.sub(tablerows_to_plot)))')
xlabel('Patch EGFP Optical Density')
ylabel('Interpatch EGFP Optical Density')
% % % % title(p(1).scalemethod)



%%% for use with old version of findpatches
%     p(c) = analyze_patchiness(allpaths.m2filt{i},logical(allpaths.patches{i}.labelmat_large),allpaths.projfilt{i},allpaths.roi{i},allpaths.baseline{i},allpaths.zoom{i},allpaths.scope{i},40); 
