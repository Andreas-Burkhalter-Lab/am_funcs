%%%% heatmap of pakan locomotion index over quantiles to show greater locomotion modulation in interpatches
% from subject 18164, session 2018-12-21
%
% updated 20/1/24

% load('F:\analyses\cor5-21.mat','filetable','roitable')

% close all

clear ops heatmap_pars

ops.meets_criteria = true(height(roitable),1); % reset meets_criteria
ops.meets_criteria = ops.meets_criteria & roitable.locm_cor_p > 0.05; 
ops.meets_criteria = ops.meets_criteria & roitable{:,'rf_sgnf'}; % if roi has detectable receptive field
ops.meets_criteria = ops.meets_criteria & strcmp(roitable.day, filetable.day{1});   %%%% only analyze particular recording sessions
ops.varname = 'pakan_loc_index';
ops.meets_criteria = ops.meets_criteria & ~isnan(roitable{:,ops.varname}); 

heatmap_pars.background_image_source = 'quant_levels_img';
heatmap_pars.colormap_for_paramvals = jet;
heatmap_pars.logval = 0;
heatmap_pars.save_image = 0; %%%% save after adjusting fig in this script, not in called function
heatmap_pars.savename = 'heatmap_over_image_fig';
heatmap_pars.resolution_dpi = 500; 

figvars = heatmap_over_image(roitable,ops.meets_criteria,'pakan_loc_index',filetable,heatmap_pars);

set(figvars.ax1,'XTick',[])
set(figvars.ax1,'YTick',[])
set(figvars.ax2,'XTick',[])
set(figvars.ax2,'YTick',[])
% figvars.cbar2.Position = [0.6    0.1100    0.0350    0.5];

%%% save figure
print(gcf, [heatmap_pars.savename,'.tif'], '-dtiffn',['-r',num2str(heatmap_pars.resolution_dpi)])
