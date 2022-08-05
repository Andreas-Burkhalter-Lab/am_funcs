%%% create plot showing the areal major/minor axis lengths of patches and...
%%%     ... create plot showing patch densities across area
%
%%% first run patchcounting_analysis.m, using the top 1/5 of pixels as patches
%
%%% updated 2022/8/4


%% set parameters
patchcounting_filelist = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\patchcounting\patchcounting_filelist.xlsx';

%%%% specify which analyses to run and save
plot_patch_shapes = 0;
    save_shapeplot_image = 0; 
        savename_shapeplot = 'axis_ratios_fig'; 
    log_axes_ratio = 0; %% makes areal differences (v1+por vs. lm+li) much more visible
    shapes_plot_as_box = 1; % boxplot rather than bargraph
    scalefactor1 = 0.9;
    patch_shapes_width = 450*scalefactor1; 
    patch_shapes_height = 350*scalefactor1;
plot_patches_per_sqmm = 1; %%% patch areal density
    save_patches_per_sqmm = 0; 
    patches_per_sqmm_width = 300;
    patches_per_sqmm_height = 340;
        savename_patches_per_sqmm = 'patches_per_sqmm_fig';
    area_density_as_box = 0; % boxplot rather than bargraph
plot_patches_per_sqdeg = 0; %%% patch retinotopic density
    save_patches_per_sqdeg = 0; 
    patches_per_sqdeg_width = 300;
    patches_per_sqdeg_height = 340;
        savename_patches_per_sqdeg = 'patches_per_sqdeg_fig';
    retin_density_as_box = 0; % boxplot rather than bargraph
plot_area_ratio_axes = 0; % plot axes with no data to put area ratio data on in illustrator    
    savename_area_ratio_axes = 'area_ratios_blank_fig'; 
    scalefactor4 = 0.95;
    ylim_area_ratios = [-6.5, 6.5];
    area_ratio_width= 350*scalefactor4 ;
    area_ratio_height = 350*scalefactor4;
 output_stats = 1; % output results of tukey's test comparing areas
    
%%% plotting parameters
show_notch = 0;
boxplot_linewidth = 2; 
bar_face_color = [0.3 0.3 0.3]; % bar graph -  bar color
bar_line_width = 2; % bar graph - bar outline line width
error_line_color = [0 0 0]; % error bars color
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
plotfont = 'Arial'; % default = Helvetica
axis_font_size = 13; 
% ylimits = [0.85, 1.89];
% xlimits = [0.3, 4.8];
resolution_dpi = 300;

% findpatches parameters
findpatches_pars.minAreaPatch_squm = 600;
findpatches_pars.maxAreaPatch_squm = inf;
findpatches_pars.blur_using_only_roi = 1; % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
findpatches_pars.diskblurradius_um = 29;
% findpatches_pars.threshMode = 'intensityQuantiles';
%     findpatches_pars.nQuantiles = 6; % number of intensity levels to divide image into (if using 'intensityQuantiles')
%     findpatches_pars.patchQuantiles = [5 6]; % % quantiles selected to be counted as patches(if using 'intensityQuantiles'); higher numbers = brighter pixels
findpatches_pars.threshMode = 'fractionOfPopulation';
    findpatches_pars.fracPixSuperthresh = 1/5; % fraction of pixels marked as above threshold (if using 'fractionOfPopulation')
findpatches_pars.include_interior_nonroi_in_roi = 1; 
findpatches_pars.show_plots = 0;
findpatches_pars.visualize_shrinking = 0;

%% patch analysis
% load('C:\Users\Burkhalter Lab\Documents\anatomy_paper\patchcounting\patchcounting_results.mat'); % unless analysis needs to be re-run, load this dataset
filelist = readtable(patchcounting_filelist); 
excel_rows = 2:1+height(filelist); % use all cases
[cases, patchtable] = patchcounting_analysis(patchcounting_filelist, excel_rows, findpatches_pars); % comment out if data already loaded

%% plot patch shape data
if show_notch
    notch_toggle = 'on';
else
    notch_toggle = 'off';
end
if plot_patch_shapes
    if log_axes_ratio
        ratio_data = log(cases.mean_patch_major_minor_ratio);
    else
        ratio_data = cases.mean_patch_major_minor_ratio;
    end
    % output values of area axes widths to be plotted as crosses in illustrator
    stats_area_axes = grpstats(cases, {'area'}, {'mean','sem'}, 'DataVars','area_major_minor_ratio');
    stats_patch_axes = grpstats(cases, {'area'}, {'mean','sem'}, 'DataVars','mean_patch_major_minor_ratio');
    fig_shapes = figure;
    if shapes_plot_as_box
        bp = boxplot(ratio_data, cases.area, 'notch',notch_toggle, 'plotstyle','traditional','Colors','k');
    else % plot as bargraph
        bg = bar(1:4,stats_patch_axes.mean_mean_patch_major_minor_ratio);
        hold on
        eb = errorbar(1:4, stats_patch_axes.mean_mean_patch_major_minor_ratio, stats_patch_axes.sem_mean_patch_major_minor_ratio,'LineStyle','none');
        hold off
    end
    plot_formatting();
    set(fig_shapes,'Renderer', 'painters', 'Position', [50, 50, patch_shapes_width, patch_shapes_height]) % length x axis
    if save_shapeplot_image 
        print(fig_shapes,savename_shapeplot, '-dtiffn','-r300') %%% save image as file
    end
    
end


%% plot patch areal density (patches per square mm)
if plot_patches_per_sqmm
    cases = movevars(cases, 'patches_per_sqmm', 'After', 'sub');
     stats_areal_density = grpstats(cases, {'area'}, {'mean','sem'}, 'DataVars','patches_per_sqmm');
    fig_area_density = figure;
    if area_density_as_box 
        bp = boxplot(cases.patches_per_sqmm, cases.area,  'plotstyle','traditional', 'notch',notch_toggle, 'Colors','k');
    else % plot as bargraph
        bg = bar(1:4,stats_areal_density.mean_patches_per_sqmm);
        hold on
        eb = errorbar(1:4, stats_areal_density.mean_patches_per_sqmm, stats_areal_density.sem_patches_per_sqmm,'LineStyle','none');
        hold off
    end
    plot_formatting();
    set(fig_area_density,'Renderer', 'painters', 'Position', [50, 50, patches_per_sqmm_width, patches_per_sqmm_height]) % length x axis
    if save_patches_per_sqmm
        print(fig_area_density,savename_patches_per_sqmm, '-dtiffn','-r300') %%% save image as file
    end
end

%% plot patch retinotopic density (patches per square degree in visual field)
if plot_patches_per_sqdeg
    cases = movevars(cases, 'patches_per_sqdeg', 'After', 'sub');
    stats_retin_density = grpstats(cases, {'area'}, {'mean','sem'}, 'DataVars','patches_per_sqdeg');
    fig_ret_density = figure;
    if retin_density_as_box
        bp = boxplot(cases.patches_per_sqdeg, cases.area, 'plotstyle','traditional', 'notch',notch_toggle, 'Colors','k');
    else % plot as bargraph
        bg = bar(1:4,stats_retin_density.mean_patches_per_sqdeg);
        hold on
        eb = errorbar(1:4, stats_retin_density.mean_patches_per_sqdeg, stats_retin_density.sem_patches_per_sqdeg,'LineStyle','none');
        hold off
    end
    plot_formatting();
    set(fig_ret_density,'Renderer', 'painters', 'Position', [50, 50, patches_per_sqdeg_width, patches_per_sqdeg_height]) % length x axis
    if save_patches_per_sqdeg
        print(fig_ret_density,savename_patches_per_sqdeg, '-dtiffn','-r300') %%% save image as file
    end
end

if plot_area_ratio_axes
    fig_area_area_ratio_axes = figure;
    blank_data = linspace(-3,3,10);
    hplot = plot(blank_data,'Color',[1 1 1]); % line white so it's invisible
    plot_formatting();
    set(gca,'YTick',[0:6])
    ylim(ylim_area_ratios)
    set(fig_area_area_ratio_axes,'Renderer', 'painters', 'Position', [50, 50, area_ratio_width, area_ratio_height]) % length x axis
    print(fig_area_area_ratio_axes,savename_area_ratio_axes, '-dtiffn','-r300') %%% save image as file
end

if output_stats
    fprintf('patches_per_sqmm \n')
    [p,tab,stats] = anova1(cases.patches_per_sqmm,cases.area,'off')
    [res, means] = multcompare(stats,'CType','hsd')
    fprintf('mean_patch_major_minor_ratio \n')
    [p,tab,stats] = anova1(cases.mean_patch_major_minor_ratio,cases.area,'off')
    [res, means] = multcompare(stats,'CType','hsd')
    fprintf('patches_per_sqdeg \n')
    [p,tab,stats] = anova1(cases.patches_per_sqdeg,cases.area,'off')
    [res, means] = multcompare(stats,'CType','hsd')
end

clear('xlimits','ylimits')