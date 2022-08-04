%%% get comparison-image vs. m2 quantiles for physiology paper figure
% use circnormed projection images 
%
%%% updated 2020-11-26 on thermaltake


% clear
% close all

%% analysis parameters
filetable_name = 'G:\rinaldo_sample_data\patchiness_filelist_sample.xlsx'; % excel file listing files and image parameters
pathway_name = 'som'; %% specify which pathway from the filetable will be analyzed
save_shuffled_proj_data = 0; 
    resampled_pix_area_squm = 5; % area of a pixel in square microns after resampling, to use in shuffled data
    shuffled_data_filename = 'shuffled_patchiness_data'; 

findpatches_pars.threshMode = 'intensityQuantiles';
    findpatches_pars.nQuantiles = 6;
    findpatches_pars.patchQuantiles = [4 5  6]; % quantiles selected to be counted as patches; higher numbers = brighter pixels; use [5 6] to match sincich and horton
findpatches_pars.diskblurradius_um = 29;
findpatches_pars.minAreaPatch_squm = 0;   % min area in square microns a blob must contain to be considered a patch
findpatches_pars.maxAreaPatch_squm = inf; 
findpatches_pars.blur_using_only_roi = 1; % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
findpatches_pars.include_interior_nonroi_in_roi = 1; % include non-roi pixels completely contained within patches as part of the surrounding patch
findpatches_pars.show_plots = 0;
findpatches_pars.save_patchdata = 0; %%% save the patchdata struct as a separate file and run save_quant_borders_image.m on all cases
    findpatches_pars.quantborders_ops.quants = [4 ];    % draw the borders between these quants and the quants directly below them

permutation_pars.npermutations = 10; % number of permutations; ?10,000 recommended to asymptote
permutation_pars.resample_pix = 1;  %%% resample pixels in proj image; set pixel area to resampled_pix_area_squm
    permutation_pars.resampled_pix_area_squm = 100; % area of a pixel in square microns after resampling
permutation_pars.normMethod ='subtractBaseline_divideByRoiMean';   
permutation_pars.custom_interpatches = 0; % use something other than all non-patch pixels as interpatch pixels
    permutation_pars.interpatchQuantiles = [1 2 3]; % quantiles selected to be counted as interpatches, if custom_interpatches==true; lower numbers = darker pixels

    %% plotting parameters
run_plotting = 1; %%% after analysis,  run physpaper_anatomy_quants_plots.m
    
plot_group_bars = 1; %% plot group means and SEM
plot_individual_bars = 1; %% plot each individual subject in a subplot
    subplot_rowcol = [3 2]; 
plot_shuffled_data = 0; %%% plot shuffled instead of original data
    
bar_colors = repmat(linspace(0.08,0.86,6)',1,3);
bar_border_width = 1; 
bar_border_color = 'none';
ebar_color = 'k'; 
ebar_linewidth = 2.5; 
axes_font_size = 22;
axes_line_width = 3; 
xlimits = [0.2, 6.8]; 
x_tick_length = [0 0];
fig_width = 520;
fig_height = 500; 
show_title = 0; 
save_quants_bar_fig = 0; %%% save output as .eps image
    savename_quants_bar_fig = ['F:\thesis\figs\fig 4 -- quants_bargraph_fig_', pathway_name];

    
    %% run analysis and plotting
quants_plotting()


