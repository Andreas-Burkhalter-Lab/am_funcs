% find patches from an M2 image and perform chi-square test of cell counts in different quantiles
%    [allresults] = analyze_patchiness_cells(m2file, cells_file, roifile, zoom, scope)
% lasted updated 2020/12/15 on thermaltake

function [res] = analyze_patchiness_cells(m2file, cells_file, roifile, zoom, scope)

%% select files
% m2file = '17185_Sec2b_Chrm2_xmicro3_2_normed51pix100um.tif';
%  cells_file = '17185_s4b_xmicro3_2_cells.png';
% roifile = '17185_s2b_xmicro3_2_v1area.png';    %%% set to [] to analyze all pixels in the image
% zoom = 3.2;
% scope = 'epimicro_old';
% scope = 'epimicro'; 



%% output and plotting - does not affect quantitative analysis
fig_pars.show_plots = 0; % show m2 and comparison images with patch and interpatch borders
fig_pars.quant_border_width = 1; %%% width of patch borders 
save_output = 0; % save the results in a filename with name specified in output_filename
    output_filename = 'temp.mat';



%% patch-finding params
% patches will be automatically found by blurring and threshold pixel intensity in the m2file image
findpatches_pars.threshMode = 'intensityQuantiles';
    findpatches_pars.nQuantiles = 6;
    findpatches_pars.topNQuantilesArePatches = 2; 
findpatches_pars.diskblurradius_um = 29;
findpatches_pars.minAreaPatch_squm = 0;   % min area in square microns a blob must contain to be considered a patch
findpatches_pars.maxAreaPatch_squm = inf; 
findpatches_pars.blur_using_only_roi = 1; % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
findpatches_pars.include_interior_nonroi_in_patches = 1; % include non-roi pixels completely contained within patches as part of the surrounding patch
findpatches_pars.show_plots = 0;
do_patch_density_analysis = 0; % get statistics on patch spacing


%% perform analyses
patchiness_pipeline_cells();
    

