% full patch analysis
%    analyze_patchiness(m2file, comparison_file, roifile, baselineFile, zoom, scope, findpatches_pars, permutation_pars, perimeter_pars)
% lasted updated 2020/9/03 on thermaltake

function [allresults] = analyze_patchiness(m2file, comparison_file, roifile, baselineFile, zoom, scope, findpatches_pars, permutation_pars, perimeter_pars)

%% select files
% m2file = vardefault('m2file','15157_s1_Venus_x3_2_cropped1.tif');
% comparison_file = vardefault('comparison_file','15157_s1_tdT_x3_2_cropped1.tif';
roifile = vardefault('roifile',[]);
baselineFile = vardefault('baselineFile',[]);
% zoom = vardefault('zoom',3.2); 
scope = vardefault('scope','epimicro');



%% output and plotting - does not affect quantitative analysis
fig_pars = vardefault('fig_pars',struct);
fig_pars.show_plots = field_default(fig_pars,'show_plots',0); % show m2 and comparison images with patch and interpatch borders
fig_pars.patch_border_width = field_default(fig_pars,'patch_border_width',20); %%% width of patch borders 
fig_pars.patch_border_color= field_default(fig_pars,'patch_border_color',[.5, 1, 1]); % use either an RGB triplet or one of the following strings: r/red, g/green, b/blue, y/yellow, m/magenta, c/cyan, w/white, k/black, none
fig_pars.interpatch_border_width = field_default(fig_pars,'interpatch_border_width',1); %%% width of patch borders 
fig_pars.interpatch_border_color= field_default(fig_pars,'interpatch_border_color',[.1 .5 .8]); 
save_output = 0; % save the results in a filename with name specified in output_filename
    output_filename = 'temp.mat';



%% patch-finding params
% patches will be automatically found by blurring and threshold pixel intensity in the m2file image
findpatches_pars = vardefault('findpatches_pars',struct);
findpatches_pars.threshMode = field_default(findpatches_pars,'threshMode','intensityQuantiles');
    findpatches_pars.nQuantiles = field_default(findpatches_pars,'nQuantiles',6);
    findpatches_pars.patchQuantiles = field_default(findpatches_pars,'patchQuantiles',[4 5  6]); % quantiles selected to be counted as patches; higher numbers = brighter pixels; use [5 6] to match sincich and horton
% findpatches_pars.threshMode = field_default(findpatches_pars,'threshMode','fractionOfMax');
%     findpatches_pars.threshFraction = field_default(findpatches_pars,'threshFraction',0.878); %%% thresh needed to put a pixel within a patch; value = fraction of intensity distance from min val within roi to max val within roi
% findpatches_pars.threshMode = field_default(findpatches_pars,'threshMode','fractionOfMaxMinusMin');
%     findpatches_pars.threshFraction = field_default(findpatches_pars,'threshFraction',0.8); %%% thresh needed to put a pixel within a patch; value = fraction of intensity distance from min val within roi to max val within roi
findpatches_pars.diskblurradius_um = field_default(findpatches_pars,'diskblurradius_um',29);
findpatches_pars.minAreaPatch_squm = field_default(findpatches_pars,'minAreaPatch_squm',0);   % min area in square microns a blob must contain to be considered a patch
findpatches_pars.maxAreaPatch_squm = field_default(findpatches_pars,'maxAreaPatch_squm',inf); 
findpatches_pars.blur_using_only_roi = field_default(findpatches_pars,'blur_using_only_roi',1); % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
findpatches_pars.include_interior_nonroi_in_roi = field_default(findpatches_pars,'include_interior_nonroi_in_roi',1); % include non-roi pixels completely contained within patches as part of the surrounding patch
findpatches_pars.show_plots = field_default(findpatches_pars,'show_plots',0);
do_patch_density_analysis = 0; % get statistics on patch spacing




%% permutation test params
% permutation test looks for significant differences in patches vs interpatches
%   by shuffling pixels in the image and comparing these randomized intensity distributions
%   to the actual patch vs interpatch distribution in the original images
permutation_pars = vardefault('permutation_pars',struct);
permutation_pars.npermutations = field_default(permutation_pars,'npermutations',100); % number of permutations; ?10,000 recommended to asymptote
permutation_pars.resample_pix = field_default(permutation_pars,'resample_pix',1);  %%% resample pixels in proj image; set pixel area to resampled_pix_area_squm
    % at epi scope x3.2, 1 pixel area == 3.84 squm.... area of a 60um-radius patch == 3600 squm
    permutation_pars.resampled_pix_area_squm = field_default(permutation_pars,'resampled_pix_area_squm',100); % area of a pixel in square microns after resampling
permutation_pars.normMethod = field_default(permutation_pars,'normMethod','subtractBaseline_divideByRoiMean');   
%     permutation_pars.normMethod = field_default(permutation_pars,'normMethod','none');
%     permutation_pars.normMethod = field_default(permutation_pars,'normMethod','divideByRoiMean');
permutation_pars.custom_interpatches = field_default(permutation_pars,'custom_interpatches',0); % use something other than all non-patch pixels as interpatch pixels
    permutation_pars.interpatchQuantiles = field_default(permutation_pars,'interpatchQuantiles',[1 2 3]); % quantiles selected to be counted as interpatches, if custom_interpatches==true; lower numbers = darker pixels


%% perimeter test params
% perimeter test will look for significant differences in pixel intensity 
%   in patches vs. in a ring around each patch (considered to be the corresponding interpatch)
%   using a paired t-test
%%%%% NB: whatever area is considered to be interpatches in previous analyses will NOT match the interpatches used in the perimeter test
perimeter_pars = vardefault('perimeter_pars',struct);
perimeter_pars.intptchThicknessUM = field_default(perimeter_pars,'intptchThicknessUM',30); % width of the ring to draw around patches to treat as interpatches
perimeter_pars.pars.min_area_patch_sqmm = field_default(perimeter_pars,'min_area_patch_sqmm',0); % minimum patch area to analyze a given patch as part of the perimeter test
% perimeter_pars.normMethod = field_default(perimeter_pars,'normMethod','subtractBaseline'); 
% perimeter_pars.normMethod = field_default(perimeter_pars,'normMethod','subtractBaseline_divideByRoiMean'); % scale pixel intensity to average intensity within roi
% perimeter_pars.normMethod = field_default(perimeter_pars,'normMethod','subtractBaseline_divideByRoiMax'); % scale pixel intensity to max intensity within roi
perimeter_pars.normMethod = field_default(perimeter_pars,'normMethod','subtractBaseline_divideByPatchesIntptchs'); % scale intensity to avg of patches and interpatches, not avg of ROI




%% perform analyses
patchiness_pipeline();
    

