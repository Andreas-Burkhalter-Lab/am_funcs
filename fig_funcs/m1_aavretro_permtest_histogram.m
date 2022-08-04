%%% for m1 aavretro in v1 case, make histogram of patch:interpatch intensity ratio permtest distribution
% updated 2020/9/23 on thermaltake

close all

% use lgn to v1 proj image as m2 image for patch finding
% both the s1 yfp and s2 tdt file were aligned to each other; don't use originals of either one
m2_file = 'F:\thesis\m1_retro_aav\19098_s1_lgn-yfp_x4_aligned-to-s2_normed62pix100um.tif';
% comparison_file = 'F:\thesis\m1_retro_aav\19098_s2_rAAVtdT_x4_aligned.tif';
    comparison_file = 'F:\thesis\m1_retro_aav\19098_s1_overlay_x4_aligned-to-s2_dendrites.png'; % use dendrite lengths as intensity image
roi_file = 'F:\thesis\m1_retro_aav\19098_s1_overlay_x4_aligned-to-s2_antlat_analysis_roi.png'; 
baseline_file = 'F:\thesis\m1_retro_aav\19098_s1_overlay_x4_aligned-to-s2_baseline.png'; 
zoom = 4;
scope = 'epimicro';
run_analysis = 0; %%% rerun analysis; alternatively, load results from previous analysis
save_image = 1; %%% save histogram figure as .svg
    image_savename = ['F:\thesis\figs\fig 5 -- ', strrep(getfname(baseline_file),'baseline','permtest_dend')];

% plotting options
xlimits = [0 1.6]; 
nbins = 50; 
bar_border_thickness = 1.7; 
edge_color = [0 0 0]; 
face_color = [0.3 0.3 0.3]; 
hist_norm_mode = 'Probability'; 
font_name = 'Arial';
axis_thickness = 2;
axis_font_bold = 1; % font bold or not bold
axis_font_size = 15;

% analysis options
findpatches_pars.threshMode = 'intensityQuantiles';
    findpatches_pars.nQuantiles = 6;
    findpatches_pars.patchQuantiles = [4 5  6]; % quantiles selected to be counted as patches; higher numbers = brighter pixels; use [5 6] to match sincich and horton
findpatches_pars.diskblurradius_um = 29;
findpatches_pars.minAreaPatch_squm = 0;   % min area in square microns a blob must contain to be considered a patch
findpatches_pars.maxAreaPatch_squm = inf; 
findpatches_pars.blur_using_only_roi = 1; % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
findpatches_pars.include_interior_nonroi_in_roi = 1; % include non-roi pixels completely contained within patches as part of the surrounding patch
findpatches_pars.show_plots = 0;
findpatches_pars.quantborders_ops.quants = [ 4]; % output patch/interpatch border

permutation_pars = vardefault('permutation_pars',struct);
permutation_pars.npermutations = 1e6; % number of permutations; D'Souza et al. used 1e5
permutation_pars.resample_pix = 1;  %%% resample pixels in proj image; set pixel area to resampled_pix_area_squm
    permutation_pars.resampled_pix_area_squm = 200; % area of a pixel in square microns after resampling
permutation_pars.normMethod = 'subtractBaseline_divideByRoiMean';   
permutation_pars.custom_interpatches = 0;

%% load or run analysis
% load('F:\thesis\m1_retro_aav\19098_s1_overlay_x4_aligned-to-s2_antlat_permtest.mat'); %%% load analysis results

if run_analysis
    patchiness_analysis = analyze_patchiness(m2_file, comparison_file, roi_file, baseline_file, zoom, scope, findpatches_pars, permutation_pars); 
end
nonshuffled_data_ratio = patchiness_analysis.permutation_test.patchInterpatchRatio; 

%% plotting
histfig = figure; 
hhist = histogram(patchiness_analysis.permutation_test.shuffledist_patchInterpatchRatio,nbins);
hhist.LineWidth = bar_border_thickness; 
hhist.EdgeColor = edge_color;
hhist.FaceColor = face_color; 
hhist.Normalization = hist_norm_mode; 
set(gca,'FontSize', axis_font_size)
set(gca,'FontName', font_name);
set(gca,'LineWidth',axis_thickness)
if axis_font_bold; set(gca,'FontWeight','bold'); end
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
set(gca,'Box','off')
xlim(xlimits)

% add a point indicating the non-shuffled patch:interpatch ratio
hold on
% hscatter = scatter(nonshuffled_data_ratio, 0.1*max(hhist.Values)); % add point near the x axis
hscatter = scatter(nonshuffled_data_ratio, 0, '.r', 'SizeData', 200); % add point on the x axis
hold off

% save image
if save_image %%% save histogram figure as .svg
    saveas(histfig, image_savename, 'svg') % don't use .eps; eps sometimes breaks up the image into pieces
end

