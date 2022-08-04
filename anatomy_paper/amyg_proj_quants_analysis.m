%%% get amyg inputs to por quantiles for anatomy paper figure 4
% use circnormed projection images 
%
%%% updated 21-3-31 on thermaltake

%% setup parameters
filetable_name = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\module_input_ratios_filelist.xlsx';
% pathway_name = 'amyg_por';
%     pathway_name = 'lp_por';
    pathway_name = 'dlgn_por';
save_analysis_results = 0; 

findpatches_pars.threshMode = 'intensityQuantiles';
    findpatches_pars.nQuantiles = 6;
    findpatches_pars.patchQuantiles = [4 5  6]; % quantiles selected to be counted as patches; higher numbers = brighter pixels; use [5 6] to match sincich and horton
findpatches_pars.diskblurradius_um = 29;
findpatches_pars.minAreaPatch_squm = 0;   % min area in square microns a blob must contain to be considered a patch
findpatches_pars.maxAreaPatch_squm = inf; 
findpatches_pars.blur_using_only_roi = 1; % if true, do not use non-roi pixels for blurring (should==1 if image has zeros/nans outside roi, such as for blood vessels/image borders)  
findpatches_pars.include_interior_nonroi_in_roi = 1; % include non-roi pixels completely contained within patches as part of the surrounding patch
findpatches_pars.show_plots = 0;
findpatches_pars.save_patchdata = 0; 

permutation_pars.npermutations = 100; % number of permutations; ?10,000 recommended to asymptote
permutation_pars.resample_pix = 1;  %%% resample pixels in proj image; set pixel area to resampled_pix_area_squm
    permutation_pars.resampled_pix_area_squm = 100; % area of a pixel in square microns after resampling
permutation_pars.normMethod ='subtractBaseline_divideByRoiMean';   
permutation_pars.custom_interpatches = 0; % use something other than all non-patch pixels as interpatch pixels
    permutation_pars.interpatchQuantiles = [1 2 3]; % quantiles selected to be counted as interpatches, if custom_interpatches==true; lower numbers = darker pixels

%% analysis
filetable = readtable(filetable_name);
filetable = filetable(strcmp(filetable.pathway, pathway_name),:); %%% keep only amyg to por subjects
nsubs = height(filetable); 
filetable.patchInterpatchRatio_projfilt = NaN(nsubs,1); % clear old results
filetable.patchInterpatchRatio_projorig = NaN(nsubs,1); % clear old results
filetable.PIRatio_baselined_normed = NaN(nsubs,1); % clear old results
filetable.patch_analysis = cell(nsubs,1);
filetable.proj_intens_by_quant = NaN(nsubs,findpatches_pars.nQuantiles); %%% mean projimage intensity in each m2 quantile
filetable.proj_intens_by_quant_normed = NaN(nsubs,findpatches_pars.nQuantiles); %%% mean projimage intensity in each m2 quantile

for isub = 1:nsubs
    filetable.patch_analysis{isub} = analyze_patchiness(filetable.m2_filt{isub}, filetable.proj_filt{isub}, filetable.roi{isub},...
        filetable.baseline{isub}, filetable.zoom(isub), filetable.scope{isub}, findpatches_pars, permutation_pars); % do main patch and projection analysis
    filetable.patchInterpatchRatio_projfilt(isub) = filetable.patch_analysis{isub}.permutation_test.patchInterpatchRatio; % store this subject's results
    %%%% converting to double results in negative values after subtracting baseline 
    projimage_baselined = double(filetable.patch_analysis{isub}.permutation_test.comparison_image) - filetable.patch_analysis{isub}.permutation_test.baseline; 
    for iquant = 1:findpatches_pars.nQuantiles % get mean input to each quantile for this subject
        this_quant_image = filetable.patch_analysis{isub}.patchdata.quantile_table.quantimage{iquant};
        filetable.proj_intens_by_quant(isub,iquant) = mean(projimage_baselined(this_quant_image)); % get mean intensity in this quant for this subject
    end
    max_intens_this_sub = max(filetable.proj_intens_by_quant(isub,:));
    filetable.proj_intens_by_quant_normed(isub,:) = filetable.proj_intens_by_quant(isub,:) ./ max_intens_this_sub;
end

%% summary stats
mean_normed_intens_by_quant = mean(filetable.proj_intens_by_quant_normed);
sem_normed_intens_by_quant = std(filetable.proj_intens_by_quant_normed) ./ sqrt(nsubs);

if save_analysis_results
    save amyg_proj_quants_analysis filetable findpatches_pars mean_normed_intens_by_quant nsubs permutation_pars sem_normed_intens_by_quant
end