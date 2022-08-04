%%% get comparison-image vs. m2 quantiles for physiology paper figure
% use circnormed projection images 
%
%%% updated 2020-9-1 on thermaltake


% clear
% close all

%% setup parameters
filetable_name = 'F:\thesis\patchiness_anatomy_filelist.xlsx';
pathway_name = 'som'; %% somatostatin expression in V1
%     pathway_name = 'm2'; % secondary motor cortex proj to V1
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

run_plotting = 1; %%% after analysis,  run physpaper_anatomy_quants_plots.m
    
%% analysis
filetable = readtable(filetable_name);
filetable = filetable(strcmp(filetable.pathway, pathway_name),:); %%% keep only subjects in this pathway
filetable.Properties.VariableNames{strcmp(filetable.Properties.VariableNames, 'analyze_Case')} = 'analyze_case'; % correct irregular varname
filetable = filetable(logical(filetable.analyze_case),:); filetable.analyze_case = []; % keep only specified cases
nsubs = height(filetable); 
filetable.patchInterpatchRatio_projfilt = NaN(nsubs,1); % clear old results
filetable.patchInterpatchRatio_projorig = NaN(nsubs,1); % clear old results
filetable.PIRatio_baselined_normed = NaN(nsubs,1); % clear old results
filetable.patch_analysis = cell(nsubs,1);
filetable.proj_intens_by_quant = NaN(nsubs,findpatches_pars.nQuantiles); %%% mean projimage intensity in each m2 quantile
filetable.proj_intens_by_quant_normed = NaN(nsubs,findpatches_pars.nQuantiles); %%% mean projimage intensity in each m2 quantile
filetable.proj_shuffled_intens_by_quant_normed = NaN(nsubs,findpatches_pars.nQuantiles); %%% mean projimage intensity in each m2 quantile


for isub = 1:nsubs
    cd(filetable.directory{isub})
    filetable.patch_analysis{isub} = analyze_patchiness(filetable.m2_filt{isub}, filetable.proj_filt{isub}, filetable.roi{isub},...
        filetable.baseline{isub}, filetable.zoom(isub), filetable.scope{isub}, findpatches_pars, permutation_pars); % do main patch and projection analysis
    filetable.patchInterpatchRatio_projfilt(isub) = filetable.patch_analysis{isub}.permutation_test.patchInterpatchRatio; % store this subject's results
    %%%% converting to double results in negative values after subtracting baseline 
    projimage_baselined = double(filetable.patch_analysis{isub}.permutation_test.comparison_image) - filetable.patch_analysis{isub}.permutation_test.baseline; 
    
    
    
    %%% get quant intens for shuffled comparison image
    original_pix_area_squm = umPerPix(filetable.zoom(isub),filetable.scope{isub})^2;
    resample_factor = sqrt(original_pix_area_squm / resampled_pix_area_squm);
    roi_downsampled = imresize(filetable.patch_analysis{isub}.patchdata.roi, resample_factor);
    projimg_downsampled = imresize(projimage_baselined, resample_factor);
    randorder = randperm(length(projimg_downsampled(:))); % generate random order in which to reassign projimg_downsampled pixels
    projimg_shuffled = reshape(projimg_downsampled(randorder), size(projimg_downsampled)); % shuffle order of pixels
    
    
    for iquant = 1:findpatches_pars.nQuantiles % get mean input to each quantile for this subject
        this_quant_image = filetable.patch_analysis{isub}.patchdata.quantile_table.quantimage{iquant};
        filetable.proj_intens_by_quant(isub,iquant) = mean(projimage_baselined(this_quant_image)); % get mean intensity in this quant for this subject
        
        % analyze shuffled data
        this_quant_image_downsampled = imresize(this_quant_image, resample_factor); % downsample this quant image
        filetable.proj_shuffled_intens_by_quant_normed(isub,iquant) = mean(projimg_shuffled(this_quant_image_downsampled)); % get mean shuffled intensity in this quant for this subject
    end
    max_intens_this_sub = max(filetable.proj_intens_by_quant(isub,:));
    max_shuffled_intens_this_sub = max(filetable.proj_shuffled_intens_by_quant_normed(isub,:));
    filetable.proj_intens_by_quant_normed(isub,:) = filetable.proj_intens_by_quant(isub,:) ./ max_intens_this_sub;
    filetable.proj_shuffled_intens_by_quant_normed(isub,:) = filetable.proj_shuffled_intens_by_quant_normed(isub,:) ./ max_shuffled_intens_this_sub;
    if save_shuffled_proj_data
        save(shuffled_data_filename, 'filetable')
    end
end
filetable = movevars(filetable,'patchInterpatchRatio_projfilt','After','subject');

%% summary stats
mean_normed_intens_by_quant = mean(filetable.proj_intens_by_quant_normed);
sem_normed_intens_by_quant = std(filetable.proj_intens_by_quant_normed) ./ sqrt(nsubs);

%% run plotting
if run_plotting
    physpaper_anatomy_quants_plots();
end