%%%% called by toplevel_quants_plotting.m
%
%%% updated 2020-11-26 on thermaltake


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
    filetable.patch_analysis{isub} = analyze_patchiness(filetable.m2_file{isub}, filetable.proj_file{isub}, filetable.roi{isub},...
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
    %%% choose shuffled or nonshuffled data
    if ~plot_shuffled_data
        projdata_by_quant = filetable.proj_intens_by_quant_normed;
    elseif plot_shuffled_data
        projdata_by_quant = filetable.proj_shuffled_intens_by_quant_normed;
    end

    % output correlation coefficient and p value for quant number vs egfp intensity
    egfp_by_quant = projdata_by_quant; 
    nquants = size(egfp_by_quant,2);
    quantvals = repmat(1:nquants,height(filetable),1);
    [r, p] = corrcoef(quantvals, egfp_by_quant); 
        fprintf(['r = ', num2str(r(2,1)), '\n']); fprintf(['p = ', num2str(p(2,1)), '\n']);

    %% make figure
    % % % % % % hbar = bar(mean(filetable.proj_intens_by_quant_normed), 'FaceColor',bar_color, 'LineWidth',bar_border_width);

    xvals = 1:6;
    if plot_group_bars
        quants_bar_fig = figure; 
        mean_proj_vals = mean(projdata_by_quant,1);
        hbar = bar(xvals, diag(mean_proj_vals), 'stacked');
        hold on
        sem = std(projdata_by_quant,[],1) ./ sqrt(nsubs);
        errorbar(mean_proj_vals,sem,'.','Color',ebar_color,'LineWidth',ebar_linewidth)
        plot_formatting; % format the plot
        set(gca,'XTick',[1:6])
        set(gca,'XTickLabels',{'1','2','3','4','5','6'})  %%%% restore x tick labels
        hax = gca;
        hax.XAxis.TickLength = x_tick_length;
        hax.FontSize = axes_font_size;
        hax.LineWidth = axes_line_width; 
        % % % % % % hbar.EdgeColor = bar_border_color;
        for ibar = 1:6
            set(hbar(ibar),'facecolor',bar_colors(ibar,:))
            set(hbar(ibar),'EdgeColor',bar_border_color);
        end
        set(quants_bar_fig,'Renderer', 'painters', 'Position', [10 10 fig_width fig_height]) % set figure dimensions
        set(gca,'Box','off')
        if plot_shuffled_data
            if show_title; title('shuffled'); end
            savename_quants_bar_fig = [savename_quants_bar_fig, '_shuffled']; 
        end
    end
    if plot_individual_bars
        figure
        for isub = 1:nsubs
            subplot(subplot_rowcol(2), subplot_rowcol(1), isub)
            hbar = bar(xvals, diag(projdata_by_quant(isub,:)), 'stacked');
            title(num2str(filetable.subject(isub)))
        end
        if plot_shuffled_data; suptitle('shuffled'); end
    end

    if save_quants_bar_fig
        saveas(quants_bar_fig,savename_quants_bar_fig, 'svg') %%% save image as file
    end
    
end