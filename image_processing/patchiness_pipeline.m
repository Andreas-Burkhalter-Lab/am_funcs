%%%%%%%%%%%% patchiness_pipeline
%%%% script called by analyze_patchiness after analysis parameters are set
% last updated 2021/7/29 on thermaltake

scope = vardefault('scope','epimicro_old');

%% get patches
roifile = vardefault('roifile',[]);
patchdata = findpatches(m2file, roifile, zoom, scope, findpatches_pars);
if do_patch_density_analysis
    dens_analysis = patchDensity(patchdata.patchimage,roifile,zoom,scope);
    patchdata = copyAllFields(patchdata,dens_analysis);
end
if save_output
    save(output_filename,'patchdata');
end
allresults.patchdata = patchdata;

%% patch vs interpatch tests
if ~isempty(comparison_file)
    % permutation test
    if permutation_pars.custom_interpatches 
        if ~strcmp(findpatches_pars.threshMode,'intensityQuantiles')
            error('If you are using custom interpatches, threshMode must be set to ''intensityQuantiles''.')
        end
        interpatchimage = false(size(patchdata.patchimage));
        for i = 1:length(permutation_pars.interpatchQuantiles) % build interpatch image from the selected quantiles
            quant = permutation_pars.interpatchQuantiles(i);
            interpatchimage = interpatchimage | patchdata.quantile_table.quantimage{quant};
        end
    end
    baselineFile = vardefault('baselineFile',[]);
    permutation_test = permutation_test_patches(patchdata.patchimage, comparison_file, patchdata.roi, baselineFile, zoom, scope, permutation_pars);
    permutation_test;
    allresults.permutation_test = permutation_test;
    
    % perimeter test
    perimeter_test = perimeter_test_patches(patchdata.patchimage, comparison_file, permutation_test.roi, permutation_test.baselineimg, zoom, scope, perimeter_pars);
    perimeter_test;
    allresults.perimeter_test = perimeter_test;
    if save_output
        save(output_filename,'permutation_test','perimeter_test','-append')
    end
end

%% plotting
if fig_pars.show_plots
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    imagesc(double(patchdata.im))
    hold on
    [bnds_patches_y, bnds_patches_x] = find(patchdata.bnds_patches');
    scatplot1 = scatter(bnds_patches_y, bnds_patches_x, '.');
    scatplot1.SizeData = fig_pars.patch_border_width;
    scatplot1.MarkerEdgeColor = fig_pars.patch_border_color;   
    title('Patches image with patch borders')
    hold off
    if ~isempty(comparison_file)
        subplot(1,2,2)
        imagesc(double(permutation_test.comparison_image)); 
        hold on
        scatplot2 = scatter(bnds_patches_y, bnds_patches_x, '.');
        scatplot2.SizeData = fig_pars.patch_border_width;
        scatplot2.MarkerEdgeColor = fig_pars.patch_border_color;
        title('Comparison image with patch borders')
        if permutation_pars.custom_interpatches
            interpatchimage = false(size(patchdata.im));
            for i = permutation_pars.interpatchQuantiles
                interpatchimage(patchdata.quantile_table.quantimage{i}) = true;
            end
            bnds_interpatches = edge(interpatchimage);
            [bnds_interpatches_y, bnds_interpatches_x] = find(bnds_interpatches);
            hold on
            scatplot2 = scatter(bnds_interpatches_x, bnds_interpatches_y, '.');
            scatplot2.SizeData = fig_pars.interpatch_border_width;
            scatplot2.MarkerEdgeColor = fig_pars.interpatch_border_color;
            hold off
            subplot(1,2,1) % add interpatches to m2 image
            hold on
            scatplot2 = scatter(bnds_interpatches_x, bnds_interpatches_y, '.');
            scatplot2.SizeData = fig_pars.interpatch_border_width;
            scatplot2.MarkerEdgeColor = fig_pars.interpatch_border_color;
        end
        hold off
    end
end

%% filled plot
fig_pars.filled_quants_plot = field_default(fig_pars, 'filled_quants_plot', false);
if fig_pars.filled_quants_plot
    figure
    imagesc(patchdata.quant_levels_img)
    colormap(fig_pars.filled_colormap)
    cbar = colorbar;
    cbar.Label.String = 'quantile';
    cbar.Label.FontSize = 13; 
    colormap(fig_pars.filled_colormap)
end
