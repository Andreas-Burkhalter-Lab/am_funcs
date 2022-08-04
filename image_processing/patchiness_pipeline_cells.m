%%%%%%%%%%%% patchiness_pipeline
%%%% script called by analyze_patchiness after analysis parameters are set
% last updated 18/7/11 on thermaltake

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
res.patchdata = patchdata;

%% chi square test of cells in different quantiles
if strcmp(findpatches_pars.threshMode, 'intensityQuantiles')
    [cells_image cells_file] = procImageInput('cells_file');
    res.cells_image = cells_image;
    res.cells_file = cells_file;
    [res.chisq_pval, res.quant_cells_table, res.cell_centers_img] = chi_square_quantiles_cells(patchdata, cells_image);
    if save_output
        save(output_filename,'res');
    end
    res.chisq_pval
else
    fprintf('\nUse findpatches_pars.threshMode==''intensityQuantiles'' to analyze cells vs. patches.\n\n')
end

%% plotting
if fig_pars.show_plots && strcmp(findpatches_pars.threshMode, 'intensityQuantiles')
    % m2 and quantiles
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1)
    imagesc(patchdata.imroi -  min(min(patchdata.imroi(patchdata.imroi>0)))); % show baselined image, roi only
    colormap('gray')
    hold on
    for i = 1:patchdata.pars.nQuantiles
        [dotsy dotsx] = ind2sub(size(patchdata.im), find(patchdata.quant_edges_img==i)); 
        scatplot = scatter(dotsx, dotsy,'.');
        scatplot.SizeData = fig_pars.quant_border_width;
    end
    title('m2 image and quantiles')
        
    % cells and quantiles
    subplot(1,2,2)
    imagesc(res.cell_centers_img); % 
    colormap('gray')
    hold on
%     quants_to_plot = 1:patchdata.pars.nQuantiles;
    quants_to_plot = patchdata.pars.nQuantiles - findpatches_pars.topNQuantilesArePatches + 1; % only plot patch borders
    for i = quants_to_plot
        [dotsy dotsx] = ind2sub(size(patchdata.im), find(patchdata.quant_edges_img==i)); 
        scatplot = scatter(dotsx, dotsy,'.');
        scatplot.SizeData = fig_pars.quant_border_width;
    end
    [dotsy dotsx] = ind2sub(size(patchdata.im), find(res.cell_centers_img));
    scatplot = scatter(dotsx, dotsy,'.','w');
    title('patch borders and cells')
    
    % quantiles bar graph
    figure
    bar(res.quant_cells_table.cells_per_sqmm)
    ylabel('cell per sq micron')
    xlabel('quantile (greater value == closer to patch center)')
end