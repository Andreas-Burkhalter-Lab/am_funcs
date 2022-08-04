%%%% cell and dendrite density analysis
%
% approximate dendrite lengths within quantiles by measuring area of manually traced dendrites
%
% manually traced dendrite lines are assumed to be 1 pixel thick, number of pixels = length of segments in pixels
% note: using 1-point stroke weight for lines in Adobe Illustrator with 1292x1040 images results in 1-pixel-thick lines
%
%%% updated 2021/3/13 on thermaltake


filelist_excel = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\dendrite_lengths_file_list.xlsx';
close all

%% plotting parameters
% analysis_type = 'cells'; 
%     path_to_plot = 'amyg_cells';
%     path_to_plot = 'mec_cells';
analysis_type = 'dendrites';
%    path_to_plot = 'amyg';
   path_to_plot = 'mec';

boxplot_linewidth = 2; 
bar_face_color = [0.3 0.3 0.3]; % bar graph -  bar color
bar_line_width = 2; % bar graph - bar outline line width
error_line_color = [0 0 0]; % error bars color
axes_line_width = 2;
axes_numbers_bold = 'normal'; %%% 'bold' or 'normal'
plotfont = 'Arial'; % default = Helvetica
axis_font_size = 18; 
% ylimits = [0.85, 1.89];
% ylimits = [0 0.4]; 2
% ytick = [0 : 0.1 : 0.4];
xlimits = [0.3, 6.6];
% if aligned with 20014 dendrites with scatter plot, use 770 width, 300 height
fig_width = 770; % pixels... if aligned with 19105_s2d_x20, use 1200 width, 430 height
output_cor_results = 0; % display results on the command line
fig_height = 300; % pixels... if aligned with 19114_ctx_s3c_x10_cropped, use 1100 width, 430 height
add_scatter_points = 1; % display all individual data points overlayed with bar graphs; remove bar fill shading
    max_x_deviation = 0.3; % maximum x spacing jitter
    scatter_color = [0.5 0.5 0.5]; 
    scatter_marker_size = 70;
save_dendrite_cells_fig = 0; 
    savename_dendrite_lengths = ['fig_', path_to_plot, '_dendrite_lengths'];
    savename_cells_by_quant = ['fig_', path_to_plot, '_by_quant'];    
perform_analysis = 1;  % turn off if data is already loaded and analyzed
    

%% findpatches settings
findpatches_pars.nQuantiles = 1;

findpatches_pars.minAreaPatch_squm = 0;  % min area in square microns a blob must contain to be considered a patch
findpatches_pars.maxAreaPatch_squm = field_default(findpatches_pars,'maxAreaPatch_squm',inf);
findpatches_pars.blur_using_only_roi = field_default(findpatches_pars,'blur_using_only_roi',1); % if true, do not use non-roi pixels for blurring 
findpatches_pars.diskblurradius_um = field_default(findpatches_pars,'diskblurradius_um',29);
findpatches_pars.threshMode = field_default(findpatches_pars,'threshMode','intensityQuantiles');
    findpatches_pars.nQuantiles = field_default(findpatches_pars,'nQuantiles',6); % number of intensity levels to divide image into (if using 'intensityQuantiles')
    findpatches_pars.patchQuantiles = field_default(findpatches_pars,'patchQuantiles',[4 5 6]); % % quantiles selected to be counted as patches(if using 'intensityQuantiles')
findpatches_pars.include_interior_nonroi_in_roi = field_default(findpatches_pars,'include_interior_nonroi_in_roi',1); 
findpatches_pars.show_plots = field_default(findpatches_pars,'show_plots',0);
findpatches_pars.save_patchdata = 0; 


%% analyze dendrite lengths for each case
if perform_analysis
%     clear
    filelist = readtable(filelist_excel); % import list of images to analyze
    filelist = filelist(strcmp(filelist.analysis_type, analysis_type), :); % keep only cases corresponding to the analysis type selected above
    filelist = filelist(strcmp(filelist.pathway, path_to_plot), :); % keep only cases corresponding to the pathway selected above    
    nfiles = height(filelist);
    for ifile = 1:nfiles
        this_roi_image = loadbw(filelist.roi_file{ifile}); % load the area roi
        % get m2 quantiles data for this case
        filelist.patchdata{ifile} = findpatches(filelist.m2_file{ifile},filelist.roi_file{ifile},filelist.zoom(ifile),filelist.scope{ifile}, findpatches_pars); 
        if strcmp(analysis_type,'dendrites')
            this_dend_image = loadbw(filelist.dendrite_file{ifile}); % load the labeled dendrites image
            this_dend_image = this_dend_image & this_roi_image; % only analyze dendrites that fall within the roi
            total_dend_pix_this_case = nnz(this_dend_image); 
            filelist.total_dend_pix_this_case(ifile) = total_dend_pix_this_case; 
            for quant = 1:findpatches_pars.nQuantiles
                this_quant_area = filelist.patchdata{ifile}.quantile_table.quantimage{quant};
                n_dend_pix_this_quant = nnz(this_dend_image & this_quant_area); %%% number of labeled dendrite pixels that overlap with this quantile
                filelist.dend_um_per_quant(ifile, quant) = umPerPix(filelist.zoom(ifile),filelist.scope{ifile}) * n_dend_pix_this_quant; % convert pix to um
                filelist.dend_proportion_per_quant(ifile, quant) = n_dend_pix_this_quant / total_dend_pix_this_case; % proportion of total labeled dendrites from this case in this quantile
            end
        elseif strcmp(analysis_type,'cells')
            
           [filelist.chisq_pval(ifile), quant_cells_table, filelist.cell_centers_img{ifile}] = chi_square_quantiles_cells(filelist.patchdata{ifile}, filelist.cells_file{ifile});
            total_cells_this_case = sum(quant_cells_table.ncells);
            for quant = 1:findpatches_pars.nQuantiles
                this_quant_area = filelist.patchdata{ifile}.quantile_table.quantimage{quant};
                filelist.cells_per_quant(ifile, quant) = quant_cells_table.ncells(quant);
                filelist.cells_proportion_per_quant(ifile, quant) = quant_cells_table.ncells(quant) / total_cells_this_case; 
            end
        end
    end  
end


%% plotting
% get data to plot
rows_to_plot = strcmp(filelist.pathway,path_to_plot);
if strcmp(analysis_type,'dendrites')
    data_to_plot = [filelist.dend_proportion_per_quant(rows_to_plot,:)]'; % dend lengths for this pathway
elseif strcmp(analysis_type,'cells')
    data_to_plot = [filelist.cells_proportion_per_quant(rows_to_plot,:)]';                 % cell proportions for this pathway
end
groupstats = table; [groupstats.mean, groupstats.sem] = grpstats(data_to_plot(:),repmat([1:size(data_to_plot,1)]',nnz(rows_to_plot),1),{'mean','sem'}); % organize data
fig_dendrites_cells = figure;
% plot proportion of cells/dendrite lengths in each quantile with standard error
nquants = height(groupstats); 
bg = bar(1:height(groupstats),groupstats.mean);
hold on

% show data from individual cases
if add_scatter_points
    mgrid = meshgrid(1:size(data_to_plot,1),1:1:size(data_to_plot,2))'; % pre-jittered indices indices
    noisegrid = 2*max_x_deviation * rand(size(data_to_plot)) - max_x_deviation; % x jitter falues
    x_indices = mgrid + noisegrid; % jittered indices
    scatplot = scatter(x_indices(:), data_to_plot(:)); 
    scatplot.SizeData = scatter_marker_size; 
    scatplot.MarkerEdgeColor = scatter_color;
    scatplot.MarkerFaceColor = scatter_color;
    bar_face_color = [1 1 1]; % white bar face
end

eb = errorbar(1:nquants, groupstats.mean, groupstats.sem,'LineStyle','none');
hold off
plot_formatting()
hax = gca;
hax.XAxis.TickLength = [0 0];
hax.XTick = 1:nquants;
hax.XTickLabel = num2str([1:nquants]');


%%% correlation stats
[~,quantmat] = meshgrid(1:size(data_to_plot,2),1:size(data_to_plot,1));
[r_corr_quants, p_corr_quants] = corrcoef(quantmat,data_to_plot);
if output_cor_results
    r_corr_quants
    p_corr_quants
end
    
%%% save figure
set(fig_dendrites_cells,'Renderer', 'painters', 'Position', [200 200 fig_width fig_height]) % set figure length and width
if save_dendrite_cells_fig
    if strcmp(analysis_type,'dendrites')
%         print(fig_dendrites_cells,savename_dendrite_lengths, '-dtiffn','-r300') %%% save image as file
        saveas(fig_dendrites_cells,savename_dendrite_lengths, 'svg') %%% save image as file
    elseif strcmp(analysis_type,'cells')
%         print(fig_dendrites_cells,savename_cells_by_quant, '-dtiffn','-r300') %%% save image as file
        saveas(fig_dendrites_cells,savename_cells_by_quant, 'svg') %%% save image as file
    end
end
