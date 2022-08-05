% fluoro_jade_expt_1.m : script which searches through the xdb for
% qualifying experiments for this set of data.  It will create a sample
% graph plotting the basics  *in progress*
%
% 2007 Julian P. Meeks
%

fjc_expts = search_xdb(@query_fluoro_jade_expt_1);
number_of_expts = size(fjc_expts, 2);

timepoints = zeros(size(fjc_expts));
regionnames = cell(1); regionnames(1) = [];    % to hold region_names
cellsinregions = zeros(3, number_of_expts); % to hold cells_in_regions
regionareas = zeros(3, number_of_expts);   % to hold region_areas
exptconds = cell(1); exptconds(1) = [];     % to hold expt_condition s
cond_list = cell(1); cond_list(1) = [];
time_list = zeros(1); time_list(1) = [];

% 1) extract # of experimental time points (minutes in condition)
for idx_expt = 1:number_of_expts
    timepoints(idx_expt) = fjc_expts{idx_expt}.fluoro_jade.expt_timepoint;
end
%1a) create a list of the time points used
for idx_expt = 1:number_of_expts
    if isempty(time_list)
        time_list(1) = timepoints(idx_expt);
    elseif (timepoints(idx_expt) == time_list(:))==0
        time_list(end+1) = timepoints(idx_expt);
    end
end
time_list = sort(time_list);
%1b index experiments by time point
timepoint_index = cell(1); timepoint_index(1)=[];
for idx_tp = 1:size(time_list,2)
    timepoint_index(idx_tp) = {find(time_list(idx_tp)==timepoints(:))};
end

% 2) extract the region_names at play
for idx_expt = 1:number_of_expts
    n_samples = size(fjc_expts{idx_expt}.hist_analysis.region_names, 2);
    for idx_samples = 1:n_samples
        n_regions = size(fjc_expts{idx_expt}.hist_analysis.region_names{idx_samples},2);
        for idx_regions = 1:n_regions
            thisregion = fjc_expts{idx_expt}.hist_analysis.region_names{idx_samples}{idx_regions};
            if isempty(regionnames)
                regionnames(1) = thisregion;
            elseif isempty(strmatch(lower(thisregion), lower(regionnames)))
                regionnames(end+1) = thisregion;
            end
        end                
    end
end

% 3) extract the condition labels for each expt
for idx_expt = 1:number_of_expts
    exptconds(idx_expt) = {fjc_expts{idx_expt}.fluoro_jade.expt_condition};
end
% 3a) make a list of the unique expt_conditions
for idx_expt = 1:number_of_expts
    if isempty(cond_list)
        cond_list(1) = exptconds(idx_expt);
    elseif isempty(strmatch(exptconds(idx_expt), cond_list))
        cond_list(end+1) = exptconds(idx_expt);
    end
end
% 3b) index by experimental condition
condition_index = cell(1); condition_index(1)=[];
for idx_cond = 1:size(cond_list,2)
    condition_index(idx_cond) = {strmatch(cond_list{idx_cond}, exptconds, 'exact')};
end

% 4) extract the total points and area for each sample, each region
for idx_expt = 1:number_of_expts
    n_samples = size(fjc_expts{idx_expt}.hist_analysis.cells_in_regions, 2);
    for idx_samples = 1:n_samples
        n_regions = size(fjc_expts{idx_expt}.hist_analysis.cells_in_regions{idx_samples},2);
        cell_tally = zeros(size(regionnames));
        area_tally = zeros(size(regionnames));
        for idx_regions = 1:n_regions
            if cell_tally(idx_regions) == 0
                cell_tally(idx_regions) = fjc_expts{idx_expt}.hist_analysis.cells_in_regions{idx_samples}(idx_regions);
            else
                cell_tally(idx_regions) = cell_tally(idx_regions)+...
                    fjc_expts{idx_expt}.hist_analysis.cells_in_regions{idx_samples}(idx_regions);
            end
            if area_tally(idx_regions) == 0
                area_tally(idx_regions) = fjc_expts{idx_expt}.hist_analysis.region_areas{idx_samples}{idx_regions};
            else
                area_tally(idx_regions) = area_tally(idx_regions)+...
                    fjc_expts{idx_expt}.hist_analysis.region_areas{idx_samples}{idx_regions};
            end
        end
        if cellsinregions(idx_regions,idx_expt) == 0;
            cellsinregions(:,idx_expt) = cell_tally;
        else
            cellsinregions(:, idx_expt) = cellsinregions(:, idx_expt)+cell_tally(:);
        end
        if regionareas(idx_regions,idx_expt) == 0;
            regionareas(:,idx_expt) = area_tally;
        else
            regionareas(:, idx_expt) = regionareas(:, idx_expt)+area_tally(:);
        end
    end
end

% 5) assign the cells per area for each experiment

% 5a) make an array of the slice thickness values
slicethickness = zeros(n_regions, 1);
for idx_expt=1:number_of_expts
    slicethickness(:,idx_expt) = fjc_expts{idx_expt}.section.thickness;
end
% 5b make an array of the area per pixel (in square microns)
pixelareas = zeros(n_regions, 1);
for idx_expt = 1:number_of_expts
    pixelareas(:,idx_expt) = fjc_expts{idx_expt}.hist_analysis.um_per_pixel{1};
end
% 5c compile cellsperarea into its final version (in approximate cells/mm3) and
% normalized for average slice thickness
cellsperarea = cellsinregions./(regionareas.*pixelareas)./slicethickness*1e+9;

% 6 calculate the average pixel intensity in each region
% 6a sum pixel intensity in masked images
pixel_intensities = zeros(n_regions,number_of_expts);
for idx_expt = 1:number_of_expts
    sample_intensity = zeros(3,1);
    masked_imgs = jpm_region_mask(fjc_expts{idx_expt});
    for idx_samples = size(fjc_expts{idx_expt}.data_locations,2)
        for idx_regions = 1:n_regions
            sample_intensity(idx_regions)=sample_intensity(idx_regions)+sum(sum(masked_imgs{idx_samples}{idx_regions},2),1);
        end
    end
    pixel_intensities(:,idx_expt)=sample_intensity(:);
end
% 6b divide by the pixel area (final result in intensity/pixel)
intensityperpixel = pixel_intensities./regionareas;
% temporary correction for a single 8-bit sample
intensityperpixel(:,5)=intensityperpixel(:,5)*(2^16)/(2^8);

% 7) Calculate mean and standard errors for each conditions (by timepoint).
for idx_cond = 1:size(cond_list,2)
    for idx_tp = 1:size(time_list,2)
        expt_idx = intersect(condition_index{idx_cond},timepoint_index{idx_tp});
        for idx_region = 1:n_regions
            mean_cellsperarea(idx_cond, idx_tp, idx_region) = mean(cellsperarea(idx_region,expt_idx));
            stderr_cellsperarea(idx_cond, idx_tp, idx_region) = std(cellsperarea(idx_region, expt_idx))...
                                                                /sqrt(size(expt_idx,1));
            mean_pixelintensities(idx_cond, idx_tp, idx_region) = mean(intensityperpixel(idx_region, expt_idx));
            stderr_pixelintensities(idx_cond,idx_tp,idx_region) = std(intensityperpixel(idx_region,expt_idx))...
                                                                  /sqrt(size(expt_idx,1));
        end
    end
end
%% Begin plotting/figuring

mainfig = figure(1); hold on;

for idx_cond = 1:size(cond_list,2)
    region_ctr = 1; % sloppy way to skip a region if desired
    for idx_region = 1:2:n_regions % skipping region 2 (LOT, granule) for now
        switch idx_region
            case 1
                symbol = 'o';
            case 2
                symbol = 'o';
                % plot(time_list(:), mean_cellsperarea(idx_cond, :, 1), symbol);
        end
        switch idx_cond
            case 1
                color = 'b';
            case 2
                color = 'r';
            case 3
                color = 'm';
            case 4
                color = 'g';
            case 5
                color = 'k';
        end
        subplot(2,2,2+region_ctr); hold on; % beware, may be changed when sloppy-plotting
        plot_handles(idx_cond,idx_region) = errorbar(time_list(:), ...
                 mean_cellsperarea(idx_cond, :,idx_region), ...
                 stderr_cellsperarea(idx_cond, :, idx_region),...
                 [color symbol]);
        set(plot_handles(idx_cond, idx_region), 'linewidth', 2);
        title_handles(idx_region) = title(gca, regionnames{idx_region}, ...
                                    'fontweight', 'bold', 'fontsize', 12);
%         subplot(3,3, 6+idx_region); hold on;
%         plot_handles(idx_cond,idx_region) = errorbar(time_list(:), ...
%                  mean_pixelintensities(idx_cond, :,idx_region), ...
%                  stderr_pixelintensities(idx_cond, :, idx_region),...
%                  [color symbol]);
%         set(plot_handles(idx_cond, idx_region), 'linewidth', 2);
%         title_handles(idx_region) = title(gca, regionnames{idx_region}, ...
%                                     'fontweight', 'bold', 'fontsize', 12);
    region_ctr = region_ctr+1;  % comment out if not using sloppy plotting
    end

end

these_axes = findobj(mainfig, 'type', 'axes');
set(these_axes, 'xtick', time_list);
xlabels = cell2mat(get(these_axes, 'xlabel'));
set(xlabels, 'string', 'time (minutes)', 'fontsize', 12, 'fontname', 'helvetica');
ylabels = cell2mat(get(these_axes, 'ylabel'));
set(ylabels, 'string', 'degenerating neurons/mm^2', 'fontname', 'helvetica', 'fontsize', 12);
legend(findobj(these_axes,'style', 'line'), cond_list, 'fontname', 'helvetica',...
               'fontsize', 8, 'location', 'southeast');
           
 set(mainfig, 'units', 'pixels', 'position', [0 0 1200 800]);
 set(mainfig, 'units', 'normalized', 'color', [1 1 1]);
 
% % Set up transcardial control image [fjc_expts{27}.data_locations{3} = 2007/10/01 samp5 sec6]
% hold on;
% img1 = imread(fjc_expts{27}.data_locations{3});
% filt_img1 = imfilter(img1, punctate_filter_jpm_1);
% median_img1 = double(median(median(filt_img1(:),2),1));
% max_img1 = double(max(max(filt_img1(:),[],2),[],1));
% control_img_axis = axes('parent', mainfig, 'units', 'normalized', 'position', [.01 .51 .45 .45]);
% control_img_handle = imshow(filt_img1, 'parent', control_img_axis);
% set(control_img_axis, 'clim', [0.8*median_img1 max_img1]);
% control_img_title = title(control_img_handle, 'string', 'Transcardial Fixation Control',...
%                                               'fontsize', 12, 'fontweight', 'bold');
% % --------end control image setup-----------------------------------------

%% Set up dissection only control image [fjc_expts{15}.data_locations{3} = 2008/04/07 sample 'e' sec7]
hold on;
img1 = imread(fjc_expts{15}.data_locations{2});
filt_img1 = imfilter(img1, punctate_filter_jpm_1);
median_img1 = double(median(median(filt_img1(:),2),1));
max_img1 = double(max(max(filt_img1(:),[],2),[],1));
filt_img1([1700:1720], [300:300+100/0.71]) = 65535;
control_img_axis = axes('parent', mainfig, 'units', 'normalized', 'position', [.01 .51 .45 .45]);
control_img_handle = imshow(filt_img1);%, 'parent', control_img_axis);
set(control_img_axis, 'clim', [0.5*median_img1 0.9*max_img1], 'xdir', 'reverse', 'visible', 'off');
control_img_title = title(control_img_handle, 'string', 'Dissection Only Control',...
                                              'fontsize', 12, 'fontweight', 'bold');
% --------end control image setup-----------------------------------------

%% Set up chamber image [fjc_expts{1}.data_locations{2} = 2007/10/06 samp2 sec5]
hold on;
img2 = imread(fjc_expts{1}.data_locations{2});
filt_img2 = imfilter(img2, punctate_filter_jpm_1);
median_img2 = double(median(median(filt_img2(:),2),1));
max_img2 = double(max(max(filt_img2(:),[],2),[],1));
filt_img2([1700:1720], [400:400+100/0.71]) = 65535;
chamber_img_axis = axes('parent', mainfig, 'units', 'normalized', 'position', [.28 .51 .45 .45]);
chamber_img_handle = imshow(filt_img2);% , 'parent', chamber_img_axis);
set(chamber_img_axis, 'clim', [0.7*median_img2 0.9*max_img2]);
chamber_img_title = title(chamber_img_handle, 'string', '6 hr. Physiology Chamber',...
                                              'fontsize', 12, 'fontweight', 'bold');
% --------end chamber image setup-----------------------------------------

%% Set up bench image [fjc_expts{12}.data_locations{2} = 2007_10_06 samp3 sec5]
hold on;
img3 = imread(fjc_expts{12}.data_locations{2});
filt_img3 = imfilter(img3, punctate_filter_jpm_1);
median_img3 = double(median(median(filt_img3(:),2),1));
max_img3 = double(max(max(filt_img3(:),[],2),[],1));
filt_img3([1750:1770], [100:100+100/0.71]) = 65535;
bench_img_axis = axes('parent', mainfig, 'units', 'normalized', 'position', [.54 .51 .45 .45]);
bench_img_handle = imshow(filt_img3);% , 'parent', bench_img_axis);
set(bench_img_axis, 'clim', [0.7*median_img3 0.9*max_img3]);
bench_img_title = title(bench_img_handle, 'string', '6 hr. Warmed Bench',...
                                              'fontsize', 12, 'fontweight', 'bold');
% --------end bench image setup-----------------------------------------

%% Align Images
align([control_img_axis, chamber_img_axis, bench_img_axis], 'distribute', 'center');

%% Add Figure marks, etc.

% Add labels to control image
glom_label = text('parent', control_img_axis, 'units', 'normalized',...
                           'string', [{'glomerular'}; {'layer'}],'fontsize', 10,...
                           'fontweight', 'bold', 'color', [1 1 1], ...
                           'position', [.27 .6], 'horizontalalignment', 'center');
mit_label = text('parent', control_img_axis, 'units', 'normalized',...
                           'string', [{'mitral cell'}; {'layer'}],'fontsize', 10,...
                           'fontweight', 'bold', 'color', [1 1 1], ...
                           'position', [.47 .4], 'horizontalalignment', 'center');
gran_label = text('parent', control_img_axis, 'units', 'normalized',...
                           'string', [{'granule cell'}; {'layer'}],'fontsize', 10,...
                           'fontweight', 'bold', 'color', [1 1 1], ...
                           'position', [.72 .6], 'horizontalalignment', 'center');
                       
% Add Panel Identifiers
panel_A_text = text('parent', control_img_axis,'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', 'A',...
                    'position', [0 1.06]);
panel_B_text = text('parent', chamber_img_axis,'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', 'B',...
                    'position', [0 1.06]);
panel_C_text = text('parent', bench_img_axis,'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', 'C',...
                    'position', [0 1.06]);
panel_D_text = text('parent', these_axes(2),'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', 'D',...
                    'position', [-.07 1.1]);
panel_E_text = text('parent', these_axes(1),'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', 'E',...
                    'position', [-.07 1.1]);
               
% Add asterisks
astk(1) = text('parent', these_axes(2),'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', '*',...
                    'position', [.63 .45]);
astk(2) = text('parent', these_axes(2),'units','normalize', 'fontname', 'helvetica',...
                    'fontsize', 18, 'fontweight', 'bold', 'string', '*',...
                    'position', [.9 .45]);
                
%% T-Tests

% Unpaired T-TEST SCRIPTS   
% (for the 240 minute timepoint: chamber vs. chamber_suction (col1 =
% glomerular/mitral, col2 = LOT, granule, col3 = granule)
[h,p] = ttest2(rot90(cellsperarea(:,intersect(condition_index{1},timepoint_index{4}))),...
              rot90(cellsperarea(:,intersect(condition_index{2},timepoint_index{4}))),...
              0.05, 'both');
% (for the 240 minute timepoint: chamber_suction vs. dissection only (col1 = glomerular/mitral, col3 = granule
[h,p] = ttest2(rot90(cellsperarea(:,intersect(condition_index{2},timepoint_index{4}))),...
              rot90(cellsperarea(:,intersect(condition_index{3},timepoint_index{1}))),...
              0.05, 'both');
% (for the 240 minute timepoint: chamber_suction vs. w. bench (col1 = glomerular/mitral, col2 = LOT, granule                
[h,p] = ttest2(rot90(cellsperarea(:,intersect(condition_index{2},timepoint_index{4}))),...
              rot90(cellsperarea(:,intersect(condition_index{4},timepoint_index{4}))),...
              0.05, 'both');
% (for the 240 minute timepoint: chamber vs. w. bench (col1 = glomerular/mitral, col2 = LOT, granule                
[h,p] = ttest2(rot90(cellsperarea(:,intersect(condition_index{1},timepoint_index{4}))),...
              rot90(cellsperarea(:,intersect(condition_index{4},timepoint_index{4}))),...
              0.05, 'both');
% (for the 240 minute timepoint: chamber vs. w. bench (col1 = glomerular/mitral, col2 = LOT, granule                
[h,p] = ttest2(rot90(cellsperarea(:,intersect(condition_index{1},timepoint_index{5}))),...
              rot90(cellsperarea(:,intersect(condition_index{4},timepoint_index{5}))),...
              0.05, 'both');