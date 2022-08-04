%%% laminar_cells_plotting: plot optical density laminar distribution to LA, MEC from POR from multiple cases
%
% interpolate values for missing sections, get means and standard error for each layer and pathway
%
% updated 2020/12/19 on thermaltake 

% close all



%     excel_rows = 2:76; % mec - 19105, 19106
%             legend_labels = {'LM-->EntM';'LI-->EntM';'P-->EntM';'POR-->EntM';'PORa-->EntM'};
            
% % %     excel_rows = 132:161; % amyg - mouse 16128
    excel_rows = 77:161; %%% amyg - mice 16127, 16128
                legend_labels = {'LA-->LM';'LA-->LI';'LA-->P';'LA-->POR';'LA-->PORa';};
            
            
analyze_data = 0; % turn off if analyzed data is already loaded
compare_por_L23_vs_L5 = 1; % run ttest on cells per section in L2/3 vs. L5
   L23_secs = [3:5]; % sections to be considered occupying Layer 2/3
   L5_secs = [9:12]; % sections to be considered occupying Layer 5
   por_pathnames = {'mecpor', 'amygpor'}; 
save_lam_density_fig = 0; 
    savename_lam_density = 'fig_cells_density';            
            
filelist = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\cellcounting_file_list.xlsx';
norm_within_path = 0; % normalize cell counts to max within mouse and path; if false, normalize to max of all paths in this mouse
% interp_method = 'spline';
%     interp_method = 'linear';
%     interp_method = 'pchip';
    interp_method = 'makima';
axis_font_size = 20; 
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
line_width = 4;
ebar_width = 1; % errorbar width
box_on_off = 'off';
section_thickness_um = 40; % for converting section number to depth
path_offset_um = 2; % shift pathways this much on plot so they don't overlap
horizontal_plot = 1;
fig_width = 1000;% pixels
fig_height = 950; % pixels
ylimits = [0 600];
xticks_hrzplot = linspace(0, 1, 5); % x ticks for horizontal plot
xticklabels_hrzplot = {'0','','0.5','','1'};
subplot_extra_spacing = 0.06; %%% space to add between subplots



colororder = [... % plot line colors
0    0    1;...
0.4660    0.6740    0.1880;...
1   0    0;...
0.4940    0.1840    0.5560;...
0.3010    0.7450    0.9330;...
0.6350    0.0780    0.1840;...
0.9290    0.6940    0.1250;...
];

sec_to_depth = @(sec,thick)thick*sec; %function for converting section to depth
%     sec_to_depth = @(sec,thick)thick*sec - thick; %function for converting section to depth

%  analyze data
if analyze_data
    filetable = countcells(filelist, excel_rows); %%% analyze and tabulate laminar density data
end

%% sort results by pathway and subject
unqpaths = unique(filetable.pathway,'stable');
npaths = length(unqpaths);
unqsubs = unique(filetable.sub,'stable');
nsubs = length(unqsubs);
lamdata = table(unqpaths,cell(npaths,1),cell(npaths,1), 'VariableNames', ...
                {'pathway', 'subjects', 'data'});
            
% % get max values within each mouse
mouse_table = table( unqsubs, NaN(nsubs,1), 'VariableNames', {'sub', 'max_cells'});
for isub = 1:nsubs
    matchrow = filetable.sub == mouse_table.sub(isub);
    mouse_table.max_cells(isub) = max(filetable.ncells(matchrow));
end

%%% analyze each path and case
for ipath = 1:npaths
    thispath = unqpaths{ipath};
    lamdata.pathway{ipath} = thispath; 
    this_path_rows = strcmp(filetable.pathway,thispath);
    thesesubs = unique(filetable.sub(this_path_rows));
    n_subs_this_path = length(thesesubs);
    lamdata.subjects{ipath} = thesesubs; 
    n_secs_this_path = max(filetable.section(this_path_rows)); % highest number section among cases in this pathway
    nancol = NaN(n_secs_this_path,1);
    nanmat = NaN(n_secs_this_path,n_subs_this_path); 
    lamdata.data{ipath} = table(nancol, nanmat,    nanmat,      nanmat,        nanmat,         nancol,               nancol,         nancol,     nancol,   'VariableNames', ...
                            {'section','depth_um','cells_normed','cells','cells_pre_interp','mean_cells_normed','sem_cells_normed','mean_cells','sem_cells'});
    lamdata.data{ipath}.section = [1:n_secs_this_path]';
    for isub = 1:n_subs_this_path
        this_sub_rows = this_path_rows & filetable.sub == thesesubs(isub); % rows match this subject and path
        this_sub_sections = filetable.section(this_sub_rows); 
        this_sub_ncells = filetable.ncells(this_sub_rows); 
        lamdata.data{ipath}.cells_pre_interp(this_sub_sections,isub) = this_sub_ncells; %%% add data to the table for this subjects' cellsities
        lamdata.data{ipath}.cells(:,isub) = round(interp1(this_sub_sections, this_sub_ncells, lamdata.data{ipath}.section, interp_method, NaN)); % interpolatemissing section cell numbers
        negrows = lamdata.data{ipath}.cells(:,isub) < 0;
        lamdata.data{ipath}.cells(negrows,isub) = 0; % set zero minimum
        max_cells_this_case_path = max(lamdata.data{ipath}.cells(:,isub)); % max value of all sections in this case in this path
        lamdata.data{ipath}.mean_cells = nanmean(lamdata.data{ipath}.cells,2);
        lamdata.data{ipath}.sem_cells = nanstd(lamdata.data{ipath}.cells,0,2) / sqrt(n_subs_this_path);
        if norm_within_path % norm within mouse and path
            lamdata.data{ipath}.cells_normed(:,isub) = lamdata.data{ipath}.cells(:,isub) ./ max_cells_this_case_path; %%% cells as a fraction of max
        else % norm within mouse
            max_cells_this_case = mouse_table.max_cells(mouse_table.sub == thesesubs(isub)); % max cells from this mouse
            lamdata.data{ipath}.cells_normed(:,isub) = lamdata.data{ipath}.cells(:,isub) ./ max_cells_this_case; %%% cells as a fraction of max
        end
        lamdata.data{ipath}.mean_cells_normed = nanmean(lamdata.data{ipath}.cells_normed,2);
        lamdata.data{ipath}.sem_cells_normed = nanstd(lamdata.data{ipath}.cells_normed,0,2) / sqrt(n_subs_this_path);
    end
    lamdata.data{ipath}.depth_um = sec_to_depth(lamdata.data{ipath}.section, section_thickness_um); % assign approximate depths
    lamdata.totalcells(ipath) = nansum(lamdata.data{ipath}.cells(:));
end

%% output L2/3 vs. L5 density
if compare_por_L23_vs_L5 % run ttest on cells per section in L2/3 vs. L5
   lamdata_por = lamdata.data{find(ismember(lamdata.pathway,por_pathnames),1)}; % get appropriate lamdata table
   L56_secs = min(L5_secs) : height(lamdata_por); % L5 sections to the deepest section
   L23_cells_normed = lamdata_por.cells_normed(L23_secs,:);
   L23_cells_normed_allsections = nanmean(L23_cells_normed,1) % combine across sections with L2/3
   L5_cells_normed = lamdata_por.cells_normed(L5_secs,:);
   L5_cells_normed_allsections = nanmean(L5_cells_normed,1); % combine across sections with L5
   L56_cells_normed = lamdata_por.cells_normed(L56_secs,:);
   L56_cells_normed_allsections = nanmean(L56_cells_normed,1) % combine across sections with L5
   [h p_L23_vs_L5] = ttest2(L23_cells_normed_allsections, L5_cells_normed_allsections);
   [h p_L23_vs_L56] = ttest2(L23_cells_normed_allsections, L56_cells_normed_allsections)
   L56_to_L23_ratios = L56_cells_normed_allsections ./ L23_cells_normed_allsections;
   L56_to_L23_ratios_mean = mean(L56_to_L23_ratios)
   L56_to_L23_ratios_sem = std(L56_to_L23_ratios)/sqrt(length(L56_to_L23_ratios))
end
   
 %% make plot           
fig_cells_density = figure;
clear hplot
hold on
for ipath = 1:npaths
%     hsubplot(ipath) = subplot(1,npaths,ipath);
    depth_offset = path_offset_um * [ipath-0.5-npaths/2]; % shift depth values this much for this path
    if horizontal_plot
        hplot(ipath,:) = herrorbar( lamdata.data{ipath}.mean_cells_normed, lamdata.data{ipath}.depth_um + depth_offset, lamdata.data{ipath}.sem_cells_normed);%
        hplot(ipath,1).Color = colororder(ipath,:);
        hplot(ipath,2).Color = colororder(ipath,:);
%         hplot(1).Color = [0 .6 0];
%         hplot(2).Color = [0 .6 0];
        hplot(ipath,1).LineWidth = ebar_width;
        hplot(ipath,2).LineWidth = line_width;
        plot_formatting()
        set(gca,'YDir','reverse')
%         xlabel('Optical density')
        set(gca,'XTick',xticks_hrzplot)
        set(gca,'XTickLabels',xticklabels_hrzplot)
%         ylabel('Depth (µm)')
        if ipath == npaths
            hleg = legend(hplot(:,2),legend_labels);
        end
    elseif ~horizontal_plot
        errorbar(lamdata.data{ipath}.depth_um + depth_offset, lamdata.data{ipath}.mean_cells_normed, lamdata.data{ipath}.sem_cells_normed, 'Color',colororder(ipath,:),'LineWidth',line_width)
        %     ylabel('Normalzied cell number')
        xlabel('Depth (µm)')
        hleg = legend(legend_labels);
    end

    set(gca,'defaultAxesColorOrder',colororder)


end

% add spacing betewen subplots
% % % % % % % % % % % % % % hsubplot(1).Position(1) = hsubplot(1).Position(1) - subplot_extra_spacing;
% % % % % % % % % % % % % % hsubplot(3).Position(1) = hsubplot(3).Position(1) + subplot_extra_spacing;


set(fig_cells_density,'Renderer', 'painters', 'Position', [50 50 fig_width fig_height]) % set figure length and width
if save_lam_density_fig
    print(fig_cells_density,savename_lam_density, '-dtiffn','-r300') %%% save image as file
end
