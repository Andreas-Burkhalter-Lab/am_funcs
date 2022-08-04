%%% module_ratio_plotting: plot patch:interpatch ratios for all cases
% updated 2020/3/1

% close all

%% parameters
filelist_name = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\module_input_ratios_filelist.xlsx';
basevalue = 1; %%% baseline value for bar graph
bar_face_color = [0.3 0.3 0.3]; % bar graph -  bar color
bar_line_width = 2; % bar graph - bar outline line width
error_line_color = [0 0 0]; % error bars color
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
y_log_or_linear = 'log'; % yscale style
plotfont = 'Arial'; % default = Helvetica
axis_font_size = 13; 
xlimits = [0.3, 15.5];
ylimits = [0.34, 6.5];
ytick = [0.5, 1, 2, 4]; 
fig_width = 550; % pixels
fig_height = 500; % pixels
plot_module_ratios = 1; 
    save_module_ratios = 0; % save figure as tiff file
        savename_module_ratios = 'fig_module_ratios';
paths_to_plot = {... % plot cases from filelist with these pathway names
%     'dlgn_v1';...
%     'lp_v1';...
    'dlgn_por';...
    'lp_por';...
    'amyg_por';...
    'v1_por';...
    'dlgn_lm';...
    'lp_lm';...
    'dlgn_li';...
    'lp_li';...
    };
%%% specify how to group pathways on the plot
path_plot_groups = {...
%                     {'dlgn_v1','lp_v1'};...
                    {'dlgn_lm','lp_lm'};...
                    {'dlgn_li','lp_li'};...
                    {'dlgn_por','lp_por','amyg_por','v1_por'}...
                    };
plot_color_indices = {...
%                       [1, 2];... %%% bargraph colors corresponding to each group/path
                      [1, 2];...
                      [1, 2];...
                      [1, 2, 3, 4];...
                     };
colororder = [... % plot line colors to plug into plot_color_indices
    0    0    1;...
    0.4660    0.6740    0.1880;...
    1   0    0;...
    0.4940    0.1840    0.5560;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840;...
    0.9290    0.6940    0.1250;...
    ];
                
%% organize data to plot
cases = readtable(filelist_name); %%% import module ratio data
pathlist = unique(cases.pathway,'stable'); 
stats_module_ratios = grpstats(cases, {'pathway'}, {'mean','sem'}, 'DataVars','patchInterpatchRatio_projfilt');
[~, rows_to_keep] = intersect(stats_module_ratios.pathway,paths_to_plot); % get the data rows to plot
rows_to_keep = sort(rows_to_keep);
stats_to_plot = stats_module_ratios(rows_to_keep,:); % discard rows we aren't plotting
npaths_to_plot = height(stats_to_plot);
ngroups = size(path_plot_groups,1); 
%%% set bargraph positioning and colors
xpos = -1; % plot x position
% assign plot groups to data
for igroup = 1:ngroups
    xpos = xpos+1; % add spacer between groups
    for index_in_group = 1:length(path_plot_groups{igroup}) % index withing group
        xpos = xpos+1; % move to next x position on bar graph
        irow = find(strcmp(path_plot_groups{igroup}{index_in_group}, stats_to_plot.pathway)); % row within stats_to_plot
        stats_to_plot.plotgroup(irow) = igroup; 
        stats_to_plot.xpos(irow) = xpos; 
        stats_to_plot.bar_color(irow,1:3) = colororder(index_in_group,:); 
    end
end
[~,roworder] = sort(stats_to_plot.xpos); 
stats_to_plot = stats_to_plot(roworder,:); % set data in order of plot position
 
%% plotting
fig_module_ratios = figure;
bg = bar(stats_to_plot.xpos, stats_to_plot.mean_patchInterpatchRatio_projfilt);
hold on
eb = errorbar(stats_to_plot.xpos, stats_to_plot.mean_patchInterpatchRatio_projfilt, stats_to_plot.sem_patchInterpatchRatio_projfilt,'LineStyle','none');
hold off
bg.BaseValue = basevalue; 
plot_formatting();
bg.FaceColor = 'flat'; %%% allow for individual bar colors
bg.CData = stats_to_plot.bar_color; % color bars by source area
set(fig_module_ratios,'Renderer', 'painters', 'Position', [10 10 fig_width fig_height]) % set figure dimensions
if save_module_ratios 
    print(fig_module_ratios,savename_module_ratios, '-dtiffn','-r300') %%% save image as file
end




