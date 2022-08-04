%%%% plot a bar graph with error bars of comparison image intensity in M2 quants
%  before running this function, run physpaper_anatomy_quants_analysis.m
%
% updated 2020-9-1 on thermaltake

close all

plot_group_bars = 1; %% plot group means and SEM
plot_individual_bars = 1; %% plot each individual subject in a subplot
    subplot_rowcol = [3 2]; 
plot_shuffled_data = 0; %%% plot shuffled instead of original data
    
bar_colors = repmat(linspace(0.08,0.86,6)',1,3);
bar_border_width = 1; 
bar_border_color = 'none';
ebar_color = 'k'; 
ebar_linewidth = 2.5; 
axes_font_size = 22;
axes_line_width = 3; 
xlimits = [0.2, 6.8]; 
x_tick_length = [0 0];
fig_width = 520;
fig_height = 500; 
show_title = 0; 
save_quants_bar_fig = 0; %%% save output as .eps image
    savename_quants_bar_fig = ['F:\thesis\figs\fig 4 -- quants_bargraph_fig_', pathway_name];

%%%% load the output of physpaper_anatomy_quants_analysis.m
% load('F:\thesis\som_m2\shuffled_patchiness_data')
