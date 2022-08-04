%%%% plot a bar graph with error bars of amyg projections to POR quants
% 
%
% 
% updated 2021-1-07 on thermaltake

close all

% % % % % % bar_color = [ 0.6 0.6 0.6]; 
% bar_colors = repmat(linspace(0.08,0.86,6)',1,3);
    bar_colors = repmat(linspace(1,1,6)',1,3);
bar_line_width = 2; 
% bar_border_color = 'none';
    bar_border_color = [0 0 0]; 
ebar_color = 'k'; 
ebar_linewidth = 2; 
axis_font_size = 16;
xlimits = [0.2, 6.8]; 
x_tick_length = [0 0];
fig_width = 520;
fig_height = 500; 
add_scatter_points = 1; % display all individual data points overlayed with bar graphs; remove bar fill shading
    max_x_deviation = 0.3; % maximum x spacing jitter
    scatter_color = [0.5 0.5 0.5]; 
    scatter_marker_size = 70;
save_amyg_plot = 0; 
    savename_amyg_fig = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\fig_amyg_to_por_quants';

%%%% load the output of amyg_proj_quants_analysis.m
% load('C:\Users\Burkhalter Lab\Documents\anatomy_paper\amyg\amyg_proj_quants_analysis')

% output correlation coefficient and p value for quant number vs egfp intensity
egfp_by_quant = filetable.proj_intens_by_quant_normed; 
nquants = size(egfp_by_quant,2);
quantvals = repmat(1:nquants,height(filetable),1);
[r, p] = corrcoef(quantvals, egfp_by_quant)

% make figure
amyg_fig = figure; 
% % % % % % hbar = bar(mean(filetable.proj_intens_by_quant_normed), 'FaceColor',bar_color, 'LineWidth',bar_border_width);
data_to_plot = filetable.proj_intens_by_quant_normed;
data_to_plot = max(data_to_plot,zeros(size(data_to_plot))); % zero floor
mean_proj_vals = mean(data_to_plot);
xvals = 1:6;
hbar = bar(xvals, diag(mean_proj_vals), 'stacked');
hold on

% show data from individual cases
if add_scatter_points
    mgrid = meshgrid(1:size(data_to_plot,2),1:1:size(data_to_plot,1)); % pre-jittered indices indices
    noisegrid = [2*max_x_deviation * rand(size(data_to_plot)) - max_x_deviation]; % x jitter falues
    x_indices = mgrid + noisegrid; % jittered indices
    scatplot = scatter(x_indices(:), data_to_plot(:)); 
    scatplot.SizeData = scatter_marker_size; 
    scatplot.MarkerEdgeColor = scatter_color;
    scatplot.MarkerFaceColor = scatter_color;
    bar_face_color = [1 1 1]; % white bar face
end


sem = std(data_to_plot) ./ sqrt(nsubs);
errorbar(mean_proj_vals,sem,'.','Color',ebar_color,'LineWidth',ebar_linewidth)
plot_formatting; % format the plot

set(gca,'XTick',[1:6])
set(gca,'XTickLabels',{'1','2','3','4','5','6'})  %%%% restor x tick labels
hax = gca;
hax.XAxis.TickLength = x_tick_length;
% % % % % % hbar.EdgeColor = bar_border_color;
for ibar = 1:6
    set(hbar(ibar),'facecolor',bar_colors(ibar,:))
    set(hbar(ibar),'EdgeColor',bar_border_color);
    set(hbar(ibar),'LineWidth',bar_line_width);
end
set(amyg_fig,'Renderer', 'painters', 'Position', [10 10 fig_width fig_height]) % set figure dimensions

if save_amyg_plot
%     print(amyg_fig,savename_amyg_fig, '-dtiffn','-r300') %%% save image as file
    saveas(amyg_fig,savename_amyg_fig, 'svg') %%% save image as file
end