%%% format appearance of box plots
% updated 2020/3/1

%% set defaults
axes_line_width = vardefault('axes_line_width',2);
axis_font_size = vardefault('axis_font_size',13); 
show_notch = vardefault('show_notch','off');
y_log_or_linear = vardefault('y_log_or_linear','linear');
axes_numbers_bold = vardefault('axes_numbers_bold','bold'); %%% 'bold' or 'normal'

%%
if show_notch
    notch_toggle = 'on';
else
    notch_toggle = 'off';
end

line_objs = findobj(gcf, 'type', 'line', 'Tag', 'Median'); % for boxplot
if ~isempty(line_objs)  % if boxplot
    set(line_objs, 'Color', 'r');
    set(bp,{'linew'},{boxplot_linewidth})
end
if exist('bg','var') && isgraphics(bg) %%% for bar graph
    bg.FaceColor = bar_face_color; 
    bg.LineWidth = bar_line_width;
    eb.LineWidth = bar_line_width; % set errorbar line width
    eb.Color = error_line_color; % set errorbar line color
end
set(gca,'YScale',y_log_or_linear)
set(gca,'Box','off')
set(gca,'XColor',[0 0 0]);
set(gca,'YColor',[0 0 0]);
set(gca,'linewidth',axes_line_width)
set(gca,'FontSize',axis_font_size)
set(gca,'FontWeight',axes_numbers_bold)
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
set(gca,'XTickLabels',{'','','','','',''})  %%%% turn off x tick labels
set(gca,'XTick',[]) % turn off x ticks
set(gca,'FontName','Arial')

if exist('xlimits','var')
    xlim([xlimits])
end
if exist('ylimits','var')
    ylim([ylimits])
end
if exist('ytick','var')
    set(gca,'YTick',ytick); 
end