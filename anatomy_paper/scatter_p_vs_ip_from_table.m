%%%% scatter patch vs interpatch ratios, with x = interpatch intensity, y = patch intensity
% input table has variables patchIntens, interpatchIntens, subject
% each subject but plotted in different color
%%%% updated 2019/9/23

% % % % % % % function scatter_p_vs_ip_from_table(input_table)

% load('C:\Users\Burkhalter Lab\Documents\anatomy_paper\amyg\amyg_P_vs_IP_perimeter_test_data')
input_table = amyg;
single_marker_color = 1;
    marker_color = [0 0 1];
marker_size = 100;
plotmin = 0.4; % min x and y values
add_legend = 0;     
axis_font_size = 20; 

xlabel_text = 'Patch EGFP optical density';
ylabel_text = 'Interpatch EGFP optical density';

subject_list = unique(input_table.subject);
nsubjects = length(subject_list);

ff = figure;
for isub = 1:nsubjects
    this_sub = subject_list(isub);
    tablerows = input_table.subject == this_sub;
    ss = scatter(input_table.patchIntens(tablerows), input_table.interpatchIntens(tablerows),'o','filled','SizeData',marker_size);
    set(gca,'FontSize',axis_font_size)
    if single_marker_color
        ss.MarkerFaceColor = marker_color;
        ss.MarkerEdgeColor = marker_color;
    end
    hold on
end

xymax = max([ylim xlim]); 
xlim([plotmin xymax]); 
ylim([plotmin xymax]);
xx = .001:.001:max([xlim ylim]);
hold on;
plot(xx,xx,'k--')
if add_legend
    legend(cellstr(num2str(subject_list)))
end
hold off

xlabel(xlabel_text)
ylabel(ylabel_text)