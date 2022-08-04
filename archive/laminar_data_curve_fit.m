%%%% plot cells by distance from injection site and section, sorted by categories patch, interpatch, border
%%%%% plotting for weiqing data
%
% input excel_file should have labeled rows for section, distance, module
% excel_rows = rows from excel table to include (default == all rows)
%
% module labels for patch, interpatch, border are p, ip, o
%
%
% updated 2020/1/14 

function data_out = laminar_data_curve_fit(excel_file, excel_rows)

line_width = 1.5;
xlimits = [0 12]; 
ylimits = [0 800]; 
interpatch_color = [0 0 0];
patch_color = [1 0 0];
border_color = [0.7 0.7 0.7];
x_axis_label = 'Section';
y_axis_label = 'Distance from injection center (um)';
module_labels = {'ip','p','o'}; % determines order of labels in legend; must match labels in excel spreadsheet
% curve_fit_model = 'smoothingspline'; % gives smoother appearance than cubic interpolation
curve_fit_model = 'cubicinterp'; % exaggerates variations between neighboring sections, causes problems with missing data

data_out = readtable(excel_file);
excel_rows = vardefault('excel_rows',[1:height(data_out)] + 1); % include all excel rows if not specified
table_rows_to_keep = excel_rows - 1; % adjust for top row with column labels in excel document
data_out = data_out(table_rows_to_keep,:);
% % % % % % % % module_labels = unique(data_out.module); 
nlabels = length(module_labels); 
section_labels = unique(data_out.section);
n_section_labels = length(section_labels);
legend_label_list = {}; % for figure legend

figure 
hold on
for ilabel = 1:nlabels
    thislabel = module_labels{ilabel};
    switch thislabel
        case 'ip'
            label_color = interpatch_color;
            legend_label_list = [legend_label_list; 'Interpatch'];
        case 'p'
            label_color = patch_color;
            legend_label_list = [legend_label_list; 'Patch'];
        case 'o'
            label_color = border_color;
            legend_label_list = [legend_label_list; 'Border'];
    end
    row_inds = strcmp(thislabel, data_out.module);
    this_x_dat = data_out.section(row_inds);
    this_y_dat = data_out.distance(row_inds);
    thisdat = data_out(row_inds,:);
    scatter(this_x_dat, this_y_dat,'MarkerFaceColor',label_color,'MarkerEdgeColor',label_color)
    module_stats = grpstats(thisdat(:,{'section','distance'}),'section');
    fitobject = fit(module_stats.section, module_stats.mean_distance, curve_fit_model);
    plotobj = plot(fitobject); 
    plotobj.Color = label_color;
    plotobj.LineWidth = line_width;
    plotobj.HandleVisibility = 'off'; % make curve invisible from legend
    xlabel(x_axis_label)
    ylabel(y_axis_label)
end

legobj = legend(legend_label_list);
xlim(xlimits)
ylim(ylimits)