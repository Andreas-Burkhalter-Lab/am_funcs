%%%%%%%%%%% perform anova and correlation tests on gcamp rois
%%%% for each table row, at least one criteria_var from each criterai var list must be true for a table row to be analyzed
%%%% ops.varname == the test statistic
%%%% ops.meets_criteria == roi table rows to be included in the analysis
%%%% ops.sorting_var == variable used to split rois into groups
%%%%%%%%%% updated on thermaltake 2021-1-25


clear pars labels_to_test values_to_test sorted_vals
clear ops
ops.meets_criteria = true(height(roitable),1); % reset meets_criteria
roitable.rf_area_sqdeg = pi * roitable.rf_size_yx_deg(:,1)/2 .* roitable.rf_size_yx_deg(:,2)/2; % sq deg; estimate as ellipse
deg_per_pix = 0.0608; %%% average taken from rf mapping
screenwidth_pix = 1920; % stats for computing RF locations
screenheight_pix = 1080;  % stats for computing RF locations
ops.newfig = field_default(ops,'newfig',1); % create a  new figure for each plot
ops.close_all = field_default(ops,'close_all',1); 

%% plotting options
ops.plot_bargraph =1; % options below
ops.plot_violin = 0;
    ops.overlay_violin_with_boxplot = 1; 
    ops.violin_median_color = [];
%     ops.violin_median_color = 'r';
    ops.violin_mean_color = 'k'; 
%     ops.violin_face_color= [1 0.5 0]; 
        ops.violin_face_color= [0.5 0.5 0.5]; 
    ops.violin_edge_color = [0 0 0 ];
    ops.violin_linewidth = 1.5;
    ops.violin_xvals = [1 1.7]; 
    ops.violin_xlim = [0.6 2.0]; 
    ops.violin_xaxis_visible = 0; 
    ops.violin_axis_thickness = 2.5;
ops.show_boxplot = 0;
    ops.boxplot_linewidth = 0.5; 
    ops.boxplot_show_notch = 1; 
    ops.boxplot_median_line = 1; % show line; if false, show median target
%     ops.boxplot_outlier_symbol = '*'; 
        ops.boxplot_outlier_symbol = ''; 
    ops.boxplot_border_color = [0.2 0.2 0.2];
ops.show_scatter = 0;
ops.show_histograms = 0; 
ops.anova_figures = 0;
ops.scatter_sessions = 0; 
    ops.scatter_session_rowcol = [6 6]; 

ops.show_refline = 0; % show line through y=0
    ops.refline_style = '--';
    ops.refline_color = [0.3 0.3 0.3]; 
    ops.refline_width = 1.2; 
    
    %% bargraph plot options
% ops.ebar_thickness = 1.5; 
    ops.ebar_thickness = 2;   % fig s1 quants
% ops.ebar_cap_size = 3.8;  % 
    ops.ebar_cap_size = 4.5;   % fig s1 quants
ops.ebar_color = [0 0 0]; 
ops.font_name = 'Arial'; 
% ops.axis_thickness = 3.5; %% does not apply to violin
    ops.axis_thickness = 4; % does not apply to violin; fig s1 quants
ops.bar_quant_colors = repmat(linspace(0.08,0.86,6)',1,3);
    ops.quant_bar_width = 0.9; %%% bar width for quantile plotting only
ops.bar_one_color = [0.5 0.5 0.5]; 
ops.bar_quant_border_color = [1 1 1]; 
ops.bargraph_line_width = 3;
ops.axis_font_bold = 1; % font bold or not bold
ops.axis_font_size = 19;


%     ops.ylim = [0.0 0.16]; 
ops.fig_outerposition = [0.45 0.15 0.22 0.60]; % [left bottom width height]... fig 2 pakanloc by module
%     ops.fig_outerposition = [0.05 0.05 0.42 0.65]; % [left bottom width height]... fig 2 orient_si, rf_sgnf
%     ops.fig_outerposition = [0.05 0.05 0.37 0.65]; % [left bottom width height]... fig 2 orient_si, ~rf_sgnf
%     ops.fig_outerposition = [0.05 0.05 0.35 0.6]; % [left bottom width height]... fig s1 quants
%     ops.fig_outerposition = [0.05 0.05 0.15 0.49]; % [left bottom width height]... fig 2B violin
    
ops.histogram_nbins = 50; 
% ops.histogram_normalization = 'probability'; % histogram y values are fraction of total rather than abs val
    ops.histogram_normalization = 'cdf'; % cumulative density
%     ops.histogram_normalization = 'counts'; % abs value of counts
ops.histogram_DisplayStyle = 'stairs';
%     ops.histogram_DisplayStyle = 'bar';



%% variables to test
% ops.varname = 'sf_pref';
% ops.varname = 'sf_hwhm';
% ops.varname = 'sf_di';
%     ops.ylim = [0.7 0.9]; 
% ops.varname = 'sf_sgnf';
% ops.varname = 'sf_anovap';
% ops.meets_criteria = ops.meets_criteria & roitable.sf_anovap < 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable.sf_anovap > 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable.sf_anovap < 0.1; 

% ops.varname = 'tf_pref';
% ops.varname = 'tf_hwhm';
% ops.varname = 'tf_di';
%     ops.ylim = [0.7 0.9]; 
% ops.varname = 'tf_sgnf';
% ops.varname = 'tf_anovap';
% ops.meets_criteria = ops.meets_criteria & roitable.tf_anovap < 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable.tf_anovap > 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable.tf_anovap < 0.1; 

% ops.varname = 'orient_hwhm';
ops.varname = 'orient_si';
%     ops.ylim = [0.35 0.55]; %%% for orient si by quant, rf_sgnf
%     ops.ylim = [0.39 0.49]; %%% for orient si by quant, ~rf_sgnf
% ops.varname = 'orient_di';
%     ops.ylim = [0.6 0.75]; 
% ops.varname = 'orient_sgnf';
% ops.varname = 'orient_pref1';
% ops.varname = 'orient_pref2';
% ops.varname = 'orient_anovap';
% ops.meets_criteria = ops.meets_criteria & roitable.orient_anovap < 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable.orient_anovap > 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable.orient_anovap < 0.1; 

ops.close_all = 1; 
if ops.close_all; close all; end

%%%%% RF variables
% ops.varname = 'rf_width_deg'; roitable.rf_width_deg = roitable.rf_size_yx_deg(:,2); 
% ops.varname = 'rf_height_deg'; roitable.rf_height_deg = roitable.rf_size_yx_deg(:,1);
% ops.varname = 'rf_area_sqdeg'; 
% ops.varname = 'rf_offcenter_pix_y'; % distance in pix of rf from screen center
% ops.varname = 'rf_offcenter_pix_x'; % distance in pix of rf from screen center
% ops.varname = 'rf_center_pix_x';
% ops.varname = 'rf_center_pix_y';
%     ops.varname = 'rf_center_deg_y'; roitable.rf_center_deg_y = deg_per_pix*[roitable.rf_center_pix_y-screenheight_pix/2];
% ops.varname = 'stim_overlap_frac';
% ops.varname = 'rf_sgnf';
% ops.varname = 'anovap_rf';
% ops.varname = 'stim_on_rf';
% ops.meets_criteria = ops.meets_criteria & roitable.rf_area_sqdeg > 270; % median ~210 among rf-tuned cells
%     ops.meets_criteria = ops.meets_criteria & roitable.rf_area_sqdeg < 270; % median ~210 among rf-tuned cells
% ops.meets_criteria = ops.meets_criteria & roitable.rf_size_yx_deg(:,1) < 21; %%% RF height requirement; median ~18
% ops.meets_criteria = ops.meets_criteria & roitable.rf_size_yx_deg(:,2) > 17; %%% RF width requirement; median ~19
%     ops.meets_criteria = ops.meets_criteria & roitable.rf_size_yx_deg(:,2) < 17; %%% RF width requirement; median ~19

%%%% locomotion and pupil variables
% ops.varname = 'locm_cor_r';
%     ops.varname = 'locm_cor_r_abs'; roitable.locm_cor_r_abs = abs(roitable.locm_cor_r);
% ops.varname = 'locm_cor_p';
% ops.varname = 'pakan_loc_abs'; roitable.pakan_loc_abs = abs(roitable.pakan_loc_index);
% ops.varname = 'pakan_loc_index';
%     ops.violin_bandwidth = 0.015; 
%     ops.ylim = [-0.5 0.75]; % fig 2B left
%     ops.ylim = [-0.7 0.75]; % fig 2B right
% ops.varname = 'pakan_loc_sgnf'; 
% ops.varname = 'pakan_loc_positive'; roitable.pakan_loc_positive = double(roitable.pakan_loc_index>0);
% ops.varname = 'pupilwidth_cor_r'; 
% ops.varname = 'pupilwidth_cor_p';
%
% ops.meets_criteria = roitable.pakan_loc_index > 0; 
% ops.meets_criteria = ops.meets_criteria & roitable.pakan_loc_sgnf; %%% significant difference in stimulus responses on locomotion trials vs. stationary trials
% ops.meets_criteria = ops.meets_criteria & roitable.locm_cor_p < 0.05; % locomotion-tuned in darkness
%   ops.meets_criteria = ops.meets_criteria & roitable.locm_cor_p > 0.05; % not locomotion-tuned in darkness
% ops.meets_criteria = ops.meets_criteria & roitable.locm_cor_r > 0.00; % positively locomotion-tuned in darkness
%   ops.meets_criteria = ops.meets_criteria & roitable.locm_cor_r < 0.00; % negatively locomotion-tuned in darkness
% ops.meets_criteria = ops.meets_criteria & roitable{:,'pupilwidth_cor_sgnf'};
 
%% criteria variables
%%%% criteria variables related to stim response
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} | roitable{:,'rspv'} ]; % if roi is responsive or tuned for sf/tf/orient
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} | roitable{:,'rf_sgnf'} ]; % if roi is tuned for sf/tf/orient/rf
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} ]; 
%   ops.meets_criteria = ops.meets_criteria & [ roitable.sf_anovap<0.1 | roitable.tf_anovap<0.1 | roitable.orient_anovap<0.1 ]; % weakly tuned for at least one of sf/tf/orient
ops.meets_criteria = ops.meets_criteria & roitable.anovap_rf < 0.05; % if roi has detectable receptive field
%   ops.meets_criteria = ops.meets_criteria & ~roitable{:,'rf_sgnf'}; % if roi does not have detectable receptive field
% ops.meets_criteria = ops.meets_criteria & roitable{:,'rspv'};

% %%%% criteria variables related to receptive field overlap with stim
ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac > 0.00; % stim on RF
% ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac == 0; % stim off RF
% ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac == 0 | ops.meets_criteria & ~roitable{:,'rf_sgnf'}; %% stim off RF or not RF-tuned
% ops.meets_criteria = ops.meets_criteria & roitable.rf_overlap_frac_aprx > 0;
% ops.meets_criteria = ops.meets_criteria & ~any(isnan(roitable.rf_size_yx_deg),2); % rf was not cut off by the screen and rf width/height were calculated

%%%% specify sorting variable
% ops.sorting_var = 'inpatch'; 
%     ops.xlim = [0.3 2.7]; % does not apply to violin
ops.sorting_var = 'quantile';
%     ops.xlim = [0.3 6.7]; % does not apply to violin
% ops.sorting_var = 'plane';
% ops.sorting_var = 'sub';
% ops.sorting_var = 'day';
% ops.sorting_var = 'sf_sgnf';    % sf-tuned vs. non-sf-tuned cells

%%%%% specify subject, layer, plane, day, quantile
% ops.meets_criteria = ops.meets_criteria & any(roitable.plane == [2], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.layer == [23], 2); 
% ops.meets_criteria = ops.meets_criteria & strcmp(roitable.day, filetable.day{32});   %%%% only analyze particular recording sessions
%     ops.meets_criteria = ops.meets_criteria & ~strcmp(roitable.day, filetable.day{25});   %%%% only analyze particular recording sessions
%     ops.meets_criteria = ops.meets_criteria & ~strcmp(roitable.day, filetable.day{26});   %%%% only analyze particular recording sessions
%     ops.meets_criteria = ops.meets_criteria & ~strcmp(roitable.day, filetable.day{28});   %%%% only analyze particular recording sessions
%     ops.meets_criteria = ops.meets_criteria & ~strcmp(roitable.day, filetable.day{32});   %%%% only analyze particular recording sessions
% ops.meets_criteria = ops.meets_criteria & [strcmp(roitable.day, filetable.day{1}) | strcmp(roitable.day, filetable.day{2})]; % analyze pair of recording sessions
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [18164], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [18165], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [19002], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [19042], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [19048], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [20003], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [20004], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [20015], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub ~= [19048], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub ~= [19089], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub ~= [20004], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.quantile == [5 6],2); %%% only analyze specific quantiles

%%% don't count cells that fall in areas excluded from patchiness analysis, which are coded as quantile==0
ops.meets_criteria = ops.meets_criteria & roitable.quantile>0;

%%% exclude rois that have value NaN for the variable of interest
ops.meets_criteria = ops.meets_criteria & ~isnan(roitable{:,ops.varname}); 



%% analysis
% get values for each group so standard error can be calculated
labelnames = unique(roitable{ops.meets_criteria,ops.sorting_var});
allvals = roitable{ops.meets_criteria, ops.varname}; %%% all variables of varname for ROIs that meet criteria
nans = NaN(size(labelnames));
sorted_vals = table(labelnames, cell(size(labelnames)), nans, nans, nans, nans, 'VariableNames', {ops.sorting_var, ops.varname,'nrois','mean','median','sem'});
for ilabel = 1:length(labelnames)
    if iscell(roitable{1,ops.sorting_var}) && isstr(roitable{1,ops.sorting_var}{:}) % if sorting var is string inside cell
        thesevals = roitable{ops.meets_criteria & strcmp(roitable{:,ops.sorting_var},labelnames{ilabel}), ops.varname};
    else %% else assume sorting var is numeric
        thesevals = roitable{ops.meets_criteria & roitable{:,ops.sorting_var}==labelnames(ilabel),ops.varname};
    end
    sorted_vals{ilabel,ops.varname} = {thesevals};
    sorted_vals.nrois(ilabel) = length(thesevals);
    sorted_vals.mean(ilabel) = mean(thesevals);
    sorted_vals.sem(ilabel) = std(thesevals) ./ sqrt(length(thesevals));
    sorted_vals.median(ilabel) = nanmedian(thesevals);
end

% do anova
values_to_test = roitable{ops.meets_criteria,ops.varname}; 
if iscell(roitable{1,ops.sorting_var}) && isstr(roitable{1,ops.sorting_var}{:}) % if sorting var is string inside cell
    labels_to_test = roitable{ops.meets_criteria,ops.sorting_var};
else %% else assume sorting var is numeric
    labels_to_test = double(roitable{ops.meets_criteria,ops.sorting_var});
    [r_corr p_corr] = corrcoef(values_to_test, labels_to_test); % only do correlation analysis if input is numeric
end
if ops.anova_figures
    anova_figures = 'on';
else
    anova_figures = 'off';
end
[anovap, anovatab, anovastats] = anova1(values_to_test, labels_to_test, anova_figures);
anovastats
fprintf(['p_corr = ', num2str(p_corr(2,1)), '... r_corr = ', num2str(r_corr(2,1)), '... anova_p = ', num2str(anovap), '\n'])

%% plotting
%% bargraph
if ops.plot_bargraph
    if ops.newfig; hfig = figure('units','normalized','outerposition',ops.fig_outerposition); end
    hax = gca; 

    % shade the bars individually if analyzing quantile
    if strcmp(ops.sorting_var,'quantile')
        nquants = nnz(unique(roitable.quantile));
        xvals = 1:nquants; 
        barhandle = bar(xvals, diag(sorted_vals.mean), 'stacked');
        for ibar = 1:nquants % for each quantile to be plotted
            set(barhandle(ibar),'facecolor', ops.bar_quant_colors(ibar,:))
            set(barhandle(ibar),'EdgeColor', ops.bar_quant_border_color);
            barhandle(ibar).LineWidth = ops.bargraph_line_width;
            barhandle(ibar).BarWidth = ops.quant_bar_width;
        end
    else
        barhandle = bar(sorted_vals.mean); 
        barhandle.FaceColor = ops.bar_one_color; 
        barhandle.LineWidth = ops.bargraph_line_width; 
        hax.XTick = []; 
    end
    
    hold on
    ebar = errorbar(1:height(sorted_vals),sorted_vals.mean,sorted_vals.sem,'.','Color',ops.ebar_color);
    hold off
    hax.FontName = ops.font_name;
    hax.FontSize = ops.axis_font_size;
    hax.LineWidth = ops.axis_thickness;
    ebar.LineWidth = ops.ebar_thickness; % set errorbar line width
    ebar.Color = ops.ebar_color; % set errorbar line color
    ebar.CapSize = ops.ebar_cap_size; 
    if isfield(ops,'ylim'); ylim(ops.ylim); end
    if isfield(ops,'xlim'); xlim(ops.xlim); end
    if ops.axis_font_bold; set(gca,'FontWeight','bold'); end
    set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
    set(gca,'Box','off')
end

%% violin and boxplot
if ops.plot_violin
    if ops.newfig; violin_fig = figure('units','normalized','outerposition',ops.fig_outerposition); end
    ops.violin_xvals = field_default(ops,'violin_xvals',[]); 
    if isfield(ops, 'violin_bandwidth') % if bandwidth was specified
        [hviolin,hleg,group_means,group_medians,bw] = violin(sorted_vals{:,ops.varname}','bw',ops.violin_bandwidth,'mc',ops.violin_mean_color,'medc',ops.violin_median_color,...
            'facecolor',ops.violin_face_color,'edgecolor',ops.violin_edge_color,'x',ops.violin_xvals);
    elseif ~isfield(ops, 'violin_bandwidth') % use default bandwidth if not specified
        [hviolin,hleg,group_means,group_medians,bw] = violin(sorted_vals{:,ops.varname}','mc',ops.violin_mean_color,'medc',ops.violin_median_color,...
            'facecolor',ops.violin_face_color,'edgecolor',ops.violin_edge_color,'x',ops.violin_xvals);      
    end
    if ops.overlay_violin_with_boxplot
        hold on
    end
    hax = gca;
    hax.FontName = ops.font_name;
    hax.FontSize = ops.axis_font_size;
    if ops.axis_font_bold; set(gca,'FontWeight','bold'); end
    hax.LineWidth = ops.violin_axis_thickness; 
    for iplot = 1:length(hviolin); hviolin(iplot).LineWidth = ops.violin_linewidth; end
    set(gca,'xtick',[])
    if ~ops.violin_xaxis_visible; hax.XAxis.Visible = 'off'; end
    hleg.Visible = 'off'; % turn legend off
    if isfield(ops,'ylim'); ylim(ops.ylim); end
    if isfield(ops,'violin_xlim'); xlim(ops.violin_xlim); end
    set(gca,'Box','off')
    set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
end
if ops.show_boxplot
    if ops.newfig && [~ops.overlay_violin_with_boxplot || ~ops.plot_violin]
        figure
    end
    if ops.boxplot_show_notch; boxplot_show_notch = 'on'; else boxplot_show_notch = 'off'; end
    if ops.boxplot_median_line; boxplot_MedianStyle = 'line'; else boxplot_MedianStyle = 'target'; end
    bp = boxplot(values_to_test,labels_to_test,'notch',boxplot_show_notch,'MedianStyle',boxplot_MedianStyle,'Symbol',ops.boxplot_outlier_symbol,...
        'Colors',ops.boxplot_border_color);
    set(bp,{'linew'},{ops.boxplot_linewidth})
    if isfield(ops,'ylim'); ylim(ops.ylim); end
    hold off
end


if ops.show_histograms
    if ops.newfig; figure; end
    valrange = [min(allvals), max(allvals)]; 
    hist_edges = linspace(valrange(1), valrange(2), ops.histogram_nbins+1); 
    for ii = 1:height(sorted_vals)
        sorted_vals.hgram{ii} = histogram(sorted_vals{ii,ops.varname}{:},hist_edges,'Normalization',ops.histogram_normalization,'DisplayStyle',ops.histogram_DisplayStyle);
        hold on
    end
    hold off
end
if ops.show_scatter
    if ops.newfig; figure; end
    scplot = scatter(labels_to_test,allvals);
end
if ops.scatter_sessions
    sessions = table;
    if ops.newfig; ff = figure('units','normalized','outerposition',[0 0 1 1]); end
    yvalrange = range(allvals);
    y_lims = [min(allvals)-0.1*yvalrange, max(allvals)+0.1*yvalrange];
    for i = 1:height(filetable)
        subplot(ops.scatter_session_rowcol(1),ops.scatter_session_rowcol(2),i)
        rr = ops.meets_criteria & strcmp(roitable.day, filetable.day{i}) & [roitable.plane == filetable.plane(i)]; % match day and plane
        warning('off','all');       sessions.day{i} = filetable.day{i};     warning('on','all')
        sessions.vals{i} = roitable{rr, ops.varname}; 
        sessions.labels{i} = double(roitable{rr,ops.sorting_var});
        [sess_anovap, sess_anovatab, sess_anovastats] = anova1(sessions.vals{i}, sessions.labels{i}, 'off');
        sessions.anovap{i} = sess_anovap;
        [sess_r_corr, sess_p_corr] = corrcoef(sessions.vals{i}, sessions.labels{i});
        sessions.r_corr{i} = sess_r_corr(length(sess_r_corr),1);
        sessions.p_corr{i} = sess_p_corr(length(sess_r_corr),1);
        scplot = scatter(sessions.labels{i},sessions.vals{i});
        title([filetable.day{i}, ', pcor=', num2str(sessions.p_corr{i})])
        ylim(y_lims)
    end
end

% reference line at y=0
if ops.show_refline && ... % show line through y=0
        any([ops.plot_bargraph,ops.plot_violin,ops.show_boxplot,ops.show_scatter,ops.show_histograms,ops.anova_figures,ops.scatter_sessions])  
    hrefline = refline(0,0);
    hrefline.LineStyle = ops.refline_style; 
    hrefline.Color = ops.refline_color; 
    hrefline.LineWidth = ops.refline_width;
end
