%%%%%%%%%%% perform anova and correlation tests on gcamp rois
%%%% for each table row, at least one criteria_var from each criterai var list must be true for a table row to be analyzed
%%%% ops.varname == the test statistic
%%%% ops.meets_criteria == roi table rows to be included in the analysis
%%%% ops.sorting_var == variable used to split rois into groups
%%%%%%%%%% updated on thermaltake 2020-01-03


% load('F:\analyses\cor5-21')
axes_line_width = 1.7;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
axis_font_size = 10; 
curve_LineWidth = 2; % width of line of plotted data


close all
ops.meets_criteria = true(height(roitable),1); % reset meets_criteria
clear labels_to_test values_to_test sorted_vals

%% plottomg options
ops.plot_bargraph = 0; 
ops.show_boxplot = 0;
ops.show_histograms = 1; 
ops.anova_figures = 0;

ops.boxplot_linewidth = 1; 
ops.histogram_nbins = 50; 
% ops.histogram_normalization = 'probability'; % histogram y values are fraction of total rather than abs val
    ops.histogram_normalization = 'cdf'; % cumulative density
%     ops.histogram_normalization = 'counts'; % abs value of counts
ops.histogram_DisplayStyle = 'stairs';
%     ops.histogram_DisplayStyle = 'bar';



%% select data and inclusion criteria
% ops.varname = 'sf_pref';
ops.varname = 'sf_hwhm';
% ops.varname = 'sf_sgnf';
% ops.varname = 'sf_anovap';
ops.meets_criteria = ops.meets_criteria & roitable{:,'sf_sgnf'};

% ops.varname = 'tf_pref';
% ops.varname = 'tf_hwhm';
% ops.varname = 'tf_sgnf';
% ops.varname = 'tf_anovap';
% ops.meets_criteria = ops.meets_criteria & roitable{:,'tf_sgnf'};

% ops.varname = 'orient_hwhm';
% ops.varname = 'orient_sgnf';
% ops.varname = 'orient_pref1';
% ops.varname = 'orient_pref2';
% ops.varname = 'orient_anovap';
% ops.meets_criteria = ops.meets_criteria & roitable{:,'orient_sgnf'};

%%%% miscellaneous variables of interest
% ops.varname = 'locm_cor_r';
% ops.varname = 'pakan_loc_index';
% ops.varname = 'pakan_loc_positive'; roitable.pakan_loc_positive = double(roitable.pakan_loc_index>0);
% ops.varname = 'rf_width_deg'; roitable.rf_width_deg = roitable.rf_size_yx_deg(:,1); 
% ops.varname = 'rf_height_deg'; roitable.rf_height_deg = roitable.rf_size_yx_deg(:,2);
% ops.varname = 'stim_overlap_frac';
% ops.varname = 'rf_sgnf';
% ops.varname = 'stim_on_rf';
% ops.varname = 'pupilwidth_cor_r'; 
% ops.varname = 'pupilwidth_cor_p';
% ops.meets_criteria = ops.meets_criteria & roitable.locm_cor_p > 0.05; 
% ops.meets_criteria = ops.meets_criteria & roitable{:,'pupilwidth_cor_sgnf'};
 
%%%% criteria variables related to stim response
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} | roitable{:,'rspv'} ]; % if roi is responsive or tuned for sf/tf/orient
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} | roitable{:,'rf_sgnf'} ]; % if roi is tuned for sf/tf/orient/rf
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} ]; 
% ops.meets_criteria = ops.meets_criteria & roitable{:,'rf_sgnf'}; % if roi has detectable receptive field
% ops.meets_criteria = ops.meets_criteria & ~roitable{:,'rf_sgnf'}; % if roi does not have detectable receptive field
% ops.meets_criteria = ops.meets_criteria & roitable{:,'rspv'};

%%%% criteria variables related to receptive field
ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac > 0.15; 
% ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac == 0; 
% ops.meets_criteria = ops.meets_criteria & roitable.rf_overlap_frac_aprx > 0;

%%%%% specify subject, layer, plane, day, quantile
% ops.meets_criteria = ops.meets_criteria & any(roitable.plane == [4 3 2 1], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.layer == [23], 2); 
% ops.meets_criteria = ops.meets_criteria & strcmp(roitable.day, filetable.day{1});   %%%% only analyze particular recording sessions
% ops.meets_criteria = ops.meets_criteria & [strcmp(roitable.day, filetable.day{1}) | strcmp(roitable.day, filetable.day{2})]; % analyze pair of recording sessions
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [18164], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [18165], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [19002], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.sub == [19042], 2); 
% ops.meets_criteria = ops.meets_criteria & any(roitable.quantile == [1 2 5 6],2); %%% only analyze specific quantiles

%%%% specify sorting variable
ops.sorting_var = 'inpatch'; 
% ops.sorting_var = 'quantile';
% ops.sorting_var = 'plane';
% ops.sorting_var = 'sub';
% ops.sorting_var = 'day';
% ops.sorting_var = 'sf_sgnf';    % sf-tuned vs. non-sf-tuned cells

%%% don't count cells that fall in areas excluded from patchiness analysis, which are coded as quantile==0
ops.meets_criteria = ops.meets_criteria & roitable.quantile>0;

%%% exclude rois that have value NaN for the variable of interest
ops.meets_criteria = ops.meets_criteria & ~isnan(roitable{:,ops.varname}); 

%% analysis
% get values for each group so standard error can be calculated
labelnames = unique(roitable{ops.meets_criteria,ops.sorting_var});
allvals = roitable{ops.meets_criteria, ops.varname}; %%% all variables of varname that meet criteria
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
    [r_corr p_corr] = corrcoef(values_to_test, labels_to_test) % only do correlation analysis if input is numeric
end
if ops.anova_figures
    anova_figures = 'on';
else
    anova_figures = 'off';
end
[anovap, anovatab, anovastats] = anova1(values_to_test, labels_to_test, anova_figures);
anovap
anovastats

%% plotting
if ops.plot_bargraph
    figure
    bar(sorted_vals.mean)
    hold on
    errorbar(1:height(sorted_vals),sorted_vals.mean,sorted_vals.sem,'.')
    hold off
end
if ops.show_boxplot
    figure;
    bp = boxplot(values_to_test,labels_to_test,'notch','on');
    set(bp,{'linew'},{ops.boxplot_linewidth})
    hold off
end
if ops.show_histograms
    histfig = figure;
    valrange = [min(allvals), max(allvals)]; 
    hist_edges = linspace(valrange(1), valrange(2), ops.histogram_nbins+1); 
    for ii = 1:height(sorted_vals)
        sorted_vals.hgram{ii} = histogram(sorted_vals{ii,ops.varname}{:},hist_edges,'Normalization',ops.histogram_normalization,'DisplayStyle',ops.histogram_DisplayStyle);
        sorted_vals.hgram{ii}.LineWidth = curve_LineWidth;
        hold on
    end
    hold off
end


x_incr = mean(diff(hist_edges));
xlim([hist_edges(1)-x_incr, hist_edges(end)])
ylim([0 1.05])
set(gca,'Box','off')
set(gca,'linewidth',axes_line_width)
set(gca,'FontSize',axis_font_size)
set(gca,'FontWeight',axes_numbers_bold)
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders

print(histfig,'sfhwhm_cdf.tif', '-dtiffn','-r300') %%% save image as file


