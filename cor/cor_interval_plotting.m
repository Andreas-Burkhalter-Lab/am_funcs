%%% plot correlations of pairs that fall in a particular distance interval
% this script should be called by cor_plotting_multi_groups.m
%
% updated 2022/3/8

% close all

clear vals_to_plot matchrows

% find pairs that fall in the specified distance interval
matchrows = cortable_out.umdist>intervalpars.umdist_minmax(1) & cortable_out.umdist<intervalpars.umdist_minmax(2); %%% apply distance interval
matchrows = matchrows & analysis_pars.include; %%% apply inclusion criteria that were specified in cor_plotting_multi_groups.m
vals_to_plot = cortable_out{matchrows, intervalpars.var_to_plot}; 
labels_to_plot = cortable_out.quantgroup_pair(matchrows); % labels = interpatch-interpatch (1), patch-interpatch (2), patch-patch (3)
labelnames = unique(labels_to_plot(~isnan(labels_to_plot)));
ngroups = nnz(unique(cortable_out.quantgroup_pair)); 

% sort cor values by quantgroup_pair
nans = NaN(size(labelnames));
sorted_vals = table(labelnames, cell(size(labelnames)), nans,       nans, nans, nans, 'VariableNames',...
            {'quantgroup_pair', intervalpars.var_to_plot,'npairs','mean','median','sem'});
for ilabel = 1:length(labelnames) % for each module pair type
    thesevals = cortable_out{matchrows & labelnames(ilabel)==cortable_out.quantgroup_pair, intervalpars.var_to_plot};
    sorted_vals{ilabel,intervalpars.var_to_plot} = {thesevals}; % cor vals from this module pair type
    sorted_vals.npairs(ilabel) = length(thesevals);
    sorted_vals.mean(ilabel) = mean(thesevals);
    sorted_vals.sem(ilabel) = std(thesevals) ./ sqrt(length(thesevals));
    sorted_vals.median(ilabel) = nanmedian(thesevals);
end

% plotting
if intervalpars.show_boxplot
    interval_fig = figure('units','normalized','outerposition',intervalpars.fig_outerposition,'WindowStyle',windowstyle); 
    boxplot(vals_to_plot, labels_to_plot, 'notch',intervalpars.notchstyle, 'symbol',intervalpars.outlier_symbol)
end

if intervalpars.show_bargraph
    interval_fig = figure('units','normalized','outerposition',intervalpars.fig_outerposition,'WindowStyle',windowstyle);
    barhandle = bar(1:ngroups, diag(sorted_vals.mean), 'stacked');
    hold on
    bargraph_ebar = errorbar(1:height(sorted_vals),sorted_vals.mean,sorted_vals.sem,'.');
    bargraph_ebar.LineWidth = intervalpars.barerror_line_width; % set errorbar line width
    bargraph_ebar.Color = intervalpars.error_line_color; % set errorbar line color
    for ibar = 1:ngroups
        set(barhandle(ibar),'facecolor',min([1 1 1], intervalpars.bar_colors{ibar}))
        set(barhandle(ibar),'EdgeColor',intervalpars.border_colors{ibar});
        barhandle(ibar).LineWidth = intervalpars.bargraph_line_width;
    end
    hold off
end

if pars.show_title; title([plot_titles{ivar}, titlestring]); end
hax = gca;
hax.FontSize = pars.axis_font_size;
hax.FontName = pars.font_name;
hax.LineWidth = intervalpars.axis_thickness; 
if pars.axis_font_bold; set(gca,'FontWeight','bold'); end
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
hax.Box = 'off'; 
hax.XLim = intervalpars.xlim; 
if isfield(intervalpars,'ylim') && ~isempty(intervalpars.ylim); hax.YLim = intervalpars.ylim; end

if pars.save_figure
    hax.Title.Visible = 'off';
    hax.XTick = []; 
% % %     interval_savename = [pars.interval_savename, strrep(thisvar,'corcoef_',''), '_', num2str(intervalpars.umdist_minmax(1)), '-', num2str(intervalpars.umdist_minmax(2)), 'um'];
    interval_savename = [pars.interval_savename, strrep(thisvar,'corcoef_',''), titlestring]; 
    saveas(interval_fig, interval_savename, 'svg')
    hax.Title.Visible = 'on';
end

[h, p_interpatch_mixed] = ttest2(sorted_vals{1,intervalpars.var_to_plot}{:}, sorted_vals{2,intervalpars.var_to_plot}{:});
[h, p_interpatch_patch] = ttest2(sorted_vals{1,intervalpars.var_to_plot}{:}, sorted_vals{3,intervalpars.var_to_plot}{:});
[h, p_mixed_patch] = ttest2(sorted_vals{2,intervalpars.var_to_plot}{:}, sorted_vals{3,intervalpars.var_to_plot}{:});
fprintf([titlestring, '...pvals... IP-mixed ', num2str(p_interpatch_mixed), '... IP-P ', num2str(p_interpatch_patch), '... mixed-P ', num2str(p_mixed_patch) '\n'])