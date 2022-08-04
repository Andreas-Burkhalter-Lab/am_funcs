%%% get average responses of cells (to preferred sf/tf/orient) with vs. without significant RF responses
%%%
%%% roitable from analyze_multi_planes.m must be loaded first
%%% updated 2020/9/23 on thermaltake

clear groupstats ops
groupstats = table; 
ops.meets_criteria = true(height(roitable),1);
% load('F:\analyses\cells20-9-21')

close all

%% parameters

% stimparam = 'sf'; 
% stimparam = 'tf'; 
stimparam = 'orient'; 

save_fig = 0;
    savename = 'F:\thesis\figs\2020-9-15\fig s1 -- center_vs_surround_response_amp';

% % % % % rf-tuned vs. non-rf-tuned cells
% groupstats.labels{1} = roitable.anovap_rf < 0.05; % rf tuned
% groupstats.labels{2} = roitable.anovap_rf > 0.05; % not rf tuned

% % % % % stim-on-rf vs stim-off-rf
groupstats.labels{1} = roitable.stim_overlap_frac > 0.00; % stim on rf
groupstats.labels{2} = roitable.stim_overlap_frac == 0.00; % stim on rf

% require significant tuning for variable of interest
% ops.meets_criteria = ops.meets_criteria & roitable{:, [stimparam, '_anovap']} < 0.05; 
ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} ]; 

% ops.meets_criteria = ops.meets_criteria & roitable.pakan_loc_sgnf; 
% ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac > 0.00; % stim on RF

ops.meets_criteria = ops.meets_criteria & roitable.inpatch; % patch ROIs
% ops.meets_criteria = ops.meets_criteria & ~roitable.inpatch; % interpatch ROIs


%% bargraph plot options
ops.newfig = 1; 
ops.ebar_thickness = 3; 
ops.ebar_cap_size = 6;
ops.ebar_color = [0 0 0]; 
ops.font_name = 'Arial'; 
ops.axis_thickness = 3.5;
ops.bar_quant_colors = repmat(linspace(0.08,0.86,6)',1,3);
    ops.quant_bar_width = 0.9; %%% bar width for quantile plotting only
ops.bar_one_color = [0.5 0.5 0.5]; 
ops.bar_quant_border_color = [1 1 1]; 
ops.bargraph_line_width = 4;
ops.axis_font_bold = 1; % font bold or not bold
ops.axis_font_size = 27;
ops.xlim = [0.3 2.7]; 

%     ops.ylim = [0.0 0.16]; 
% ops.fig_outerposition = [0.05 0.05 0.22 0.60]; % [left bottom width height]... fig 2 pakanloc by module
%     ops.fig_outerposition = [0.05 0.05 0.42 0.65]; % [left bottom width height]... fig 2 orient_si, rf_sgnf
    ops.fig_outerposition = [0.05 0.05 0.37 0.65]; % [left bottom width height]... fig 2 orient_si, ~rf_sgnf

    %%% don't count cells that fall in areas excluded from patchiness analysis, which are coded as quantile==0
ops.meets_criteria = ops.meets_criteria & roitable.quantile>0;


    
%% analysis
isTableCol = @(t, thisCol) ismember(thisCol, t.Properties.VariableNames);
if ~isTableCol(roitable,'sf_dff_epoch_mean') % only compute if this variable is not already in the table
    pars_to_analyze = {'sf', 'tf', 'orient'};
    for ipar = 1:length(pars_to_analyze)
        thispar = pars_to_analyze{ipar};
        dff_epochs_mean_full = [];
        for epoch = 1:3
            dff_epoch = roitable{:,[thispar, '_timecourse_mean_best']};
            dff_epoch = dff_epoch(:,epoch); 
            dff_epoch_mean = cellfun(@(x)mean(x),dff_epoch);
            dff_epochs_mean_full = [dff_epochs_mean_full, dff_epoch_mean];
        end
        try roitable{:,[thispar '_dff_epoch_mean']} = []; end %% clear variable if necessary
        roitable{:,[thispar '_dff_epoch_mean']} = dff_epochs_mean_full;
    end
end

values_to_test = []; % for anova
labels_to_test = []; % for anova
for igroup = 1:height(groupstats)
    matchrows = groupstats.labels{igroup} & ops.meets_criteria;
    dff_epochs_mean_full = roitable{:,[stimparam '_dff_epoch_mean']}; 
    groupstats{igroup,[stimparam '_dff_prestim']} =  {dff_epochs_mean_full(matchrows,1)}; % 1st column = dff before stim
    groupstats{igroup,[stimparam '_dff_during_stim']} =  {dff_epochs_mean_full(matchrows,2)}; % 2nd column = dff during stim
    dff_during_stim_minus_prestim = dff_epochs_mean_full(matchrows,2) - dff_epochs_mean_full(matchrows,1); % dff during stim minus dff prestim
    groupstats.response_mean(igroup) = mean(dff_during_stim_minus_prestim);
    groupstats.response_sem(igroup) = std(dff_epochs_mean_full(matchrows,2)) / sqrt(nnz(matchrows));
    values_to_test = [values_to_test; dff_during_stim_minus_prestim];
    labels_to_test = [labels_to_test; igroup*ones(nnz(matchrows),1)];
end

%% plotting
if ops.newfig; hfig = figure('units','normalized','outerposition',ops.fig_outerposition); end
hax = gca; 

% shade the bars individually if analyzing quantile
barhandle = bar(groupstats.response_mean); 
barhandle.FaceColor = ops.bar_one_color; 
barhandle.LineWidth = ops.bargraph_line_width; 
hax.XTick = []; 
hold on
ebar = errorbar(1:height(groupstats),groupstats.response_mean,groupstats.response_sem,'.','Color',ops.ebar_color);
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

if save_fig
    saveas(hfig, savename, 'svg') %%% save image as file
end

%% stats
[anovap, anovatab, anovastats] = anova1(values_to_test, labels_to_test, 'off');
anovap
anovastats