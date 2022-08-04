%%% plot pairwise correlations while varying multiple parameters
%
%%% cortable must be loaded first
%%% include_temp will be determined in this upper level script
%
%%% % updated 2022/3/9 on thermaltake

if ~exist('cortable','var')
    load('F:\analyses\cor2022-03-14_for_plotting.mat')
    cortable = cortable_out;
end

close all
clear pars include_temp matchrows labels_to_plot vals_to_plot ebar
include_temp = true(height(cortable),1); % initialize
titlestring = ''; 

%%
plot_modules_overlaid = 1; % plot P-P, P-IP, IP-IP on the same axis... leave on if running interval plotting
    pars.plot_noise = 0; 
    pars.plot_rest = 1; 
    %%%%%%% run_sorting_function must be turned on for changes to include_temp to take effect
    pars.run_sorting_function = 1; %% run pairwise_cor_plotting.m (slow); if false, use the data already in the workspace
    pars.overlay_subplots = 0; % plot overlaid axes in subplots
    pars.curves_colorlist = {[0 0.6 0] [0 0 0] [0.8 0 1] [1 1 0] [1 0 1] [ 0 1 1]};
    pars.star_vspace_proportion = 0.05; % space between error bar top and signficance stars, as proportion of range(ylim)
    pars.star_color = [0 0 0]; 
    pars.star_font_size = 25; 
    pars.star_font_name = 'Arial'; 
    pars.curve_thickness = 7; % LineWidth
    pars.ebar_thickness = 3; 
    pars.distcurve_axis_thickness = 3;
    pars.fig_outerposition = [0.05 0.05 0.9 0.9]; % [left bottom width height]
%     pars.rest_ylim = [0.01 0.115]; 
%     pars.noise_ylim = [0.01 0.175]; 
%         pars.noise_ylim = [0.0 0.125]; 
%         pars.noise_ylim = [0.01 0.3]; 
        pars.noise_ylim = [0.0 0.22]; 
    pars.save_figure = 0; 
        pars.distcurve_savename = 'F:\thesis\figs\fig 3 -- cor_by_distance_'; % rest/noise will be appended
        pars.interval_savename = 'F:\thesis\figs\fig 3 -- cor_bargraph_'; % rest/noise will be appended

% interval plotting options 
run_interval_plotting = 1; %%% plot aggregated cor values from a specific distance interval (make sure plot_modules_overlaid==true)
    intervalpars.axis_thickness = 5; 
    intervalpars.show_boxplot = 0; 
        intervalpars.notchstyle = 'marker';
    %     intervalpars.notchstyle = 'on';
    %     intervalpars.notchstyle = 'off';
        intervalpars.outlier_symbol = 'rx';
    intervalpars.show_bargraph = 1;
        intervalpars.bargraph_line_width = 5;
        intervalpars.barerror_line_width = 3; % set errorbar line width
        intervalpars.error_line_color = [0 0 0];    %%%%% repmat({[0 0 0]},1,length(pars.curves_colorlist)); % set errorbar line color
        intervalpars.umdist_minmax = [250 425]; 
%         intervalpars.bar_colors = cellfun(@(x) x*1.6, pars.curves_colorlist,'UniformOutput',0); % use lightened curve colors for bars
        intervalpars.bar_colors = [pars.curves_colorlist(1), {0.2*ones(1,3)}, pars.curves_colorlist(3)];
        intervalpars.border_colors = repmat({[0 0 0]},1,length(pars.curves_colorlist));
        intervalpars.xlim = [0.25 3.6];
        intervalpars.ylim = [0 0.2]; 
        intervalpars.fig_outerposition = [0.05 0.05 0.4 0.7]; % [left bottom width height]
    % intervalpars.umdist_minmax = [150 180]; 
    intervalpars.var_to_plot = 'corcoef_rest'; % NB: this field doesn't have much effect when running cor_plotting_multi_groups
%     intervalpars.var_to_plot = 'corcoef_noise';% NB: this field doesn't have much effect when running cor_plotting_multi_groups

% general options
pars.show_title = 1; 
pars.font_name = 'Arial';
pars.figs_docked = 0;
pars.axis_font_bold = 1; % font bold or not bold
pars.axis_font_size = 25;

% analysis options
analysis_pars.plot_bins = 1;
analysis_pars.nbins = 17;
analysis_pars.grouping_mode = 'quantile_groups'; % group quantiles together then plot based on distances between these groups
    analysis_pars.sort_by_quant_group = 1; %%% sort plot results by the quant groups specified (recommended for clarity)... if false, sort by group distances
    analysis_pars.quantile_groups_to_plot = {[1 2 3]  [4 5 6]}; % top half vs bottom half
analysis_pars.use_umdist_threshold = 1;  % if true, make sure that all category distances are drawing from the same range of umdists
analysis_pars.umdist_max = 425; %%% do not analyze or plot pairs that are with a greater distance in um separating them
analysis_pars.hold_plots =0; 
analysis_pars.newfig = 0; 

% separate plots using patches and interpatches as seed quantile groups... function not maintained
plot_modules_separately = 0; 
    subplot_rowcol = [2 2]; 

% include_temp = include_temp & all(cortable.sf_sgnf,2); titlestring = [titlestring, ' sf_sgnf']; % require sf tuned
% include_temp = include_temp & all(cortable.tf_sgnf,2); titlestring = [titlestring, ' tf_sgnf']; % require tf tuned
% include_temp = include_temp & all(cortable.orient_sgnf,2); titlestring = [titlestring, ' orient_sgnf']; % require orient tuned
% include_temp = include_temp & all(cortable.orient_anovap<0.1,2); titlestring = [titlestring, ' orient p<0.01']; % require orient weakly tuned
include_temp = include_temp & all(cortable.sf_sgnf | cortable.tf_sgnf | cortable.orient_sgnf | cortable.rspv_pval<0.05, 2); % tuned for sf/tf/orient or general responsiveness
% include_temp = include_temp & all(cortable.sf_sgnf | cortable.tf_sgnf | cortable.orient_sgnf | cortable.rf_sgnf | cortable.rspv_pval<0.05, 2); % tuned for sf/tf/orient/rf or general responsiveness
% include_temp = include_temp & all(cortable.rf_sgnf,2); [titlestring, ' rf_sgnf']; % require both cells to be rf tuned

% include_temp = include_temp & all(cortable.stim_on_rf,2); titlestring = [titlestring, ' stim_on_rf']; %%% both cells had the stimulus overlapping with their RFs
%     include_temp = include_temp & all(~cortable.stim_on_rf,2); titlestring = [titlestring, ' ~stim_on_rf']; %%% neither cell had the stimulus overlapping with their RF
% include_temp = include_temp & all(cortable.locm_sgnf,2); titlestring = [titlestring, ' locm_sgnf']; %%% both cells in pair locomotion responsive in darkness
%   include_temp = include_temp & all(~cortable.locm_sgnf,2); titlestring = [titlestring, ' ~locm_sgnf']; %%% both cells in pair locomotion non-responsive in darkness
% include_temp = include_temp & all(cortable.rf_width_deg<19,2);  titlestring = [titlestring, ' rf_size']; 

%% run cor plotting for patches, interpatches, rest cor, stim noise cor
if pars.figs_docked; windowstyle = 'docked'; else windowstyle = 'normal'; end
analysis_pars.include = include_temp; %%%% apply inclusion criteria

% plot 2 axes overlaying p-p, p-ip, ip-ip
if plot_modules_overlaid
    pars.hold_plots = 1; 
    pars.plot_bins = 0; 
    ylims = {field_default(pars,'rest_ylim',[]), field_default(pars,'noise_ylim',[])}; 
    plot_titles = {'rest', 'noise'}; 
    varnames = {'corcoef_rest','corcoef_noise'}; 
    quant_groups_temp = analysis_pars.quantile_groups_to_plot;
    if pars.plot_rest | pars.plot_noise
        overlay_fig(1) = figure('units','normalized','outerposition',pars.fig_outerposition,'WindowStyle',windowstyle); 

    end

    % run sorting function if required
    if pars.run_sorting_function
        clear plotdat data_out cortable_out
        for ivar = 1:2 % for corcoef_rest and corcoef_noise
            analysis_pars.var_to_plot = varnames{ivar}; 
            analysis_pars.seedquant = [1 2 3]; % interpatches
            [cortable_out,data_out{1}] = pairwise_cor_plotting(cortable,analysis_pars);
            patch_patch_pairs = [cortable_out.roi1group==2 & cortable_out.roi2group==2]; 
            cortable_out.quantgroup_pair(patch_patch_pairs) = 3; %%% fill in patch-patch pair labels to the cortable
            analysis_pars.seedquant = [4 5 6]; % patches
            [~,data_out{2}] = pairwise_cor_plotting(cortable,analysis_pars);
            plotdat{ivar} = data_out{1}.cor_by_umdist; % copy the data for IP-IP and P-IP pairs
            plotdat{ivar}.bin_centers = data_out{ivar}.bin_centers';   plotdat{ivar} = movevars(plotdat{ivar},'bin_centers','Before','umdist');
            plotdat{ivar}.cor(:,3) = data_out{2}.cor_by_umdist.cor(:,2); %%% add p-p to ip-ip and p-ip data
            plotdat{ivar}.cormean(:,3) = data_out{2}.cor_by_umdist.cormean(:,2); %%% add p-p to ip-ip and p-ip data
            plotdat{ivar}.cor_sem(:,3) = data_out{2}.cor_by_umdist.cor_sem(:,2); %%% add p-p to ip-ip and p-ip data
            for ibin = 1:height(plotdat{ivar}) %%% get statistical significance for each bin
                [h, plotdat{ivar}.pval_patch_interpatch(ibin)] = ttest2(plotdat{ivar}.cor{ibin,1}, plotdat{ivar}.cor{ibin,3});
            end
            % assign significance stars
            plotdat{ivar}.sgnf_patch_interpatch = zeros(height(plotdat{1}),1);
            plotdat{ivar}.sgnf_patch_interpatch(plotdat{ivar}.pval_patch_interpatch < 0.05) = 1;
            plotdat{ivar}.sgnf_patch_interpatch(plotdat{ivar}.pval_patch_interpatch < 0.01) = 2;
            plotdat{ivar}.sgnf_patch_interpatch(plotdat{ivar}.pval_patch_interpatch < 0.001) = 3;
            plotdat{ivar}.sgnf_patch_interpatch(plotdat{ivar}.pval_patch_interpatch < 0.0001) = 4;  
        end
        plotdat{length(plotdat)+1} = analysis_pars; %%% record how the analysis was performed
    end
    % do plotting
    nplots = nnz([pars.plot_rest, pars.plot_noise]); 
    for ivar = 1:2
        thisvar = varnames{ivar};
        if [pars.plot_rest && ivar==1] || [pars.plot_noise && ivar==2] %%% plot only the specified variables
            if pars.overlay_subplots
                subplot(1, nplots, min(ivar, nplots))
            elseif ivar > 1 && nplots > 1
                overlay_fig(ivar) = figure('units','normalized','outerposition',pars.fig_outerposition,'WindowStyle',windowstyle);
            end
            hplot{ivar} = plot(repmat(plotdat{ivar}.bin_centers,1,size(plotdat{ivar}.cor,2)), plotdat{ivar}.cormean);
            hold on;
            ebar{ivar} = errorbar(repmat(plotdat{ivar}.bin_centers,1,size(plotdat{ivar}.cor,2)), plotdat{ivar}.cormean, plotdat{ivar}.cor_sem, plotdat{ivar}.cor_sem,'LineStyle','none'); 
            linewidths = num2cell(repmat(pars.curve_thickness,1,length(hplot{ivar}))); [hplot{ivar}.LineWidth] = linewidths{:}; % set line width
            ebarlinewidths = num2cell(repmat(pars.ebar_thickness,1,length(ebar{ivar}))); [ebar{ivar}.LineWidth] = ebarlinewidths{:}; % set line width
            [ebar{ivar}.Color] = pars.curves_colorlist{1:length(ebar{ivar})}; 
            [hplot{ivar}.Color] = pars.curves_colorlist{1:length(hplot{ivar})};
            hax = gca;
            hax.FontSize = pars.axis_font_size;
            hax.FontName = pars.font_name;
            hax.LineWidth = pars.distcurve_axis_thickness;
            if pars.axis_font_bold; set(gca,'FontWeight','bold'); end
            if ~isempty(ylims{ivar}); ylim(ylims{ivar}); else ylims{ivar} = ylim; end %%% set ylimits if specified
            xlim([0, analysis_pars.umdist_max+5]);
            if pars.show_title; title([plot_titles{ivar}, titlestring]); end
            % add significance stars
            maxbars = max(plotdat{ivar}.cormean + plotdat{ivar}.cor_sem, [], 2); % highest errorbar in each bin
            stars_vspace = pars.star_vspace_proportion * range(ylims{ivar}); 
            stars_yloc = maxbars + stars_vspace; % significance stars location
            for istar = 1:4
                starrows = plotdat{ivar}.sgnf_patch_interpatch == istar;
                nstars = nnz(starrows);
                textscat = textscatter(plotdat{ivar}.bin_centers(starrows), stars_yloc(starrows), repmat({repmat('*',1,istar)},nstars,1), 'TextDensityPercentage',100);
                textscat.ColorData = pars.star_color; 
                textscat.FontSize = pars.star_font_size;
                textscat.FontName = pars.star_font_name; 
            end
            set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
            set(gca,'Box','off')
            if pars.save_figure
                hax.Title.Visible = 'off';
                distcurve_savename = [pars.distcurve_savename, strrep(thisvar,'corcoef_','')];
                saveas(overlay_fig, distcurve_savename, 'svg')
                hax.Title.Visible = 'on';
            end
            hold off
            npairs_per_group = sum(cellfun(@(x)length(x),plotdat{2}.cor)); % total IP-IP, P-IP, P-P pairs plotted
            fprintf([thisvar, titlestring, '...n pairs... IP-IP ', num2str(npairs_per_group(1)), '... IP-P ', num2str(npairs_per_group(2)), '... P-P ', num2str(npairs_per_group(3)) '\n'])
        end
    end
end

if run_interval_plotting
    cor_interval_plotting();
end

% get cor data when pairs are combined without regard to distance
plotted_ipip_rest = cell2mat(plotdat{1}.cor(:,1));
plotted_pip_rest = cell2mat(plotdat{1}.cor(:,2));
plotted_pp_rest = cell2mat(plotdat{1}.cor(:,3));
plotted_ipip_noise = cell2mat(plotdat{2}.cor(:,1));
plotted_pip_noise = cell2mat(plotdat{2}.cor(:,2));
plotted_pp_noise = cell2mat(plotdat{2}.cor(:,3));

cortable_alldist_rest = table([plotted_ipip_rest; plotted_pip_rest; plotted_pp_rest],...
    [repmat({'ipip'},length(plotted_ipip_rest),1); repmat({'pip'},length(plotted_pip_rest),1); repmat({'pp'},length(plotted_pp_rest),1)],...
    'VariableNames',{'cor','pairtype'});
cortable_alldist_noise = table([plotted_ipip_noise; plotted_pip_noise; plotted_pp_noise],...
    [repmat({'ipip'},length(plotted_ipip_noise),1); repmat({'pip'},length(plotted_pip_noise),1); repmat({'pp'},length(plotted_pp_noise),1)],...
    'VariableNames',{'cor','pairtype'});

[p_rest,tbl_rest,stats_rest] = anova1(cortable_alldist_rest.cor, cortable_alldist_rest.pairtype,'on');
[p_noise,tbl_noise,stats_noise] = anova1(cortable_alldist_noise.cor, cortable_alldist_noise.pairtype,'on');


%%%% plot 4 axes without overlaying p-p with ip-ip... may no longer be functional
if plot_modules_separately
    cor_multiplotting_separate_groups()
end

