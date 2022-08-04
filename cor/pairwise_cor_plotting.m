%%% plot umdist against corcoef
% cortable_out = pairwise_cor_plotting(cortable)
%%%%% updated 2020-02-27 on thermaltake

function [cortable_out,data_out] = pairwise_cor_plotting(cortable, pars)

% close all
pars = vardefault('pars', struct); 
include_temp = true(height(cortable), 1); % initialize (can be overwritten by pars.include)

% pars.var_to_plot = field_default(pars,'var_to_plot','corcoef_rest'); % dff cor during dark period
pars.var_to_plot = field_default(pars,'var_to_plot','corcoef_noise'); % noise cor during stim
pars.plot_points = field_default(pars,'plot_points',0); 
pars.plot_curves = field_default(pars,'plot_curves',0);
pars.plot_bins = field_default(pars,'plot_bins',1);
    pars.nbins = field_default(pars,'nbins',25);
pars.hold_plots = field_default(pars,'hold_plots',0);
pars.newfig = field_default(pars,'newfig',0); % make new figure

% include_temp = include_temp & all(cortable.sf_sgnf,2); % require sf tuned
% include_temp = include_temp & all(cortable.tf_sgnf,2); % require tf tuned
% include_temp = include_temp & all(cortable.orient_sgnf,2); % require orient tuned
% include_temp = include_temp & all(cortable.sf_sgnf | cortable.tf_sgnf | cortable.orient_sgnf | cortable.rspv_pval<0.05, 2); % tuned for sf/tf/orient or general responsiveness
include_temp = include_temp & all(cortable.sf_sgnf | cortable.tf_sgnf | cortable.orient_sgnf | cortable.rf_sgnf | cortable.rspv_pval<0.05, 2); % tuned for sf/tf/orient/rf or general responsiveness
% include_temp = include_temp & all(cortable.rf_sgnf,2); % require both cells to be rf tuned

% include_temp = include_temp & all(cortable.stim_on_rf,2); %%% both cells had the stimulus overlapping with their RFs
    include_temp = include_temp & all(~cortable.stim_on_rf,2); %%% neither cell had the stimulus overlapping with their RF
include_temp = include_temp & all(cortable.locm_sgnf,2); %%% both cells in pair locomotion responsive in darkness
%   include_temp = include_temp & all(~cortable.locm_sgnf,2); %%% both cells in pair locomotion non-responsive in darkness

% pars.grouping_mode = field_default(pars,'grouping_mode','quantdist_groups'); % plot groups of quantdists rather than individual quantdists
%         pars.quantdist_groups_to_plot = field_default(pars,'quantdist_groups_to_plot',{[0 1], [2 3 4 5]});
%     pars.quantdist_groups_to_plot = field_default(pars,'quantdist_groups_to_plot',{[0 1 2], [3 4 5]});
%     pars.quantdist_groups_to_plot = field_default(pars,'quantdist_groups_to_plot',{[0], [1 2 3 4 5]});
% %     pars.quantdist_groups_to_plot = field_default(pars,'quantdist_groups_to_plot',{[0 1 2 3] [4 5]});
% pars.grouping_mode = field_default(pars,'grouping_mode','single_quantdists'); % plot invididual quantdists
%     pars.quantdists_to_plot = field_default(pars,'quantdists_to_plot',[0: 5]);
pars.grouping_mode = field_default(pars,'grouping_mode','quantile_groups'); % group quantiles together then plot based on distances between these groups
    pars.sort_by_quant_group = 1; %%% sort plot results by the quant groups specified (recommended for clarity)... if false, sort by group distances
%     pars.quantile_groups_to_plot = field_default(pars,'quantile_groups_to_plot',{[1 2] [3 4] [5 6]}); % upper mid and lower quants
%     pars.quantile_groups_to_plot = field_default(pars,'quantile_groups_to_plot',{[1 2 3 4] [5 6]}); % reversing sincich and horton, assuming mouse interpatches = monkey blobs
%      pars.quantile_groups_to_plot = field_default(pars,'quantile_groups_to_plot',{[1 2] [3 4 5 6]}); % matching sincich and horton, assuming mouse interpatches = monkey blobs
     pars.quantile_groups_to_plot = field_default(pars,'quantile_groups_to_plot',{[1 2 3]  [4 5 6]}); % top half vs bottom half
%     pars.quantile_groups_to_plot = field_default(pars,'quantile_groups_to_plot',{[1] [2 3 4 5 6]}); % interpatch centers vs everything else

pars.use_umdist_threshold = field_default(pars,'use_umdist_threshold',1); % if true, make sure that all category distances are drawing from the same range of umdists
pars.umdist_max = field_default(pars,'umdist_max',700); %%% do not analyze or plot pairs that are with a greater distance in um separating them
pars.use_seedquant = field_default(pars,'use_seedquant',1); %%% if true, only look at roi pairs with one roi from one quant
%     pars.seedquant = field_default(pars,'seedquant',[1 2]); %%% only look at roi pairs with at least one roi from these quants
%     pars.seedquant = field_default(pars,'seedquant',[3 4 5 6]);
%     pars.seedquant = field_default(pars,'seedquant',[1 2 3]);
     pars.seedquant = field_default(pars,'seedquant',[4 5 6]);
%     pars.seedquant = field_default(pars,'seedquant',[1]);
pars.colorlist = field_default(pars,'colorlist',{[1 0 0] [0 1 0] [0 0 1] [1 1 0] [1 0 1] [ 0 1 1]});
pars.errorbar_color = field_default(pars,'errorbar_color',[0 0 0]); 
pars.linestyle = field_default(pars,'linestyle','-');
pars.ylimits = field_default(pars,'ylimits',[-0.05 0.1]); % coercoef rest only
% pars.xlimitright = max(cortable.umdist);
pars.xlimitright = field_default(pars,'xlimitright',pars.umdist_max); 
% pars.xlimitright = field_default(pars,'',500);
pars.add_title = field_default(pars,'add_title',1); 


%% session inclusion criteria 
% include_temp = include_temp & strcmp(cortable.day,... 
% ...   {'2018-12-21'}...
% ...    {'2018-12-22'}...
% ...    {'2019-01-16'}...
% ...    {'2019-01-17'}...
% ...    {'2018-12-30'}...
% ...    {'2019-01-04'}...
% ...    {'2019-01-07'}...
% ...    {'2019-01-18'}...
% ...    {'2019-01-19'}...
% ...    {'2019-01-20'}...
% ...    {'2019-01-21'}...
% ...    {'2019-03-30'}...
%  ...      {'2020-03-20'}...
%  ...        {'2020-03-22'}...
% ...        {'2020-03-24'}...
%         {'2020-03-25'}...
% ...    {'2020-04-01'}...
% ...     {'2020-04-07'} ...
% );

% include_temp = include_temp & cortable.sub ==...
% ...      18164 ...
%        18165 ...
% ...       19002 ...
% ...       19042 ...
% ...        19048 ...
% ...        19089 ...
% ...        20003 ...
% ...        20004 ...
% ...     20015 ...
% ;


%%
    
% if pars.include not specified in function input, use the vector built by defaults in this function
%       if pars.include was specified in function input, none of the inclusion criteria specified by include_temp will be applied
pars.include = field_default(pars,'include',include_temp);

if pars.use_seedquant  %%% only look at roi pairs with at least one roi from the listed seed quantiles
    matches_seedquant = [any(cortable.quantile1==pars.seedquant,2), any(cortable.quantile2==pars.seedquant,2)]; 
else
    matches_seedquant = true(height(cortable),2);
end
seedquant_check = any(matches_seedquant,2); 
pars.include = pars.include & seedquant_check; % exclude pairs not including the seedquant(s)

npairs = height(cortable);
switch pars.grouping_mode
    case 'quantdist_groups'
        ncurves_to_plot = length(pars.quantdist_groups_to_plot);
        category_labels = cortable.quantdist; 
        unq_catlabels = cell2mat(pars.quantdist_groups_to_plot); 
    case 'single_quantdists'
        ncurves_to_plot = length(pars.quantdists_to_plot); 
        category_labels = cortable.quantdist; 
        unq_catlabels = pars.quantdists_to_plot; 
    case 'quantile_groups'
        ncurves_to_plot = length(pars.quantile_groups_to_plot); 
        % find out which quant group(s) the seedquant(s) belong to
        for igroup = 1:length(pars.quantile_groups_to_plot)
            matchgroup(igroup) = ~isempty(intersect(pars.seedquant, pars.quantile_groups_to_plot{igroup}));
        end
        seedgroup = find(matchgroup); 
        cortable.roi1group = NaN(npairs,1);
        cortable.roi2group = NaN(npairs,1);
        cortable.quantgroup_pair = NaN(npairs,1);
        for igroup = 1:ncurves_to_plot
            these_quantiles = pars.quantile_groups_to_plot{igroup};
            roi1_is_this_group = arrayfun(@(x)any(x == these_quantiles),cortable.quantile1);
            roi2_is_this_group = arrayfun(@(x)any(x == these_quantiles),cortable.quantile2);
            cortable.roi1group(roi1_is_this_group) = igroup;
            cortable.roi2group(roi2_is_this_group) = igroup;
        end
        quant_group_pairs = [cortable.roi1group, cortable.roi2group]; 
        both_roi_seedgroup = quant_group_pairs(:,1) == seedgroup & quant_group_pairs(:,2) == seedgroup; % both rois are from the seedgroup
        neither_roi_seedgroup = quant_group_pairs(:,1) ~= seedgroup & quant_group_pairs(:,2) ~= seedgroup; % neither roi is from the seedgroup
        quant_group_pairs(neither_roi_seedgroup,:) = NaN; % neither roi is from the seedgroup, do not plot
        quant_group_pairs(quant_group_pairs==seedgroup) = NaN; % label the mixed groups by the roi NOT from the seedgroup
        cortable.quantgroup_pair = max(quant_group_pairs,[],2); % label the mixed groups by the roi NOT from the seedgroup
        cortable.quantgroup_pair(both_roi_seedgroup) = seedgroup; % label this pair with the seedgroup index
        cortable.groupdist = abs(cortable.roi1group - cortable.roi2group); % difference between group indices between each member of an roi pair
        category_labels = cortable.groupdist; 
        unq_catlabels = 0:ncurves_to_plot-1;
end
    
if pars.use_umdist_threshold
    for icat = 1:length(unq_catlabels)
        thiscatlabel = unq_catlabels(icat);
        minvals(icat) = min(cortable.umdist(category_labels==thiscatlabel & pars.include)); % min umdist in this category
    end
    umdist_thresh = max(minvals); % find the largest minimum umdist from all categories (from rois included in the analysis)
else
    umdist_thresh = 0;
end
umdist_check_min = cortable.umdist >= umdist_thresh; %% pairs that aren't below umdist
umdist_check_max = cortable.umdist <= pars.umdist_max; % pairs that don't exceed the specified maximum
pars.include = pars.include & umdist_check_min; % exclude pairs below umdist min
pars.include = pars.include & umdist_check_max; % exclude pairs above umdist max
xlimits = [umdist_thresh-10 pars.xlimitright];


%% make plots 
if pars.plot_curves
    % initialize plot
    set(gca,'XLim',xlimits)
    set(gca,'YLim',pars.ylimits)
end
if ~pars.hold_plots
    hold off
else
    hold on
end
if pars.newfig
    corfig = figure;
end
    
matchrows = false(height(cortable), ncurves_to_plot);
cor_by_umdist = table(NaN(pars.nbins,ncurves_to_plot), cell(pars.nbins,ncurves_to_plot), NaN(pars.nbins,ncurves_to_plot), NaN(pars.nbins,ncurves_to_plot),...
                        'VariableNames', {'umdist','cor','cormean','cor_sem'});
for icurve = 1:ncurves_to_plot
    clear binmeans
    switch pars.grouping_mode
        case 'quantdist_groups'
            these_quantdists = pars.quantdist_groups_to_plot{icurve}; 
            matchrows(:,icurve) = any(cortable.quantdist == these_quantdists, 2);
        case 'single_quantdists'
            this_quantdist = pars.quantdists_to_plot(icurve);
            matchrows(:,icurve) = cortable.quantdist == this_quantdist;
        case 'quantile_groups'
            if pars.sort_by_quant_group % sort by quant group of the non-seed pair member
                this_quantgroup = icurve;
                matchrows(:,icurve) = cortable.quantgroup_pair == this_quantgroup;
            elseif ~pars.sort_by_quant_group % sort by group distance within the pair
                this_groupdist = icurve-1;
                matchrows(:,icurve) = cortable.groupdist == this_groupdist;
            end
    end
    matchrows = matchrows & pars.include; % integrate previous exclusionary criteria
    umdist_vals{icurve} = cortable.umdist(matchrows(:,icurve)); 
    corcoef_vals{icurve} = cortable{matchrows(:,icurve),pars.var_to_plot};
    data_out.matchrows = matchrows;
    data_out.umdist_vals = umdist_vals;
    data_out.corcoef_vals = corcoef_vals;

    [xData, yData] = prepareCurveData( umdist_vals{icurve}, corcoef_vals{icurve} );
    clear h
    if pars.plot_points
        h = plot( xData, yData, '.');
        hold on
    end
    if pars.plot_curves
        ft = fittype( 'exp2' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.0375031660988949 -0.0131607246108098 0.000453282061531754 0.00728078156135307];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        h = plot(fitresult, 'LineStyle',pars.linestyle);
        hold on
        legend(bins_ax, num2str([1:ncurves_to_plot]'))
    end
    
    %%% get bin information even if we're not going to plot it
    [umdist_disc, edgevals] = discretize(umdist_vals{icurve},pars.nbins); % discretize distances
    cor_by_umdist.umdist(:,icurve) = edgevals(2:end); % umdist in the table cor_by_umdist is the largest umdist of pairs in this bin
    bin_centers = edgevals(1:end-1) + mean(diff(edgevals))/2; 
    for ibin = 1:pars.nbins
        cor_by_umdist.cor{ibin,icurve} = corcoef_vals{icurve}(umdist_disc==ibin); %% cor of all pairs that belong in this umdist bin
        cor_by_umdist.cor_sem(ibin,icurve) = nanstd(cor_by_umdist.cor{ibin,icurve}) / sqrt(length(cor_by_umdist.cor{ibin,icurve})); % standard error of the mean for this bin and curve           
        binmeans(ibin) = nanmean( cor_by_umdist.cor{ibin,icurve} ); % mean of corcoefs for pairs with distances within this bin
        ebar = errorbar(bin_centers(ibin), binmeans(ibin), cor_by_umdist.cor_sem(ibin,icurve), cor_by_umdist.cor_sem(ibin,icurve)); 
        ebar.Color = pars.errorbar_color;
        hold on
    end
    cor_by_umdist.cormean(:,icurve) = binmeans';
    if pars.plot_bins
        bins_ax(icurve) = plot(bin_centers, binmeans, 'LineStyle',pars.linestyle);
        hold on
        legend(bins_ax, num2str([1:ncurves_to_plot]'))
    end
    bins_ax(icurve).Color = pars.colorlist{icurve};

    %     scatter(cortable.umdist(matchrows), cortable.cortable{matchrows,pars.var_to_plot})
end

data_out.cor_by_umdist = cor_by_umdist;
data_out.bin_centers = bin_centers; 
data_out.mean = cellfun(@nanmean,data_out.corcoef_vals);
if length(data_out.umdist_vals) == 2 % if we're comparing exactly two groups [usually patches vs. interpatches]
    [data_out.sgnf, data_out.pval] = ttest2(data_out.corcoef_vals{1}, data_out.corcoef_vals{2});
    title_append = [', p = ' num2str(data_out.pval)];
else
    title_append = '';
end


if strcmp(pars.var_to_plot,'corcoef_rest')
    set(gca,'YLim',pars.ylimits)
end

if pars.add_title
    title([pars.var_to_plot, ' ', pars.grouping_mode, title_append])
end
if ~pars.hold_plots
    hold off
end
cortable_out = cortable;
xlabel('distance (um)')
ylabel('correlation coefficient')