%%%% updated 20/1/25

% load('F:\analyses\cor5-21.mat','roitable','filetable') %%%% load this dataset before running the commands below

%% plot overall pakan loc index
% need to modify script slightly to override presets
ops.plot_bar_with_sem = 1; 
ops.plot_anova = 0;

%%%% miscellaneous variables of interest
% ops.varname = 'locm_cor_r';
ops.varname = 'pakan_loc_index';

%%%% criteria variables related to stim response
ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} | roitable{:,'rspv'} ]; % if roi is responsive or tuned for sf/tf/orient
% ops.meets_criteria = ops.meets_criteria & [ roitable{:,'sf_sgnf'} | roitable{:,'tf_sgnf'} | roitable{:,'orient_sgnf'} | roitable{:,'rf_sgnf'} ]; % if roi is tuned for sf/tf/orient/rf
ops.meets_criteria = ops.meets_criteria & roitable{:,'rf_sgnf'}; % if roi has detectable receptive field
% ops.meets_criteria = ops.meets_criteria & roitable{:,'rspv'};

%%%% criteria variables related to receptive field
% ops.meets_criteria = ops.meets_criteria & roitable{:,'stim_on_rf'};
ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac > 0.0; 
% ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac == 0; 
% ops.meets_criteria = ops.meets_criteria & roitable.rf_overlap_frac_aprx > 0.1;

%%%% specify sorting variable
ops.sorting_var = 'inpatch';

tuning_stats;
% ylabel('Locomotion modulation index')
% set(gca,'XTickLabel',{'Interpatch','Patch'})
title('')
xlabel('')
ylabel('')
set(gca,'Box','off')
hbar = findall(gcf,'Type','bar'); 
hbar.FaceColor = [0.5 0.5 0.5];
set(gca,'XTick',[])
set(gca,'linewidth',1) %%% axis line width
set(gca,'LooseInset',get(gca,'TightInset'))
print(gcf,'pakan_locmod_bargraph.tif', '-dtiffn','-r300')

%% plot locomotion mod examples
% quantile 1, high loc modulation; xy prereg = 34, 283.... pakan loc mod index = 0.2714; original row = 22
load('F:\18164\2019-01-16\tuningdat')
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
axis_font_size = 15; 
high_locmod_fig = figure;
[a, b] = sort(tuningdat.tuning.quantile);
tuningdat.tuning = tuningdat.tuning(b,:); %%% sort by quantile to make patch cells easier to find
plot_tuning_curve(tuningdat,7,'sf',0.003) % .003 = locomotion threshold
% title('Interpatch cell')
title('')
xlabel('')
ylabel('')
hLeg=findobj(gcf,'type','legend');
    set(hLeg,'visible','off')
set(gca,'Box','off')
set(gca,'linewidth',axes_line_width)
set(gca,'FontSize',axis_font_size)
set(gca,'FontWeight',axes_numbers_bold)
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
print(high_locmod_fig,'locmod-high_sfcurve.tif', '-dtiffn','-r300')


% quantile 5, low loc modulation; xy prereg = 118, 228..... pakan loc mod index = 0.0749
%%%%% other option: xy prereg = 26, 310; session 19-1-18
load('F:\19002\2019-01-18\tuningdat')
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
axis_font_size = 15; 
low_locmod_fig = figure;
[a, b] = sort(tuningdat.tuning.quantile);
tuningdat.tuning = tuningdat.tuning(b,:); %%% sort by quantile to make patch cells easier to find
plot_tuning_curve(tuningdat,67,'sf',.0001)
% title('Patch cell')
title('')
xlabel('')
ylabel('')
hLeg=findobj(gcf,'type','legend');
    set(hLeg,'visible','off')
set(gca,'Box','off')
set(gca,'linewidth',axes_line_width)
set(gca,'FontSize',axis_font_size)
set(gca,'FontWeight',axes_numbers_bold)
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
print(low_locmod_fig,'locmod-low_sfcurve.tif', '-dtiffn','-r300')


% plot timecourse of response for preferred sfreq on high loc vs low loc trials
%%%% code copied from get_tuning_curves, load_f_behavior_stim_data, get_f_behavior_stim
locm_thresh = 0.003;% %% use 0.0005 for 18165\2019-12-30 cellind 2; 6 subthresh, 4 superthresh trials
cellind =2; % % 18165\2019-12-30 high loc modulation; xy prereg = 73 13
line_width = 2; %%%% line for response timecourse
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
axis_font_size = 15;
% load('F:\18165\2018-12-30\tuningdat')
%%% other possible cells to use: % row 22 18164\2019-01-16, xy prereg = 34, 283; 18165 12/30 row 7, 18164 12/21 row 5; 18164 1/171 rows 11, 21
% % % % % extract SF trials data for this cell for its best tested sfreq
thiscell_sfbest_trials = tuningdat.tuning.sf_trials{cellind}(tuningdat.tuning.sf_trials{cellind}.sf==tuningdat.tuning.sfbest_tested_val(cellind),:);
sfbest_trials_loc = tuningdat.tuningpars.locm_trials.sf(tuningdat.tuning.sf_trials{cellind}.sf==tuningdat.tuning.sfbest_tested_val(cellind),:);
superthresh_trials = find(sfbest_trials_loc.locm_forw_mps > locm_thresh);
subthresh_trials = find(sfbest_trials_loc.locm_forw_mps <= locm_thresh);
thissf_stimpar_sets = tuningdat.tuningpars.stimpar_sets(tuningdat.tuningpars.stimpar_sets.sfreq==tuningdat.tuning.sfbest_tested_val(cellind),:);
ntrials_sfbest = height(thissf_stimpar_sets);
prairie_trigger_sampling_rate_hz = 10000; 
trig_sampling_interval_us = 1e6 / prairie_trigger_sampling_rate_hz;
scope_events = tuningdat.scope_events;
stimpar_sets = tuningdat.tuningpars.stimpar_sets; 
stimpars = tuningdat.stimpars; 
stimdur_us = stimpars.stimdur * 1e6; % convert from s to microseconds
stimdur_scans = round(stimdur_us / trig_sampling_interval_us);
stimdur_scope_events = tuningdat.tuningpars.stimdur_scope_events;
isi_scope_events = tuningdat.tuningpars.isi_scope_events;
scope_sampling_interval_sec = median(diff(tuningdat.scope_events.onset)) / prairie_trigger_sampling_rate_hz; % for timecourse x axis
for itrial = 1:ntrials_sfbest
        sevents_in_range_stim = find( scope_events.onset>thissf_stimpar_sets.stim_onset(itrial) & scope_events.onset<[thissf_stimpar_sets.stim_onset(itrial)+stimdur_scans] );
        sevents_in_range_stim = sevents_in_range_stim(1:end-1); % last image duration will be cut off by stim ending
        stim_first_scan = min(sevents_in_range_stim); % row index in scope_events
        sevents_in_range_stim_fixed_duration = stim_first_scan : stim_first_scan+stimdur_scope_events-1; % stim events in this trial using a fixed stim duration
        stim_last_scan = max(sevents_in_range_stim); % row index in scope_events
        % use a fixed number of scans for prestim, stim, and poststim periods to make peristimulus timecourse plots consistent
        sevents_in_range_prestim = stim_first_scan-isi_scope_events : stim_first_scan-1; 
        sevents_in_range_poststim = stim_last_scan+1 : stim_last_scan+isi_scope_events; % scope events during isi between this trial and next trial
        thissf_stimpar_sets.timecourse_prestim(itrial,:) = num2cell(scope_events.dff(sevents_in_range_prestim,cellind), 1);
        thissf_stimpar_sets.timecourse_stim(itrial,:) = num2cell(scope_events.dff(sevents_in_range_stim_fixed_duration,cellind), 1);
        thissf_stimpar_sets.timecourse_poststim(itrial,:) = num2cell(scope_events.dff(sevents_in_range_poststim,cellind), 1);
end
timecourses_subthresh = [cell2mat(thissf_stimpar_sets.timecourse_prestim(subthresh_trials)')', cell2mat(thissf_stimpar_sets.timecourse_stim(subthresh_trials)')',...
    cell2mat(thissf_stimpar_sets.timecourse_poststim(subthresh_trials)')'];
timecourses_superthresh = [cell2mat(thissf_stimpar_sets.timecourse_prestim(superthresh_trials)')', cell2mat(thissf_stimpar_sets.timecourse_stim(superthresh_trials)')',...
    cell2mat(thissf_stimpar_sets.timecourse_poststim(superthresh_trials)')'];
timecourse_x_axis = ( 0:scope_sampling_interval_sec:scope_sampling_interval_sec*[size(timecourses_subthresh,2)-1] ) - tuningdat.stimpars.isi;
%%% shift x axis back slightly for 18165\2019-12-30 cellind 2... onset and offset appear to come early (maybe due to cropped dff timepoints at end of timecourse segmeents)
timecourse_x_axis = timecourse_x_axis + scope_sampling_interval_sec;
close all
timecourse_locomotion_fig = figure;
subthresh_plot = plot(timecourse_x_axis,mean(timecourses_subthresh),'k');
subthresh_plot.LineWidth = line_width;
hold on
superthresh_plot = plot(timecourse_x_axis,mean(timecourses_superthresh),'r');
superthresh_plot.LineWidth = line_width;
% legend({'low locomotion','high locomotion'})
hLeg=findobj(gcf,'type','legend');
    set(hLeg,'visible','off')
set(gca,'Box','off')
set(gca,'linewidth',axes_line_width)
set(gca,'FontSize',axis_font_size)
set(gca,'FontWeight',axes_numbers_bold)
set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
print(timecourse_locomotion_fig,'timecourse_locomotion_fig.tif', '-dtiffn','-r300')

%% plot sf hwhm patches vs interpatches
% need to modify script slightly to override presets
ops.plot_bar_with_sem = 1; 
axes_line_width = 2;
axes_numbers_bold = 1; %%% make axis labels bold
ops.varname = 'sf_hwhm';
ops.meets_criteria = ops.meets_criteria & roitable{:,'sf_sgnf'};
ops.meets_criteria = ops.meets_criteria & roitable{:,'rf_sgnf'}; % if roi has detectable receptive field
ops.meets_criteria = ops.meets_criteria & roitable.stim_overlap_frac > 0.0; 
ops.sorting_var = 'inpatch';
tuning_stats
% ylabel('HWHM (cycle/deg)')
% set(gca,'XTickLabel',{'Interpatch', 'Patch'})
title('')
xlabel('')
ylabel('')
set(gca,'Box','off')
hbar = findall(gcf,'Type','bar'); 
hbar.FaceColor = [0.5 0.5 0.5];
set(gca,'XTick',[])
set(gca,'linewidth',axes_line_width)
set(gca,'FontWeight',axes_numbers_bold)
set(gca,'LooseInset',get(gca,'TightInset'))
print(gcf,'sf hwhm bargraph.tif', '-dtiffn','-r300')


%%
% % % %% plot timecourses by quantile
% % % % first need to run tuning_stats to get meets_criteria, requiring significant sf tuning and that stimulus is on RF
% % % load('cor5-21')
% % % 
% % % figure
% % % q = cell2mat(cellfun(@(x)x(1:8),roitable.sftimecourse_mean_best(:,2),'UniformOutput',0));
% % % sf12 = q(ops.meets_criteria & [roitable.quantile==1|roitable.quantile==2]& roitable.sf_anovap<1e-6,:);
% % % sf34 = q(ops.meets_criteria & [roitable.quantile==3|roitable.quantile==4]& roitable.sf_anovap<1e-6,:);
% % % sf56 = q(ops.meets_criteria & [roitable.quantile==5|roitable.quantile==6]& roitable.sf_anovap<1e-6,:);
% % % sfp = q(ops.meets_criteria & roitable.inpatch & roitable.sf_anovap<1e-6,:);
% % % sfip = q(ops.meets_criteria & ~roitable.inpatch & roitable.sf_anovap<1e-6,:);
% % % hold on
% % % s= sf12; errorbar(1:8,mean(s),std(s)./sqrt(size(s,1)))
% % % s= sf34; errorbar(1:8,mean(s),std(s)./sqrt(size(s,1)))
% % % s= sf56; errorbar(1:8,mean(s),std(s)./sqrt(size(s,1)))
% % % 
% % % legend({'1-2','3-4','5-6'})

% % % %% plot pakan loc index by day
% % % plot_anova = 'off';
% % % for i = 1:13
% % %     other_criteria = strcmp(roitable.day, filetable.day{i});
% % %     tuning_stats
% % %     subplot(4,4,i)
% % %     bar(sorted_vals.mean)
% % %     hold on
% % %     errorbar(1:height(sorted_vals),sorted_vals.mean,sorted_vals.sem,'.')
% % %     title([filetable.day{i}, ', p=' num2str(p(2,1))])
% % %     hold off
% % % end








