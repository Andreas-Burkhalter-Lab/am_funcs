% plot timecourse of response for preferred sfreq on high loc vs low loc trials
%%%% code copied from get_tuning_curves, load_f_behavior_stim_data, get_f_behavior_stim
%
%%%% updated 20/9/25

clear pars

locm_thresh = 0.003;% %% use 0.0005 for 18165\2019-12-30 pars.cellind 2; 6 subthresh, 4 superthresh trials

pars.timecourse_line_width = 1.8; %%%% line for response timecourse
pars.ebar_linewidth = 0.8; 
pars.ebar_cap_size = 3; % errorbar cap length
pars.axes_line_width = 2;
pars.axis_font_bold = 0; %%%
pars.axis_font_size = 15;
pars.font_name = 'Arial'; 
pars.fig_outerposition = [0.05 0.05 0.27 0.35]; % [left bottom width height]
pars.fig_crop_distances = [0 0 0.005 0]; % [left bottom width height]
pars.show_legend = 0; 
pars.line_colors = {[0 0 0], [1 0 0]}; % subthresh, superthresh
% pars.xtick = [-4:2:8]; 
    pars.xtick = [-4:4:8]; 
pars.ytick = [0:0.4:0.8]; 
pars.xlim = [-4 8]; 
pars.ylim = [-0.1, 0.9];
pars.save_timecourse_fig = 1; 
    pars.save_prepend = ['F:\thesis\figs\fig 2 -- timecourse_']; 

% specify which cell to plot
pars.subject = 18165; 
pars.session_date = '2018-12-30';
pars.cellind = 2; % % 18165\2019-12-30 high loc modulation; xy prereg = 73 13

% load(['F:',filesep,num2str(pars.subject),filesep,pars.session_date,filesep,'tuningdat'])

close all

%%% other possible cells to use: % row 22 18164\2019-01-16, xy prereg = 34, 283; 18165 12/30 row 7, 18164 12/21 row 5; 18164 1/171 rows 11, 21
% % % % % extract SF trials data for this cell for its best tested sfreq
thiscell_sfbest_trials = tuningdat.tuning.sf_trials{pars.cellind}(tuningdat.tuning.sf_trials{pars.cellind}.sf==tuningdat.tuning.sf_best_tested_val(pars.cellind),:);
sfbest_trials_loc = tuningdat.tuningpars.locm_trials.sf(tuningdat.tuning.sf_trials{pars.cellind}.sf==tuningdat.tuning.sf_best_tested_val(pars.cellind),:);
superthresh_trials = find(sfbest_trials_loc.locm_forw_mps > locm_thresh);
subthresh_trials = find(sfbest_trials_loc.locm_forw_mps <= locm_thresh);
thissf_stimpar_sets = tuningdat.tuningpars.stimpar_sets(tuningdat.tuningpars.stimpar_sets.sfreq==tuningdat.tuning.sf_best_tested_val(pars.cellind),:);
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
        thissf_stimpar_sets.timecourse_prestim(itrial,:) = num2cell(scope_events.dff(sevents_in_range_prestim,pars.cellind), 1);
        thissf_stimpar_sets.timecourse_stim(itrial,:) = num2cell(scope_events.dff(sevents_in_range_stim_fixed_duration,pars.cellind), 1);
        thissf_stimpar_sets.timecourse_poststim(itrial,:) = num2cell(scope_events.dff(sevents_in_range_poststim,pars.cellind), 1);
end
n_timepoints = length([thissf_stimpar_sets.timecourse_prestim{1}; thissf_stimpar_sets.timecourse_prestim{1}; thissf_stimpar_sets.timecourse_prestim{1}]); 
timecourses_subthresh = table;     timecourses_superthresh = table;
timecourses_subthresh.dff(1:n_timepoints,:) = [cell2mat(thissf_stimpar_sets.timecourse_prestim(subthresh_trials)')', cell2mat(thissf_stimpar_sets.timecourse_stim(subthresh_trials)')',...
    cell2mat(thissf_stimpar_sets.timecourse_poststim(subthresh_trials)')']'; % get dff from all trials
timecourses_subthresh.dff_mean = mean(timecourses_subthresh.dff,2); 
timecourses_subthresh.dff_sem = std(timecourses_subthresh.dff,0,2) ./ size(timecourses_subthresh.dff,2); 
timecourses_superthresh.dff(1:n_timepoints,:) = [cell2mat(thissf_stimpar_sets.timecourse_prestim(superthresh_trials)')', cell2mat(thissf_stimpar_sets.timecourse_stim(superthresh_trials)')',...
    cell2mat(thissf_stimpar_sets.timecourse_poststim(superthresh_trials)')']'; % get dff from all trials
timecourses_superthresh.dff_mean = mean(timecourses_superthresh.dff,2); 
timecourses_superthresh.dff_sem = std(timecourses_superthresh.dff,0,2) ./ size(timecourses_superthresh.dff,2); 
n_trials_lowloc = size(timecourses_subthresh.dff,1); 
n_trials_highloc = size(timecourses_superthresh.dff,1); 
timecourse_x_axis = ( 0:scope_sampling_interval_sec:scope_sampling_interval_sec*[length(timecourses_subthresh.dff_mean)-1] ) - tuningdat.stimpars.isi;
%%% shift x axis back slightly for 18165\2019-12-30 pars.cellind 2... onset and offset appear to come early (maybe due to cropped dff timepoints at end of timecourse segmeents)
timecourse_x_axis = timecourse_x_axis + scope_sampling_interval_sec;

timecourse_locomotion_fig = figure('units','normalized','outerposition',pars.fig_outerposition);
subthresh_plot = plot(timecourse_x_axis,timecourses_subthresh.dff_mean,'Color',pars.line_colors{1});
hold on
subthresh_plot.LineWidth = pars.timecourse_line_width;
subthresh_ebar = errorbar(timecourse_x_axis, timecourses_subthresh.dff_mean, timecourses_subthresh.dff_sem,'Color',pars.line_colors{1}); 
subthresh_ebar.LineWidth = pars.ebar_linewidth; 
subthresh_ebar.CapSize = pars.ebar_cap_size; 

superthresh_plot = plot(timecourse_x_axis,timecourses_superthresh.dff_mean,'Color',pars.line_colors{2});
hold on
superthresh_plot.LineWidth = pars.timecourse_line_width;
superthresh_ebar = errorbar(timecourse_x_axis, timecourses_superthresh.dff_mean, timecourses_superthresh.dff_sem,'Color',pars.line_colors{2}); 
superthresh_ebar.LineWidth = pars.ebar_linewidth; 
superthresh_ebar.CapSize = pars.ebar_cap_size; 

% 
hLeg=findobj(timecourse_locomotion_fig,'type','legend');
if ~pars.show_legend 
    set(hLeg,'visible','off');     
elseif pars.show_legend 
    legend({['Low Locomotion (N = ',num2str(n_trials_lowloc),' trials)'],  ['High Locomotion (N = ',num2str(n_trials_highloc),' trials)']})
end
hax = gca; 
ylim(pars.ylim); 
set(gca,'Box','off')
set(gca,'LineWidth',pars.axes_line_width)
set(gca,'FontSize',pars.axis_font_size)
if pars.axis_font_bold; set(gca,'FontWeight','bold'); end
set(gca,'LooseInset',get(gca,'TightInset')+pars.fig_crop_distances) % crop borders
hax.FontName = pars.font_name;
hax.XTick = pars.xtick; 
hax.XLim = pars.xlim; 
if isfield(pars,'ytick'); hax.YTick = pars.ytick; end
if pars.save_timecourse_fig
    savename = [pars.save_prepend, num2str(pars.subject), '_', pars.session_date, '_cell', num2str(pars.cellind)]; 
    saveas(timecourse_locomotion_fig, savename, 'svg')
end


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








