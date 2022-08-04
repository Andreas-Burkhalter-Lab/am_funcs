%%% timecourse_response_plotting
% plot timecourse of response to best stimpar value, from tuning table saved by get_tuning_curves
%%% updated 19-6-03
 
function timecourse_response_plotting(roitable,tablerow,stimparname)
 
timecourse_color = [0 0 0];
onset_offset_color = [.5 .5 1];
onset_offset_linestyle = '--';
stim_duration_sec = 4; % seconds

timecourse_mean_varname = [stimparname, 'timecourse_mean_best']; 
timecourse_sem_varname = [stimparname, 'timecourse_sem_best']; 

%%% seconds per timestep - may need to adjust based on which recording is being plotted
timestep_sec = stim_duration_sec / length(roitable{tablerow,timecourse_mean_varname}{1}); % assume that timestep is stim duration divided by n timesteps during stim period

timecourse_mean = cell2mat(roitable{tablerow,timecourse_mean_varname});
timecourse_sem = cell2mat(roitable{tablerow,timecourse_sem_varname});
xticks = linspace(0, timestep_sec*length(timecourse_mean), length(timecourse_mean)); 
stim_onset = timestep_sec * length(roitable{tablerow,timecourse_mean_varname}{1}); % for marking stim onset
stim_offset = timestep_sec * [length(roitable{tablerow,timecourse_mean_varname}{1}) + length(roitable{tablerow,timecourse_mean_varname}{2})]; % for marking stim offset
 
tcourse_plot = errorbar(xticks, timecourse_mean, timecourse_sem);
set(tcourse_plot,'Color',timecourse_color)
xlabel('Time (s)')
ylabel('dF/F')
ylimits = get(gca,'Ylim'); 
hold on
onset_plot = plot([stim_onset, stim_onset], ylimits, 'LineStyle', onset_offset_linestyle, 'Color', onset_offset_color);
offset_plot = plot([stim_offset, stim_offset], ylimits, 'LineStyle', onset_offset_linestyle, 'Color', onset_offset_color);
ylim(ylimits)
hold off