%%% get locomotion correlation for each roi in a table created by get_tuning_curves
% analyze locomotion correlations only from scope images from the period specified in the 'scopetiming_dark' file
%%% updated 2019-2-1 on thermaltake

function tuningdat_out = locomotion_analysis(tuningdat_in)

tuningdat_out = tuningdat_in; 
load(tuningdat_out.input_files.scopetiming_file_dark,'scopetiming');
tuningdat_out.scopetiming_dark = scopetiming;
nrois = size(tuningdat_out.scope_events.dff,2);
tuningdat_out.tuning.locm_cor_r = NaN(nrois,1); % correlation coefficient of dff with forward locomotion
tuningdat_out.tuning.locm_cor_p = NaN(nrois,1); 

events_to_analyze = tuningdat_out.scope_events.onset > scopetiming.scope_start_timepoint & tuningdat_out.scope_events.onset < scopetiming.scope_stop_timepoint; 
% note - the check below does not look at the whole trigdata timecourse; it can only find overlap in the specified periods designed as rf test, dark, and tuning testing
if any(tuningdat_out.scope_events.stim_present(events_to_analyze)) % if we find stim events during the specified dark period
    error('stim events found during specified dark period')
end
locm_timecourse = tuningdat_out.scope_events.locomotion_forw_mpersec(events_to_analyze); % locomotion forward during the dark period
  
for iroi = 1:nrois
    [loc_r, loc_p] = corrcoef(locm_timecourse, tuningdat_out.scope_events.dff(events_to_analyze,iroi));
    tuningdat_out.tuning.locm_cor_r(iroi) = loc_r(2,1);
    tuningdat_out.tuning.locm_cor_p(iroi) = loc_p(2,1);
end