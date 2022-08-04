% pakan_locomotion_modulation
%      get the locomotion modulation index described in Pakan et al. 2016 - Behavioral-state modulation of inhibition is context-dependent and cell type specific in mouse visual cortex
%
% tuning_out = pakan_locomotion_modulation(tuning_in, scope_events, scopetiming, pars)
%
% locomotion moduation index = (R_L – R_s)/(R_L + R_s), where R_L = df/f during locomotion, R_s=df/f during rest
%
%%%  called by get_tuning_curves
%
%%%%% updated 2020/9/21 on thermaltake


function tuning_out = pakan_locomotion_modulation(tuning_in, scope_events, scopetiming, pars)

%%% threshold speed in meters/s for considering locomotion to be occurring (.001 in Pakan 2016)
pars = vardefault('pars',struct);
pars.locomotion_threshold_mps = field_default(pars, 'locomotion_threshold_mps', 0.001); 

% only analyze scope events within the limits specified in scopetiming
scope_events_cut = scope_events(scope_events.onset > scopetiming.scope_start_timepoint & scope_events.onset < scopetiming.scope_stop_timepoint, :);

% sort scope events by whether there was locomotion during them
scope_events_cut.locomotion_on = scope_events_cut.locomotion_forw_mpersec > pars.locomotion_threshold_mps;


% get modulation index during stim
dff_stim_and_locomotion = scope_events_cut.dff(scope_events_cut.stim_present & scope_events_cut.locomotion_on,:);
dff_stim_and_rest = scope_events_cut.dff(scope_events_cut.stim_present & ~scope_events_cut.locomotion_on,:);
dff_loc_mean = mean(dff_stim_and_locomotion); % average over relevant scope events
dff_rest_mean = mean(dff_stim_and_rest); 
[loc_sgnf, loc_pval] = ttest2(dff_stim_and_locomotion, dff_stim_and_rest); %% ttest of differnce between rest and loc
tuning_out = tuning_in;
tuning_out.pakan_loc_index = ( [dff_loc_mean - dff_rest_mean] ./ [dff_loc_mean + dff_rest_mean] )';
tuning_out.pakan_loc_pval = loc_pval';
tuning_out.pakan_loc_sgnf = loc_sgnf';
tuning_out = movevars(tuning_out,{'pakan_loc_pval','pakan_loc_sgnf'},'After','pakan_loc_index'); % rearrange in case these vars were analyzed separately