%%%%% get_f_behavior_stim: align 2-photon scope frames with simultaneously collected locomotion
%%%%%   data and visual stimulus data
%%%%% last updated 2020/02/04  on thermaltake

save_results = 0;
    save_append = '_scope_events';
metersPerSec_per_ballVolt = 0.49; % han lab recordings only: meters per second of locomotion per V in forward and yaw  abf channels (measured by Moi) 

% get scope events corresponding to this plane; frames are acquired in order of plane name (iplane) 
    %%%% note: if multiple subsessions were concatenated, this method of picking scope events matched to iplane is wrong if the # of images acquired 
    %%%%    in each sub-session within the concatenated run is not divisible by nplanes; scope event picking should start over at the beginning of each sub-session
if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
    if ops1{1}.dobidi % bidirectional
        fprintf('Processing scope timing for bidirectional acquisition...\n')
        scope_events = scope_events_allplanes(iplane:nplanes_in_session_total:end,:); 
    elseif ~ops1{1}.dobidi % unidirectional
        if ~pars.skip_all_checks
            go_on = input('Unidirectional acquisition - are you sure? (y=yes) ','s');
            if ~strcmp(go_on,'y')
                error('quitting this function')
            end
        end
        scope_events = scope_events_allplanes(iplane:2*nplanes_in_session_total:end,:);  % skip every other scope event because filler scope events are added for unidirectional acquisition
    end
end
if strcmp(stimpars.computer,'ANDREWLAB-PC') %% if from 4th floor scope, assume unidirectional acquisition
    %%% in case there are leftover scope pulses at the end of a preceding session (without matching scope images) from multi-plane acquisition, ...
    %%%     ... load all session timing data and eliminate extra pulses so that scope image timing in this session is not incorrectly offset
    if ~isempty(pars.all_scopetiming_files) && nplanes_in_session_total > 1 % if preceding sessions are listed and multiple planes were acquired
        n_events_within_a_session = 0; % add up event rows within each session, warn if there are events found outside any session
        
        %% section below attempts to correct for extra pulses not matching frames acquired during multi plane acquisition...
        %%%% .... but incorrectly assumes planes are always acquired as a complete stack...
        %%%% ... appears to insert buffer planes even when there is no mismatch between pulses and planes
        %%%%%% .... in future, maybe just manually check that there isn't a mismatch
        
        for isession = 1:length(pars.all_scopetiming_files)
            scopetiming_for_pulse_correction = load(pars.all_scopetiming_files{isession},'scopetiming'); scopetiming_for_pulse_correction = scopetiming_for_pulse_correction.scopetiming; 
            scope_event_rows_this_session =... %%% scope event rows matching the session being checked
                find(scope_events_allplanes.onset > scopetiming_for_pulse_correction.scope_start_timepoint & scope_events_allplanes.onset < scopetiming_for_pulse_correction.scope_stop_timepoint);
             %%% # scope events per sub-session must be divisible by number of planes because planes are always acquired as a complete stack...
             %%%        ... assign a complete stack to remainder pulses at the end of the session, do not delete remainder pulses
            n_remainder_events = mod(length(scope_event_rows_this_session), nplanes_in_session_total); % events not divisible by n_planes_in_session_total
            n_events_to_add = nplanes_in_session_total - n_remainder_events; %%% number of scope frames whose pulses were cut off by end of session
            mean_event_spacing = mean(diff(scope_events_allplanes.onset(scope_event_rows_this_session))); % scope event duration in this session
            last_event_table = scope_events_allplanes(scope_event_rows_this_session(end),:); 
            buffer_events_table = repmat(last_event_table,n_events_to_add,1); %%% events to insert to match scope frames whose pulses were cut off by end of session
            buffer_events_table.onset = last_event_table.onset + [floor(mean_event_spacing:mean_event_spacing:n_events_to_add*mean_event_spacing)]'; % assign estimated timing for unmatched scope frames
            % add the buffer events for unmatched scope frames
            scope_events_allplanes = [scope_events_allplanes(1:scope_event_rows_this_session(end),:); buffer_events_table; scope_events_allplanes(scope_event_rows_this_session(end)+1:end,:)];
            n_events_within_a_session = n_events_within_a_session + length(scope_event_rows_this_session) + n_events_to_add;
        end
        
        %% 
        if height(scope_events_allplanes) > n_events_within_a_session %%% check that all scope events are accounted for by scopetiming files
            error('Found scope events outside of any indicated recording session; these events may cause timing offsets if they''re not equal in number to number of planes')
        end
    end
    scope_events = scope_events_allplanes(iplane:nplanes_in_session_total:end,:);
end

% downsample if applicable
if pars.dff_downsample_factor == 1
    pars.dff_down_sample_factor = []; 
else
   fprintf(['Downsampling by keeping only 1 out of every ' num2str(pars.dff_downsample_factor) ' scope frames.\n'])
   old_scope_events_height = height(scope_events);
   scope_events = scope_events(1:pars.dff_downsample_factor:end,:); % downsample scope events 
   dff_data.dff = dff_data.dff(1:pars.dff_downsample_factor:end,:); % downsample dff timepoints
   nimages = size(dff_data.dff,1); % recount number of dff values
   stimdur_scope_events = stimdur_scope_events / pars.dff_down_sample_factor; % there are now fewer scope events per stim duration
   isi_scope_events = isi_scope_events / pars.dff_down_sample_factor; % there are now fewer scope events per isi
end

ntrials = height(stimpar_sets);

%% align scope and dff data
%   -to start, ignore manually picked scopetiming limits, just check whether n scope pulses in the entire timecourse == n total dff values provided
%   **** align pulses with dff values from the beginning of the WHOLE TIMECOURSE, NOT from the beginning of the epoch specified by manual scopetiming limits
%      * assume that first pulse == first dff value, and non-aligned extra pulses or dff values occur only at the end
%      * cut pulses and dff events from the whole timecourse based on manual scopetiming limits after both have been combined
%%%% analyze the whole timecourse because multiple copies of it (extracting different variables) will later be combined in full_sesion_analysis.m
nscopeevents_orig = height(scope_events);
if nimages < nscopeevents_orig % if scope images are missing
    if ~pars.skip_all_checks
        go_on = input(['Found fewer scope images (' num2str(nimages) ') than pulses (' num2str(nscopeevents_orig) ') in scope channel - ignore extra pulses? (y=yes) '],'s');
        if ~strcmp(go_on,'y')
            error('quitting this function')
        end
    end
    scope_events = scope_events(1:nimages,:); % ignore scope events that don't have a corresponding image; assuming images were cut off from end, not beginning of recording
elseif nimages > nscopeevents_orig % if scope trigger events are missing
    if ~pars.skip_all_checks
        go_on = input(['Found more scope images (' num2str(nimages) ') than pulses (' num2str(nscopeevents_orig) ') in scope channel - ignore extra images? (y=yes) '],'s');
        if ~strcmp(go_on,'y')
            error('quitting this function')
        end
    end
    dff_data.dff = dff_data.dff(1:height(scope_events),:); % ignore scope images without a corresponding trigger event; assuming events were cut off from end, not beginning of recording
end

nevents = height(scope_events);
scope_events.locomotion_forw_mpersec = NaN(nevents,1);
scope_events.locomotion_yaw_mpersec = NaN(nevents,1);
scope_events.stim_present = false(nevents,1);
if ~isempty(pupil_data_file) % if pupil recording data was provided
    scope_events.pupil_centeryx = NaN(nevents,2);
    scope_events.pupil_keepframe = false(nevents,1); % whether or not to analyze the frame based on pupil data
    timedif_lim = 1.2*max(diff(scope_events.onset)); % scope event and pupil frame must be within this time difference to match them
end
switch stimpars.experiment_type
    case 'sf_tf_orient_diam'
        scope_events.stim_orient_deg = NaN(height(scope_events),1);
        scope_events.stim_sf_cyclesperdeg = NaN(height(scope_events),1);
        scope_events.stim_tf_hz = NaN(height(scope_events),1);
        scope_events.stim_diam_deg = NaN(height(scope_events),1);
    case 'rf_mapping'
        scope_events.stim_row = NaN(height(scope_events),1);
        scope_events.stim_column = NaN(height(scope_events),1);
end
scope_events = [scope_events table(dff_data.dff,'VariableNames',{'dff'})]; % align scope pulses with dff values - see above for aligment assumptions
% delete scope evets AND corresponding dff values found outside the range specified by scope timing file
deleteevents = scope_events_allplanes.onset < scopetiming.scope_start_timepoint | nansum([scope_events_allplanes.onset,scope_events_allplanes.duration],2) > scopetiming.scope_stop_timepoint;
scope_events_allplanes(deleteevents,:) = [];
seconds_per_scope_event = range(scope_events_allplanes.onset) * trig_sampling_interval_us/1e6 /  [height(scope_events_allplanes)-1]; 
stimdur_scope_events = floor(stimpars.stimdur / seconds_per_scope_event); % scope events per stimulus duration
res.stimdur_scope_events = stimdur_scope_events; 
isi_scope_events = floor(stimpars.isi / seconds_per_scope_event); % scope events per interstimulus interval
res.isi_scope_events = isi_scope_events; 

%% get locomotion, stim, and pupil data for each scope event for this plane
tpoints_frame_duration = floor(mean(diff(scope_events_allplanes.onset))); %tpoints between scope frames
for ievent = 1:nevents
    end_tpoint = min([scope_events.onset(ievent) + tpoints_frame_duration, height(trigdata)]); % take the last time point if it occurs before a full trial duration has elapsed
    thisevent_tpoints = scope_events.onset(ievent) : end_tpoint;  % time points covered by these scope events
    % locomotion
    if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
        loc_forw_v = mean(abf_timepoints(thisevent_tpoints,locm_forw_chan)); % take mean of locomotion over this frame
        loc_yaw_v = mean(abf_timepoints(thisevent_tpoints,locm_yaw_chan)); % take mean of locomotion over this frame
        scope_events.locomotion_forw_mpersec(ievent) = loc_forw_v * metersPerSec_per_ballVolt;
        scope_events.locomotion_yaw_mpersec(ievent) = loc_yaw_v * metersPerSec_per_ballVolt;
    elseif strcmp(stimpars.computer,'ANDREWLAB-PC')  %% if recorded on 4th floor scope
        scope_events.locomotion_forw_mpersec(ievent) = mean(trigdata.locm_forw_mps(thisevent_tpoints)); % get mean forward locomotion over this timeframe; yaw will be NaN
    end
    % stimuli
    coinciding_stim_ind = find(stimpar_sets.stim_onset < scope_events.onset(ievent)  &  stimpar_sets.stim_onset+stimpar_sets.stim_duration > scope_events.onset(ievent)); 
    if ~isempty(coinciding_stim_ind)
        scope_events.stim_present(ievent) = true;
        switch stimpars.experiment_type
            case 'sf_tf_orient_diam'
                scope_events.stim_orient_deg(ievent) = stimpar_sets.orient(coinciding_stim_ind);
                scope_events.stim_sf_cyclesperdeg(ievent) = stimpar_sets.sfreq(coinciding_stim_ind);
                scope_events.stim_tf_hz(ievent) = stimpar_sets.tfreq(coinciding_stim_ind);
                scope_events.stim_diam_deg(ievent) = stimpar_sets.diam(coinciding_stim_ind);
            case 'rf_mapping'
                scope_events.stim_row(ievent) = stimpar_sets.row(coinciding_stim_ind);
                scope_events.stim_column(ievent) = stimpar_sets.column(coinciding_stim_ind);
        end
    end
    % pupil data
    if ~isempty(pupil_data_file) % if pupil recording data was provided
        [timedif, rowmatch] = min(abs(scope_events.onset(ievent) - pupil.onset)); % rowmatch = row within pupil table of closest pupil frame to this scope event
        if timedif < timedif_lim  % if pupil frame and the nearest scope event are close enough in time
            scope_events.pupil_centeryx(ievent,:) = pupil.centeryx(rowmatch,:);
            scope_events.pupil_keepframe(ievent) = pupil.keepframe(rowmatch);
        end
    end
end
        

% first frame collected for plane 1 will often be from the wrong plane in han lab recordings, so
% delete it (see first few frames of raw .sbx files for examples)
if strcmp(stimpars.computer,'IMAGING_VR-PC') && iplane == 1 && ops1{1}.nplanes>1 && ~deleteevents(1) % only if we have multiple planes and didn't already delete the first event
    scope_events = scope_events(2:end,:);
end

% % % % % %%%% han lab recordings only (comment in to use): if badframes file is present, cut out from scope_events all frames that occurred before the last or after the first bad frame
% % % % %  %%%%% generally do not use this section if we are interested in trial responses (which might get cut by the following function);...
% % % % %  %%%%%       ...only use this section if we are interested in the entire activity timecourse and not trial responses
% % % % % badframes_data = struct; % results struct
% % % % % if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
% % % % %     [badframes_data, scope_events] = cutBadFramesBlock(badframes_data, scope_events, scopetiming);
% % % % % end

%% output and save results
scope_events.chan_val = [];
scope_events.duration = [];

if save_results
    savename = [getfname(dff_file), save_append];
    save(savename,'scope_events')
end