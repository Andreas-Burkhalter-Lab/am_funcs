%%%%% get_f_behavior_stim: align 2-photon scope frames with simultaneously collected locomotion
%%%%% data and visual stimulus data acquired in the Ed Han lab
% 
%      res = get_f_behavior_stim(dff_file, abf_file, stimdata_file, scopetiming_file, regops_file)
%

% Instructions:
%       1. Run new_main.m on '*_plane1.mat','*_plane2.mat' (etc.) files to manually  
%       select neurons/neurites you want to use and save the results, producing
%       '_proc.mat' files.
%       2. Run convertSuite2Pdata_AM(1) on your '_proc.mat' files (may take a
%       few minutes to process).
%       3. Run get_f_behavior_stim (this function), either with the specified inputs or using the
%       input prompts. Your results will be saved in a table containing dF/F
%       traces along with locomotion speed and visual stimulus parameters.
%
%%% Inputs: 
% 1. dff_file = dF/F data for each plane generated by convertSuite2Pdata_AM.m (not _proc files from new_main)   
% 2. abf_file = behavioral + microscope timing file
% 3. stimdata_file = data on stimulus parameters for each trial
% 4. scopetiming_file = file listing start and stop points of scope within abf file generated with getAbfStartStop.m   
% 5. regops_file = file generated during registration
% 6. pupil_recording_file = file with data from pupil video
%
%%% Outputs:
% results struct containing the following fields:
%   scope_events: table with dF/F, stimulus, and locomotion data for this plane   
%       -in 'F' variable, columns are cells/neurites and rows are time points 
%   iplane: imaging plane index (lower = deeper under pia)
%   nplanes_session_total: number of planes acquired this session
%   meanImage: time-averaged image of this plane
%   stimpar_sets: responses to each stimulus event
%   stimpars: general stimulus info about this session 

%%%%% last updated 2018/11/27 by Andrew Meier on thermaltake
function res = get_f_behavior_stim_func(dff_file, triggertiming_file, stimdata_file, scopetiming_file, regops_file, pupil_recording_file)

save_results = 0;
    save_append = '_scope_events';
metersPerSec_per_ballVolt = 0.49; % han lab recordings only: meters per second of locomotion per V in forward and yaw  abf channels (measured by Moi) 

load_f_behavior_stim_data()

% get scope events from trigger timing file
scope_events_allplanes = find_timingchan_events(abf_timepoints(:,scopetiming.scope_chan)');
% delete events found outside the range specified by scope timing file
deleteevents = scope_events_allplanes.onset < scopetiming.scope_start_timepoint | nansum([scope_events_allplanes.onset,scope_events_allplanes.duration],2) > scopetiming.scope_stop_timepoint;
scope_events_allplanes(deleteevents,:) = [];
% get scope events corresponding to this plane; frames are acquired in order of plane name (iplane) 
if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
    if ops1{1}.dobidi % bidirectional
        fprintf('Processing scope timing for bidirectional acquisition...\n')
        scope_events = scope_events_allplanes(iplane:nplanes_in_session_total:end,:); 
    elseif ~ops1{1}.dobidi % unidirectional
        go_on = input('Unidirectional acquisition - are you sure? (y=yes) ','s');
        if ~strcmp(go_on,'y')
            error('quitting this function')
        end
        scope_events = scope_events_allplanes(iplane:2*nplanes_in_session_total:end,:);  % skip every other scope event because filler scope events are added for unidirectional acquisition
    end
elseif strcmp(stimpars.computer,'')  %% if recorded on 4th floor scope
    frpintf('Assuming single plane acquisition because recording is from 4th floor scope...\n')
    scope_events = scope_events_allplanes; 
end
ntrials = height(stimpar_sets);

%% find stim events
scope_events = scope_events(1:nimages,:); % ignore scope events that don't have a corresponding image
nevents = height(scope_events);
scope_events.locomotion_forw_mpersec = NaN(nevents,1);
scope_events.locomotion_yaw_mpersec = NaN(nevents,1);
scope_events.stim_present = false(nevents,1);
if ~isempty(pupil_recording_file) % if pupil recording data was provided
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
scope_events = [scope_events table(dff_data.dff,'VariableNames',{'dff'})];

% get locomotion, stim, and pupil data for each scope event for this plane
tpoints_frame_duration = floor(mean(diff(scope_events_allplanes.onset))); %tpoints between scope frames
for ievent = 1:nevents
    thisevent_tpoints = scope_events.onset(ievent) : scope_events.onset(ievent) + tpoints_frame_duration;  % time points covered by this scope events
    % locomotion
    loc_forw_v = mean(abf_timepoints(thisevent_tpoints,locm_forw_chan)); % take mean of locomotion over this frame
    loc_yaw_v = mean(abf_timepoints(thisevent_tpoints,locm_yaw_chan)); % take mean of locomotion over this frame
    scope_events.locomotion_forw_mpersec(ievent) = loc_forw_v * metersPerSec_per_ballVolt;
    scope_events.locomotion_yaw_mpersec(ievent) = loc_yaw_v * metersPerSec_per_ballVolt;
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
    if ~isempty(pupil_recording_file) % if pupil recording data was provided
        [timedif, rowmatch] = min(abs(scope_events.onset(ievent) - pupil.onset)); % rowmatch = row within pupil table of closest pupil frame to this scope event
        if timedif < timedif_lim % if pupil frame and scope events are close enough in time
            scope_events.pupil_centeryx(ievent,:) = pupil.centeryx(rowmatch,:);
            scope_events.pupil_keepframe(ievent) = pupil.keepframe(rowmatch);
        end
    end
end
        

% first frame collected for plane 1 will often be from the wrong plane, so
% delete it (see first few frames of raw .sbx files for examples)
if iplane == 1
    scope_events = scope_events(2:end,:);
end

%%%% han lab recordings only: if badframes file is present, cut out from scope_events all frames that occurred before the last or after the first bad frame
res = struct; % results struct
if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
    [res, scope_events] = cutBadFramesBlock(res, scope_events, scopetiming);
end


%% output and save results
scope_events.chan_val = [];
scope_events.duration = [];
scope_events.Properties.VariableNames{'onset'} = 'onset_ms';
res.scope_events = scope_events;
res.stimpars = stimpars;
res.iplane = iplane;
res.nplanes_in_session_total = nplanes_in_session_total;
res.stimpar_sets = stimpar_sets;
res.meanImage = dff_data.meanImage;

if save_results
    savename = [getfname(dff_file), save_append];
    save(savename,'res')
end