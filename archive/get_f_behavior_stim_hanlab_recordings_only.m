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

%%%%% last updated 2018/09/27 by Andrew Meier on thermaltake
function res = get_f_behavior_stim(dff_file, abf_file, stimdata_file, scopetiming_file, regops_file)

save_results = 0;
    save_append = '_scope_events';

% meters per second of locomotion per V in forward and yaw  abf channels (measured by Moi) 
metersPerSec_per_ballVolt = 0.49;

load_f_behavior_stim_data()

% get scope events from abf file
scope_events_allplanes = find_stimchan_events(abf_timepoints(:,scopetiming.scope_chan)');
% delete events found outside the range specified by scope timing file
deleteevents = scope_events_allplanes.onset < scopetiming.scope_start_timepoint | nansum([scope_events_allplanes.onset,scope_events_allplanes.duration],2) > scopetiming.scope_stop_timepoint;
scope_events_allplanes(deleteevents,:) = [];
% get scope events corresponding to this plane; frames are acquired in order of plane name (iplane)   
scope_events = scope_events_allplanes(iplane:nplanes_in_session_total:end,:); 

% get scope events corresponding to this plane; frames are acquired in order of plane name (iplane) 
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

%% find stim events
scope_events = scope_events(1:nimages,:); % ignore scope events that don't have a corresponding image
scope_events.locomotion_forw_mpersec = NaN(height(scope_events),1);
scope_events.locomotion_yaw_mpersec = NaN(height(scope_events),1);
scope_events.stim_present = false(height(scope_events),1);
scope_events.stim_Angle_deg = NaN(height(scope_events),1);
scope_events.stim_sf_cyclesperdeg = NaN(height(scope_events),1);
scope_events.stim_tf_hz = NaN(height(scope_events),1);
scope_events.stim_diam_deg = NaN(height(scope_events),1);
scope_events = [scope_events table(dff_data.dff,'VariableNames',{'dff'})];

% get locomotion and stim data for each scope event for this plane
tpoints_frame_duration = floor(mean(diff(scope_events_allplanes.onset))); %tpoints between scope frames
for ievent = 1:height(scope_events)
    thisevent_tpoints = scope_events.onset(ievent) : scope_events.onset(ievent) + tpoints_frame_duration; 
    loc_forw_v = mean(abf_timepoints(thisevent_tpoints,locm_forw_chan));
    loc_yaw_v = mean(abf_timepoints(thisevent_tpoints,locm_yaw_chan));
    scope_events.locomotion_forw_mpersec(ievent) = loc_forw_v * metersPerSec_per_ballVolt;
    scope_events.locomotion_yaw_mpersec(ievent) = loc_yaw_v * metersPerSec_per_ballVolt;
    coinciding_stim_ind = find(par_sets.stim_onset < scope_events.onset(ievent)  &  par_sets.stim_onset+par_sets.stim_duration > scope_events.onset(ievent)); 
    if ~isempty(coinciding_stim_ind)
        scope_events.stim_present(ievent) = true;
%         scope_events.stim_Angle_deg(ievent) = par_sets.Angle(coinciding_stim_ind);
        scope_events.stim_orient_deg(ievent) = par_sets.orient(coinciding_stim_ind);
        scope_events.stim_sf_cyclesperdeg(ievent) = par_sets.sfreq(coinciding_stim_ind);
        scope_events.stim_tf_hz(ievent) = par_sets.tfreq(coinciding_stim_ind);
        scope_events.stim_diam_deg(ievent) = par_sets.diam(coinciding_stim_ind);
    end
end
        

% first frame collected for plane 1 will often be from the wrong plane, so
% delete it (see first few frames of raw .sbx files for examples)
if iplane == 1
    scope_events = scope_events(2:end,:);
end

%%%% if badframes file is present, cut out from scope_events all frames that occurred before the last or after the first bad frame
res = struct; % results struct
[res, scope_events] = cutBadFramesBlock(res, scope_events, scopetiming);

%% output and save results
scope_events.chan_val = [];
scope_events.duration = [];
scope_events.Properties.VariableNames{'onset'} = 'onset_ms';
res.scope_events = scope_events;
res.stimpars = pars;
res.iplane = iplane;
res.nplanes_in_session_total = nplanes_in_session_total;
res.stimpar_sets = par_sets;
res.meanImage = dff_data.meanImage;

if save_results
    savename = [getfname(dff_file), save_append];
    save(savename,'res')
end