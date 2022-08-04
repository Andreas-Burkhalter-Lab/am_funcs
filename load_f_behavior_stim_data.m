%%%%%%%%%% load data from files for tuning curve processing
% part of get_tuning_curves
% updated 2020/4/21 on thermaltake

prairie_stim_chan_name = 'Input2'; % name of input channel for visual stim timing
        prairie_trigger_sampling_rate_hz = 10000;     %% this value is not necessarily correct for all recordings... need to find out how to get sampling rate from trigger file
han_lab_stim_chan = 'Psycho_In'; % name of input channel for visual stim timing
han_lab_scope_pupil_chan = 'Lick'; % name of input channel for pupil recording acqusition timing
    pupil_timing_pars.pulse.Vstimon = 5; % V turns to 5 at beginning of frame acquisition
    pupil_timing_pars.pulse.Vstimoff = 0;
    pupil_pulses_per_frame = 2; 
pupil_pulse_duration_range_fraction_limit = 0.25; % this value > range(pupil_frame_events.duration) / nanmean(pupil_frame_events.duration), OK to correct for missed frames
pupil_framepulse_dif_fraction_limit = 0.01; % if abs(frames-pulses)/npulses is less than this value, OK to correct for missed frames

%% load files, set params
if ~exist('dff_file','var') || isempty(dff_file)
    dff_file = uigetfile('*.mat','Select dff .mat file.','MultiSelect','off');
end
dff_data = load(dff_file,'dff','meanImage','masks','fName');
res.meanImage_pre_rotate = dff_data.meanImage;

if ~exist('stimdata_file','var') || isempty(stimdata_file)
    stimdata_file = uigetfile('*stimdata*.mat','Select stimdata .mat file.');
end
load(stimdata_file,'-mat');
stimpars.experiment_type = field_default(stimpars, 'experiment_type', 'sf_tf_orient_diam'); %%% assume it's sf_tf_orient_diam if not specified
% stimpars.stim_center_yx = [stimpars.stim_center(2), stimpars.stim_center(1)]; %%% put in format consistent with rest of analysis
% stimpars = rmfield(stimpars, 'stim_center');

if ~exist('scopetiming_file','var') || isempty(scopetiming_file)
    scopetiming_file = uigetfile('*_scopetiming.mat','Select scopetiming .mat file.');
end
load(scopetiming_file);
res.scopetiming = scopetiming;

if ~exist('regops_file') || isempty(regops_file)
    regops_file = uigetfile('*regops*.mat','Select regops .mat file.');
end
load(regops_file);
res.regops_file = regops_file;

% number of planes
if isfield(ops1{1},'nplanes')
    nplanes_in_session_total = ops1{1}.nplanes;
else 
    nplanes_in_session_total = ops1{1}.numPlanes;
end
res.nplanes_in_session_total = nplanes_in_session_total;

%% load han lab data
if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope
    % trigger timing
    if ~exist('triggertiming_file','var') || isempty(triggertiming_file)
        triggertiming_file = uigetfile('*.abf','Select .abf file.');
    end
    try
        [abf_timepoints, abf_sampling_interval_us, abfinfo] = abfload(triggertiming_file); % abf_sampling_interval_us = sampling interval microseconds
    catch abf_mexception % if we had to convert from abf to mat file
        fprintf('     ABF file not readable as .abf; loading as .mat file instead.\n')
        load(triggertiming_file,'-mat')
    end
    trig_sampling_interval_us = abf_sampling_interval_us;
        % get stim events from abf file
    stimchan = find(strcmp(abfinfo.recChNames,han_lab_stim_chan));
    locm_forw_chan = find(strcmp(abfinfo.recChNames,'Ball_forw'));
    locm_yaw_chan = find(strcmp(abfinfo.recChNames,'Ball_yaw'));
    stim_events = find_timingchan_events(abf_timepoints(:,stimchan)',stimpars,abf_sampling_interval_us);
    
    % get scope events from trigger timing file
    scope_events_allplanes = find_timingchan_events(abf_timepoints(:,scopetiming.scope_chan)');
    
    % rotate mean image so that anterior is upward; anterior should already be upward on 4th floor scope recordings
    res.meanImage_ant_up = rot90(dff_data.meanImage,-1); % rotate so that anterior==up
        
%% load 4th floor scope data
elseif strcmp(stimpars.computer,'ANDREWLAB-PC')  %% if recorded on 4th floor scope
        fprintf(['\n .....Assuming sample rate for voltage recording = ' num2str(prairie_trigger_sampling_rate_hz)])
    if ~exist('triggertiming_file','var') || isempty(triggertiming_file)
        triggertiming_file = uigetfile('*.csv','Select .csv file.');
    end
    trigdata = load_prairie_triggertiming(triggertiming_file);
    trig_sampling_interval_us = 1e6 / prairie_trigger_sampling_rate_hz; % convert from hz to microseconds
    stim_events = find_timingchan_events(trigdata{:,prairie_stim_chan_name} , stimpars); 
    valid_events = stim_events.onset > scopetiming.scope_start_timepoint & stim_events.onset < scopetiming.scope_stop_timepoint; % events specified by scopetiming
    stim_events = stim_events(valid_events,:); % only include events within the bounds defined by scopetiming_file
    % get scope events from trigger timing file
    scope_events_allplanes = find_timingchan_events(trigdata{:,scopetiming.scope_chan_name});
    delete_events = scope_events_allplanes.duration < 500; %%% scope events of this duration or less are probably the extra 'tseries end' pulse and do not denote an extra image acquired
    scope_events_allplanes(delete_events,:) = []; % delete tseries-end scope events
    scope_events = scope_events_allplanes; 
else
    error('Recording computer not recognized')
end

%% pupil data
% load pupil frame data, get timing
pupil_data_file = vardefault('pupil_data_file',[]);
if ~isempty(pupil_data_file) %%% if a pupil recording file was specified
    load(pupil_data_file) % should contain table 'pupil' and struct 'pupilinfo'
    if strcmp(stimpars.computer,'IMAGING_VR-PC') %% if recorded on han lab scope    
        pupilchan =  find(strcmp(abfinfo.recChNames,han_lab_scope_pupil_chan));
        pupil_frame_events = find_timingchan_events(abf_timepoints(:,pupilchan)', pupil_timing_pars);
        pupil_frame_events = pupil_frame_events(1:pupil_pulses_per_frame:end,:); % downsample so that there is one pulse per frame
        n_pupil_frames_this_session = pupilinfo.nframes; % assume all pupil frames were from this session (pupil data is not from concatenated sessions)
    elseif strcmp(stimpars.computer,'ANDREWLAB-PC')  %% if recorded on 4th floor scope 
        % get the session index being analyzed by comparing scopetiming filename to the pupil recording filenames
        session_index = find(strcmp(cellfun(@getfname,pupilinfo.scopetiming_filename,'UniformOutput',0), getfname(scopetiming_file)));
        pupil_frames_before_this_session = sum(pupilinfo.nframes_raw(1:session_index-1)); %%% =0 if session_index==1
        n_pupil_frames_this_session = pupilinfo.nframes_raw(session_index); 
        pupil = pupil(pupil_frames_before_this_session+1: pupil_frames_before_this_session+n_pupil_frames_this_session,:); %%% cut out pupil frames not coinciding with this session
        % assume that there is one pupil frame per scope event from this session
        matching_frames = scope_events.onset > scopetiming.scope_start_timepoint & scope_events.onset < scopetiming.scope_stop_timepoint;
        pupil_frame_events = scope_events(matching_frames,{'onset','duration'});
    end


    %%%%%% if video frames and pulses don't match up but pulse-frame difference is less than some fraction of number of pulses, 
    %%%         and range(pupil_frame_events.duration) / nanmean(pupil_frame_events.duration) < some fraction (interpulse time is not too variable),
    %%%         then get first and last pulses and evenly distribute frame timing using frames from avi
    res.excess_pupil_frames = height(pupil_frame_events) - n_pupil_frames_this_session; % number of pulses with no matching frame
    if abs(res.excess_pupil_frames) > 0 %%% if some pulses did not result in an acquired frame
        pupilinfo.nframes_equals_npulses = false;
        if range(pupil_frame_events.duration) / nanmean(pupil_frame_events.duration) < pupil_pulse_duration_range_fraction_limit &&... % check pulse duration variability
                abs(res.excess_pupil_frames/height(pupil_frame_events)) < pupil_framepulse_dif_fraction_limit % check proportion unaccounted-for pulses or frames
            fprintf(['\n Number of pupil video frames exceeds number of pulses by '...
                num2str(res.excess_pupil_frames) '(' num2str(res.excess_pupil_frames/pupilinfo.nframes) ' of frames)\n'])
            firstpulse = pupil_frame_events.onset(1);
            lastpulse = pupil_frame_events.onset(end);
            new_pulse_times = round(linspace(firstpulse, lastpulse, n_pupil_frames_this_session))'; % one pulse per frame, fitting within time between first and last pulse
            pupil_frame_events = table(new_pulse_times, 'VariableNames', {'onset'}); % evenly distribute adjusted onset times
        else 
            error('Pupil frame vs. pulse discrepancy...')
        end
    else
        pupilinfo.nframes_equals_npulses = true;
    end
    pupil = [pupil_frame_events, pupil]; % associate frame timing with pupil data

end


%% reconstruct stimulus aperture mask if it was not saved
    % recreate the mask only if there was only one stim size
    if strcmp(stimpars.experiment_type,'sf_tf_orient_diam') && stimpars.n_diams == 1 && ~isfield(stimpars,'stim_area') && strcmp(stimpars.warping,'warped')
        eye2screen_center = sqrt( stimpars.eye2screen_top_bottom^2 - (stimpars.screen_height/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
        eye2screen_center_pix = eye2screen_center * (stimpars.screenstats.height / stimpars.screen_height); % x0 in Marshel et al; convert to pix
        scrnperpx = round(stimpars.screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
        scrnperpy = round(stimpars.screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
        perp2stimcenterx = scrnperpx - stimpars.stim_center_yx(2); % x pixels between screen center and stim center
        perp2stimcentery = scrnperpy - stimpars.stim_center_yx(1); % y pixels between screen center and stim center
        [xcentermesh ycentermesh] = meshgrid(1:stimpars.screenstats.width,1:stimpars.screenstats.height);
        perp2meshx = scrnperpx-xcentermesh;
        perp2meshy = scrnperpy-ycentermesh;
        stimpars.stim_area = deg2rad(stimpars.diam_minmax(1)/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
            sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
        delete([getfname(stimdata_file), '.mat'])
        save(stimdata_file,'stimpars','stimpar_sets') %%% save aperture into stimdata file
    end

%% stim timing table
stim_events.Properties.VariableNames{find(strcmp('onset',stim_events.Properties.VariableNames))} = 'stim_onset';
stim_events.Properties.VariableNames{find(strcmp('duration',stim_events.Properties.VariableNames))} = 'stim_duration';
stim_events.chan_val = [];
stimpar_sets = [stimpar_sets stim_events]; clear stim_events
 

%% plane number
if exist([getfname(dff_data.fName) '.mat'],'file')
    planeinfo = load([getfname(dff_data.fName) '.mat'],'ops');
    iplane = planeinfo.ops.iplane; % plane number for aligning with TTL pulses; if iplane==1, first frame might be wrong
else 
    iplane = input('Enter plane number: ');
end


%% final steps

stimdur_us = stimpars.stimdur * 1e6; % convert from s to microseconds
stimdur_scans = round(stimdur_us / trig_sampling_interval_us);
isi_us = stimpars.isi * 1e6; % convert from s to microseconds
isi_scans = round(isi_us / trig_sampling_interval_us);
ntrials = height(stimpar_sets);
nimages = size(dff_data.dff,1);
ncells = size(dff_data.dff,2);
res.stimpars = stimpars;
res.iplane = iplane;

