%%%%% This is an old version of the receptive field mapping analysis
%%%%% function, which includes orientation mapping in addition to receptive
%%%%% field center location mapping. A number of sections involving
%%%%% actually determining orientation preference are incomplete; if
%%%%% complete, they should use the von mises fitting function circ_vmpar
%%%%% and, if a prediction from the fit is to be made, the function
%%%%% circ_vmpdf and maybe circ_vmrnd. The sections for rf center location mapping
%%%%% are not updated. 

function [varargout] = analyze_rf_responses(loc_resp_file, orient_resp_file,...
    rec_chans, stim_chan, savedir, save_append, overwrite_ok, avg_method)
%%% Last updated 4/29/2015
%ANALYZE_RF_RESPONSES Determine stimulus preferences of recorded cells. 
% This function is to be run on the recording computer in conjunction with
% running rf_mapping on the stimulus computer. This function constitutes Step 4
% in the outline below, after response data has been acquired via Merec2. 
% To skip location-preference or orientation-preference mapping, set
% the 'loc_resp_file'/'orient_resp_file' arguments to non-strings, and press
% Enter when prompted for filenames. 
% 'rec_chans' is a vector of the electrode channels to analyze for RF preferences.
%       These values are fed into autosnip2.
% 'stim_chan' is the channel through which stimulus timing information 
%       was sent by rf_mapping (assumed to be the same for loc and orient).
% 'save_append' is a string (the subject/session ID) to be appended to the
%     names of files to be created; if empty, will append the date
% avg_method is one of the following and determines how multiple channels are analyzed:
%   'fitall' (default): fit gaussians to all channels, then average
%   'meanall': average responses from all channels, then fit Gaussian
%   'meannormall': average channel-normalized responses from all channels, then fit Gaussian   
%   'onechan': only map rf preferences for the first channel in 'rec_chans'

%% Outline (from rf_mapping)
% [square brackets] indicate whether action is performed on recording computer 
%   (via ssh/nomachine from stimulus computer) or stimulus computer (via this program). 
% 
% 1. [stim] start this program, map approximate RF by ear, save location.
% 2. [rec] start merec2 recording with duration greater than duration of RF
%   mapping stimuli, saving file to 'rf_location_response_[subject/experiment ID]'.
% 3. [stim] send stimuli in grid along with distinct trigger-channel
%   signals for each stimulus parameter value.
% 4. [rec] using matlab program on [rec], analyze responses to each
%   parameter value, then save response information into a .mat file 
%   'rf_map_spikecounts_[subject/experiment ID]' and to a generic ascii file.
% 5. [stim] import spike count information via 'ssh.... head [generic file]'. 
% 6. repeat steps 2-5 for other next stimulus parameters, replacing 'location' with 
%   appropriate parameter name (orientation)

%% Parameters 
pulse.Vbase = 0;    %%% must match the value of this variable in rf_mapping
max_vsteps = 10;    %% maximum number of values of positive target voltages on the stim chan; usually 10
generic_dir = '/home/illya/Andrew/recordings/for_stimulus_computer'; %directory for files for stim comp to find; must match expected dir in rf_mapping

if ~exist(stim_chan,'var')
    stim_chan = input('Stimulus channel?');
end
if numel(stim_chan)>1
    error('Only one channel may be assigned as the stimulus channel.')
end

if ~exist(avg_method,'var') 
    avg_method = 'fitall'; 
elseif ~any(strcmp(avg_method,{'fitall' 'meanall' 'meannormall' 'onechan'})) % we don't recognize avg_method
    disp('Specified channel averaging method not recognized; defaulting to ''fitall''.')
    avg_method = 'fitall';
end

if ~exist(rec_chans,'var')
    rec_chans = input('Recording channel(s) to analyze for RF preferences?');
end
if avg_method == 'onechan'
    if numel(rec_chans) > 1
    warning(sprintf(['More than one recording channel was specified while using averaging method ''onechan''.',
        'Only channel %g will be analyzed.'],rec_chans(1)))
    end
    rec_chans = rec_chans(1);
end

if ~exist(save_append,'var')
    save_append = input('ID to append to file? (Leave blank to append today''s date.)');
    if isempty(save_append)
        save_append = date; % may not be identical to generic file date, but will be identical across orient and loc
    end
end



% RF location-mapping parameters

%%% if using a photodiode, could perform checking of photodiode times...
%%% send flash at actual stimulus onset, then determine whether the
%%% corresponding stim channel pulse arrived at a significantly different
%%% time; provide warning if disparity is above a threshold (maybe perform
%%% check of this timing disparity before running anything), maybe if below
%%% a threshold, set the stimchan event time to the nearest photodiode time
%%% .... would need a simple function for appending to the first frame...
%%% unfulfilled if statement take 3usecs, so could add 'if frame==1, show
%%% the photodiode circle' to the ptb frames loop without significantly
%%% approaching ifi

loc_opt_autosnip2.resume = 0;
loc_opt_autosnip2.feedback_channels = 63;
loc_opt_autosnip2.do_filter = 0;
loc_opt_autosnip2.compare_hz60 = 0;
loc_opt_autosnip2.snip_range = [-25 75]; % adjusted from standard [10 -30] because assuming that scanrate = 25khz.... maybe should be shorter
loc_opt_autosnip2.time4data_segment = 5;
loc_opt_autosnip2.snip_channels = rec_chans;
loc_opt_autosnip2.polarity = 0;

locmap_stim.iterations = 1; %%% must match the value of this variable in rf_mapping
locmap_stim.rows = 4;%%% must match the value of this variable in rf_mapping
locmap_stim.columns = 4; %%% must match the value of this variable in rf_mapping
loc_resp_window = 2; %% time in seconds in which to look for spikes after a location-mapping stimulus
% % %      below not needed if not predicting z
% % % loc_gauss_interp = 5; % number of spaces to interpolate between contiguous grid spaces when predicting z from fmgaussfit

if ischar(loc_resp_file) % if name of RF location response file to analyze was specified
    loc_opt_autosnip2.files = loc_resp_file;
end

% Orientation preference-mapping parameters
orient_opt_autosnip2.resume = 0;
orient_opt_autosnip2.feedback_channels = 63;
orient_opt_autosnip2.do_filter = 0;
orient_opt_autosnip2.compare_hz60 = 0;
orient_opt_autosnip2.snip_range = [-25 75]; % adjusted from standard [10 -30] because assuming that scanrate = 25khz.... maybe should be shorter
orient_opt_autosnip2.time4data_segment = 5;
orient_opt_autosnip2.snip_channels = rec_chans;
orient_opt_autosnip2.polarity = 0;

orientmap_stim.iterations = 1; %%% must match the value of this variable in rf_mapping
orientmap_stim.nAngles = 16; %%% must match the value of this variable in rf_mapping
orientmap_stim.angle_range = [0  360*(1-1/orientmap_stim.nAngles)]; %must match the value of this variable in rf_mapping

orient_resp_window = 2; %% time in seconds in which to look for spikes after a orient-mapping stimulus
% % %      below not needed if not predicting z
% % % orient_gauss_interp = 30; %% number of orientations to interpolate between presented orienations when predicting orientation response

% Promt to check that specified parameters match what was presented from the stimulus computer.
input(sprintf(['If the following parameters match the stimulus script, press Enter:\n',...
    'Location-mapping stimulus rows = %g \n Location-mapping columns = %g\n',...
    'Orientation-mapping # angles = %g \n Stimulus channel = %g'],...
    locmap_stim.rows, locmap_stim.columns, orientmap_stim.nAngles, stim_chan));

if exist(orient_resp_file,'var') % if name of orientation response file to analyze was specified
    orient_opt_autosnip2.files = orient_resp_file;
end

%% Assign save directory and perform analysis. 
if ~exist('savedir','var')  % if save directory not specified, save everything in a new folder in the working directory
    savedir = pwd;
end

if ~isdir(savedir) %% if savedir does not exist yet
    mkdir(savedir); %% make a directory in which to save analysis results 
elseif ~(exist(overwrite_ok, 'var') && overwrite_ok)   %% unless specified by 'overwrite_ok', warn if overwriting an existing directory
    go_on = input('Save directory already exists. OK to overwrite?');
    if ~(strcmp(go_on,'y') || strcmp(go_on,'yes'))
        disp('Will not overwrite directory; quitting analyze_rf_responses...')
        return
    end
    disp('Overwriting save directory.')
else
    disp('Overwriting save directory.')
end

compute_rf_center(loc_opt_autosnip2, loc_resp_window, locmap_pars, pulse,... %% use function below for RF location mapping
    max_vsteps, stim_chan, savedir, generic_dir, save_append, loc_gauss_interp); 
disp('RF location analysis complete. Press Enter on the stimulus computer.') % next, present orientation-preference mapping stimuli
input('Press Enter after orientation-preference-response file is acquired.') % once orientation-preference mapping stimuli have been presented
compute_orient_pref(orient_opt_autosnip2, orient_resp_window, orientmap_pars, pulse,... %% use function below for orientation preference mapping
    max_vsteps, stim_chan, generic_dir, savedir, generic_dir, save_append, orient_gauss_interp);
disp('Oriention preference analysis complete. Press Enter on the stimulus computer.') % proceed with the experiment








%% Function for computing RF center location
function [loc_done] = compute_rf_center(opt_autosnip2, resp_window, stim_pars, pulse_pars,...
        max_vsteps, stim_chan, savedir, generic_dir, save_append, avg_method)
    % input locmap_stim for stim_pars, pulse for pulse_pars
    if ~isfield(opt_autosnip2,'files')
        opt_autosnip2.files = input('RF location-mapping response filename? (Press Enter to skip RF center location mapping.)');
    end
    
    if isempty(opt_autosnip2.files) %% if no location response filename was provided when prompted, skip RF location mapping
        disp('No RF location response filename provided... skipping RF center location mapping...')
        loc_done = 1;
        return
    end
    
    disp(sprintf('Computing RF center location (post-stimulus response window = %gs)...',resp_window))
    
    if strcmp(opt_autosnip2.files(end-5:end), '.merec')  
        opt_autosnip2.files = opt_autosnip2.files(1:end-7);  %% remove extension
    end
    
    merec_obj = merecmm(strcat(opt_autosnip2.files,'.merec'));
    stimwave = merec_obj([stim_chan], [1:merec_obj.nscans]);
    window_scans = resp_window * merec_obj.scanrate;    % convert window size from seconds to scans
    event_bins = extract_event_timing(stimwave, pulse_pars.Vbase, window_scans);  %% get stimulus timing information from the stimulus channel

    % Checks: is #events correct, is each row/column combo found once per
    % iteration, and are iterations identical? 
    %%%%%%%%%%%%%%%% maybe add a check that stim duration is the expected length (do in event parsing function)    
    if stim_pars.iterations*(stim_pars.rows * stim_pars.columns) ~= size(event_bins,2) % check that event_bins is the expected size
        error(['Event timing-parsing error: number of detected stim events (%g)',...
            ' does not equal of iterations (%g) x specified rows (%g) x specified columns (%g) = %g.'],...
            size(event-timing,2), stim_pars.iterations, stim_pars.rows, stim_pars.columns, stim_pars.rows*stim_pars.columns);
    else  % if the number of events is correct
        first_iter = event_bins(:, 1 : (size(event_bins,2)-1)/stim_pars.iterations ); % first iteration of stimuli
        for row = stim_pars.rows
            for column = stim_pars.columns
                if numel(find( 10*row + column == first_iter(3,:))) ~= 1 % if this row/column combination isn't represented exactly once
                    error(['Event timing-parsing error: in each iteration, each row/column combination',...
                        'must be represented exactly once on the stim channel.'])
                end
            end
        end
        if event_bins(3,:) ~= repmat(first_iter(3,:), 1, stim_pars.iterations); % if each iteration isn't identical
            error('Event timing-parsing error: each iteration must be identical')
        end
    end
    
    % Get spike snippets from the electrode data. 
    autosnip2(opt_autosnip2);
    autosort(snipfile2sortheader(strcat(opt_autosnip2.files,'.ssnp'), savedir));
    load('overview.mat')
    snipdata = sorthead_from_raw(sorthead); %% maybe show random examples of snips to visually check that snips are good
    
    % Find mean number of spike snippets in the response window
    % following each stimulus type. All times should be measured in scans. 
    mean_spike_counts = count_snips(snipdata.sniptimes, event_bins(2,:),...
        stim_pars.iterations, resp_window); % use the spike counting function below

    % Fit a 2D Gaussian to the data to compute point of peak response (RF center).  
    switch avg_method
    %%%%%%%    under the options 'meanall' and 'meannormall,' we assume that rf centers are close enough for all
    %%%%%%%    selected (ipsishank) sites - could then check this
    %%%%%%%    assumption offline afterwards and maybe only analyse
    %%%%%%%    (for size tuning etc) sites/cells with rf centers
    %%%%%%%   closely matching the site-averaged computed rf
    %%%%%%%    center (also serves as a check on penetration orthogonality to
    %%%%%%%    cortical surface)
        
    %%% could use points from all trials, rather than just the means, for
    %%% fitting the Gaussian
    %%%%%     probably should include some checking of how good the
    %%%%%     Gaussian fit is and how different the peak is from the
    %%%%%     surrounding values
% %     %%%%% maybe exclude channels from the averaging below if they
% %     %%%%% don't show much spiking response and/or don't show clear
% %     %%%%% location preference based on statistical test(anova, or
% %     %%%%% test that we gain information by fitting the gaussian
% %     %%%%% over just taking the average) or peak amplitude of the 2d Guassian 
        case 'onechan'
            [RF_CENTER amps gaussfit] = loc_gaussfit(mean_spike_counts,...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps); 
        case 'fitall'
            [xypeaks amps gaussfit] = loc_gaussfit(mean_spike_counts,...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps); 
            %%% check for significant differences in rf center between chans
            RF_CENTER = bsxfun(@times,xypeaks,amps)/mean(amps); %weight by amps: minimize effect of unresponsive chans  
        case 'meanall'
            counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
            [RF_CENTER amps gaussfit] = loc_gaussfit(mean(counts_mat),...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps);
        case 'meannormall'
            counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
            counts_mat =  bsxfun(@rdivide,counts_mat,mean(counts_mat,2)); % normalize to mean response on each channel
            [RF_CENTER amps gaussfit] = loc_gaussfit(mean(counts_mat),...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps);
    end

    % Save the necessary variables. 
    %%% in generic file, maybe add: order of stimuli, approximate stim/stimchan timing, Vbase
    save(strcat(savedir,'/rf_location_mapping_',save_append)); %% save all rf location-mapping variables into a permanent directory
    savethis = [];
    savethis(2,:) = clock; % for checking that this file was created in the same session as session when stimuli were presented
    savethis(1,:) = RF_CENTER;
    delete(strcat(generic_dir,'/loc_response_generic')); % remove old generic files
    save(strcat(generic_dir,'/loc_response_generic'),savethis,'-ascii'); %save summary file for stim comp. to find
    
    % Delete everything except the results and generic files. 
    delete(strcat(opt_autosnip2.files,'.env'))
    delete(strcat(opt_autosnip2.files,'.vlv'))
    delete(strcat(opt_autosnip2.files,'.ssnp'))
    delete('overview.mat')
    delete(strcat('chan',str2double(opt_autosnip2.snip_channels),'/autosort_info.mat'))
    if length(extractfield(dir(strcat('chan',str2double(opt_autosnip2.snip_channels))),'name')) < 3; % if chan[recording_chan(1)] is now empty
        rmdir(strcat('chan',str2double(opt_autosnip2.snip_channels)));
    end
    
    loc_done = 1; % output indicating that we don't need to do RF location mapping again
end








%% Function for computing orientation preference
function [orient_done] = compute_orient_pref(opt_autosnip2, resp_window, stim_pars, pulse_pars,...
        max_vsteps, stim_chan, savedir, generic_dir, save_append, avg_method)
    % input orientmap_stim for stim_pars, pulse for pulse_pars
    if ~isfield(opt_autosnip2,'files')
        opt_autosnip2.files = input('Orientation response-mapping filename? (Press Enter to skip orientation preference mapping.)');
    end
    
    if isempty(opt_autosnip2.files) %if no orientation response filename was provided when prompted, skip orient preference mapping
        disp('No orientation response filename provided... skipping orientation preference mapping...')
        orient_done = 1;
        return
    end
    
    disp(sprintf('Computing orientation preference (post-stimulus response window = %gs)...',resp_window))
    
    if strcmp(opt_autosnip2.files(end-5:end), '.merec')  
        opt_autosnip2.files = opt_autosnip2.files(1:end-7);  %% remove extension
    end
    
    merec_obj = merecmm(strcat(filename,'.merec'));
    stimwave = merec_obj([stim_chan], [1:merec_obj.nscans]);
    window_scans = resp_window * merec_obj.scanrate;    % convert window size from seconds to scans
    event_bins = extract_event_timing(stimwave, pulse_pars.Vbase, max_vsteps, window_scans);  %% get stimulus timing information from the stimulus channel 
    
    % Checks: is #events correct, is each orientation index found once per
    % iteration, and are iterations identical? 
    if length(event_bins(3,:))~=stim_pars.nAngles % check that event_bins is the expected size
        error(['Event timing-parsing error: number of stim events represented on stim channel (%g)',...
                ' do not match number of angles specified in this script (%g).'],...
                length(event_bins(3,:)), stim_pars.nAngles)
    else  % if the number of events is correct
        first_iter = event_bins(:, 1 : (size(event_bins,2)-1)/stim_pars.iterations ); % first iteration of stimuli
        if sort(event_bins(3,:)) ~= 11:11+stim_pars.nAngles %if each angle isn't found exactly once in iteration 1
            error('Event timing-parsing error: each stimulus must be represented exactly once in each iteration.')
        elseif event_bins(3,:) ~= repmat(first_iter(3,:), 1, stim_pars.iterations); % if each iteration isn't identical
            error('Event timing-parsing error: each iteration must be identical.')
        end
    end
        
    % Get spike snippets from the electrode data. 
    autosnip2(opt_autosnip2);
    autosort(snipfile2sortheader(strcat(opt_autosnip2.files,'.ssnp'), savedir));
    load('overview.mat')
    snipdata = sorthead_from_raw(sorthead); 

    % Find mean number of spike snippets in the response window
    % following each stimulus type. All times should be measured in scans.  
    %%%%%%%%%%%%%%%%%%%%%%% the 'switch' section is copy-pasted from the
    %%%%%%%%%%%%%%%%%%%%%%% similar section in rf location mapping and not
    %%%%%%%%%%%%%%%%%%%%%%% properly customized for orientation mapping yet
    switch avg_method
        case 'onechan'
            [RF_CENTER amps gaussfit] = orient_vmfit(mean_spike_counts,...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps); 
        case 'fitall'
            [xypeaks amps gaussfit] = orient_vmfit(mean_spike_counts,...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps); 
            RF_CENTER = bsxfun(@times,xypeaks,amps)/mean(amps); %weight by amps to minimize effects of unresponsive chans  
        case 'meanall'
            counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
            [RF_CENTER amps gaussfit] = orient_vmfit(mean(counts_mat),...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps);
        case 'meannormall'
            counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
            counts_mat =  bsxfun(@rdivide,counts_mat,mean(counts_mat,2)); % normalize to mean response on each channel
            [RF_CENTER amps gaussfit] = orient_vmfit(mean(counts_mat),...
                first_iter(3,:),stim_pars,opt_autosnip2.snip_channels,max_vsteps);
    end
    
    %%%%% this section not customized to work with the multiple options
    %%%%% for avg_method yet
    [junk orient_ind] = ismember(first_iter(3,:),sort(first_iter(3,:))); %ordinal value of event labels in the list of all labels
    mean_spike_counts = count_snips(snipdata.sniptimes, event_bins(2,:), stim_pars.iterations,...
        resp_window); % use the spike counting function below
    mean_spike_counts = [first_iter(3,:);orient_ind;mean_spike_counts];...%labels (row1); orientations (row2); spikes (row3)
    sorted_counts = mean_spike_counts(3, sort(mean_spike_counts(2,:))); % get spike counts in order of orient index
    sorted_angles = deg2rad(linspace(stim_pars.angle_range(1), stim_pars.angle_range(2), stim_pars.nAngles));

    angles_as_dist = []; % initialize 
    for ind = 1:length(sorted_counts) %%turn spike counts into distribution for von Mises model fitting
        angles_as_dist = [angles as dist; sorted_angles(ind)*ones(sorted_counts(ind),1)];
    end
    
    %%%%%     probably should include some checking of how good the
    %%%%%     Gaussian fit is and how different the peak is from the
    %%%%%     surrounding values... AB says ANOVA should be done (on the
    %%%%%     different orientations presented) to determine that there is
    %%%%%     orientation selectivity at all; really need to do sum of 2 von
    %%%%%     Mises dists at 180deg separation for single-unit analysis
    %%%%%     (what Gao did)
    [preferred_angle kappa] = circ_vmpar(angles_as_dist); % fit von Mises function to orientation responses
    preferred_angle = rad2deg(preferred_angle); % convert preferred angle to degrees
    %%% optionally, predict a finer distribution with circ_vmpdf for plotting
    
    % Save the necessary variables. 
    save(strcat(savedir,'/rf_orientation_mapping_',save_append)); %% save all rf location-mapping variables into a permanent directory
    savethis = [];
    savethis(2,:) = clock; % for checking that this file was created in the same session as session when stimuli were presented
    savethis(1,:) = preferred_angle;
    delete(strcat(generic_dir,'/orient_response_generic')); % remove old generic files
    save(strcat(generic_dir,'/orient_response_generic'),savethis,'-ascii'); %save summary file for stim comp. to find
  
    % Delete everything except the results and generic files.
    delete(strcat(opt_autosnip2.files,'.env'));
    delete(strcat(opt_autosnip2.files,'.vlv')); 
    delete(strcat(opt_autosnip2.files,'.ssnp')); 
    delete('overview.mat');
    delete(strcat('chan',str2num(opt_autosnip2.snip_channels),'/autosort_info.mat'));
    if length(extractfield(dir(strcat('chan',str2num(opt_autosnip2.snip_channels))),'name')) < 3; % if chan[recording_chan(1)] is now empty
        rmdir(strcat('chan',str2num(opt_autosnip2.snip_channels)));
    end
    

    
    orient_done = 1; % output indicating that we don't need to do orientation mapping again
end










%% Function for extracting and checking stim channel signals
% This function assumes the following: V is held at Vbase when no stimuli
% are present, then shifts to a negative integer at stimulus onset, then
% shifts to a positive integer and is held there until stim offset, at which
% time V returns to Vbase; also, that no stim are present at t=0 and Vbase
% is an integer. 
% Input the merec stimchan wave, the events-off voltage, and response window in
% scans. The output 'stim_info' contains event onsets in row 1, event
% durations in row 2, and compressed event labels in row 3. 
function stim_info = extract_event_timing(wave, Vbase, max_vsteps, window) 
    dig_wave = round(wave);  %% convert the stimchan wave to digital (nearest integer)
    transitions(1,:) = [1, find(~diff(dig_wave)==0) + 1]; % event timings, starting with initial state (1)... +1 is to find post-transition time
    transitions(2,:) = dig_wave(transitions(1,:)); %% get stim channel values after transitioning from one value to another
    
    % Delete all event transitions immediately following the preceding transition, 
    % in case it takes more than one timestep to complete a transition. 
    trans_diff = [diff(transitions(1,:)) 0]; %% get time length in time steps between events, append 0 for last timestep
    count = 0;
    for i = 1:size(transitions,2)
        if trans_diff(i)==1 %% if the following transition is 1 timestep away from this transition
            count = count+1;
        elseif count  %% if this transition is 1 away from the preceding but not the following transition
            transitions(1,i-count+1:i) = -1;        %% mark problematic transitions with -1 for deletion
            count = 0;
        end
    end
    transitions = transitions(:,transitions(1,:)>0);  %% eliminate all transitions one timestep after the preceding transition
    
    % Check that stim events follow the prescribed order of of Vbase,negative,positive,Vbase....positive,Vbase
    if ~rem(length(transitions)-1,3)
        error('Event timing-parsing error: number of events must equal 3 * nStim + 1')
    elseif any(transitions(2,1:3:end) ~= Vbase ||...    % if Vbase values are inccorrect
           any(transitions(2,2:3:end) >= 0 || ...       % ... or negative values are incorrect
           any(transitions(2,3:3:end) <= 0              % .... or positive values are incorrect
       error('Event timing-parsing error: events must follow this sequence: [Vbase, -, +, Vbase, -.... +, Vbase].'
    end
    
    % Compress labels: reverse sign of negative events and multiply by max_vsteps, then add the following positive event, 
    % then delete the positive events (positive events only encode information about the preceding negative event). 
    transitions(2,2:3:end) = -max_vsteps*transitions(2,2:3:end) + transitions(2,3:3:end); % negative event = tens place, positive event = ones place
    transitions(:,3:3:end) = []; % eliminate positive-event columns
    
    % Convert to event-block form as function output: first row is stim onset, 
    % second row is duration, third row is event label. 
    stim_info = [transitions(1,:); [diff(transitions(1,:)) 0]; transitions(3,:)];%add 2nd row: event durations (ignore last nonstim event duration)
    stim_info(:,1:2:end) = [];    %% eliminate stim-absent (Vbase) events
    
    % Check that response window is shorter than the time between stim onsets. (This will not the check final stim vs. end of recording.) 
    if any(diff(stim_info(1,:) < window))
        error('Event timing-parsing error: response window to analyze must be shorter than time between stimulus onsets.')
    end
end

%% Function for counting stimulus-evoked spike snippets after each stimulus event
% tsnips is the cell array of column vectors of spike times
% tevents is the vector of stimulus onset times
% iters is the number of iterations of the full cycle of stimuli (stim_pars.iterations)  
% window is the amount of time after each stimulus onset in which to look for spikes responding to this stimulus 
%%% counts_mean is a vector of length = # unique event labels listing the 
%%%     number of spike following  the corresponding event within time 
%%%     'window', averaged over all iterations
function counts_mean = count_snips(tsnips_cellarray, tevents, iters, window)
    events_per_iter = length(tevents)/iters;
    counts_mean = cell(size(tsnips_cellarray)); % row array
    for chan = 1:length(tsnips_cellarray)
        clear tsnips
        tsnips = tsnips_cellarray{chan};
        countvector = NaN(size(tevents)); % initialize
        for event = 1:length(tevents)  %% for each event, get number of snips following that event 
            countvector(event) = length(find(tsnips > tevents(event) &  tsnips <= tevents(event)+window));
        end
        counts_mean{1,chan} = mean(reshape(countvector,iters,events_per_iter));  %% mean of each event label over all iterations
    end
end

%% Function for fitting 2D Gaussian a grid of spiking responses
% 'counts_in' is a cell array of row vectors of length locmap_stim.rows x locmap_stim.columns, 
%     where each cell contains the grid-spiking responses from one channel.
% 'event_labels' is a row vector of the event labels in the order they appear
%     in the cells of 'grid_in'
% 'pars' is the stim_pars struct argument from compute_rf_center
% 'max_vsteps' is maximum number of values of positive target voltages on the stim chan
%%% 'fit_out' contains the outputs from fmgaussfit, where the elements in
%%%         each field are the results for the channel named by the  
%%%         corresponding element in the 'channel' field
%%% 'rf_center_xy' a vector of the x,y coordinates of the fitted peak [fitresult(5:6)]; rows are chans    
%%% 'amp' is a vector of the peak values of the fitted 2d Gaussian [fitresult(1)]; rows are chans 
% See added notes in fmgaussfit line 73 for description of fitresult values. 
function [rf_center_xy amp fitout] = loc_gaussfit(counts_in, event_labels, pars, channels, max_vsteps)
    fitout = struct('channels', channels, 'meanspikes_grid_z', {}, 'fitresult', {},...
        'zfit', {}, 'fiterr', {}, 'zerr', {}, 'resnorm', [], 'rr', []); 
    [meanspikes_grid_x meanspikes_grid_y] = meshgrid(1:pars.columns, 1:pars.rows); 
    for chan = 1:length(counts_in)
        mean_spike_counts{chan} = [event_labels; max_vsteps*floor(event_labels/max_vsteps);...%event labels (row1); grid-rows (row2)
           rem(event_labels-1,max_vsteps)+1; counts_in{chan}]; %% grid-columns (row3); spike counts (row4)
        fit_out.meanspikes_grid_z{chan} = NaN(stim_pars.rows, stim_pars.columns); % initialize
            for i = 1:size(mean_spike_counts{chan},2)         %% put the event-response spike data into grid form
                meanspikes_grid_z(mean_spike_counts{chan}(2,i), mean_spike_counts{chan}(3,i))...
                    = mean_spike_counts{chan}(4,i);
            end
        [fitout.fitresult{chan}, fitout.zfit{chan}, fitout.fiterr{chan}, fitout.zerr{chan},...%fitresult=row vector
            fitout.resnorm(chan), fitout.rr{chan}] = ...  % fit the grid-spike data with a 2D Gaussian
            fmgaussfit(meanspikes_grid_x, meanspikes_grid_y, fitout.meanspikes_grid_z{chan});
        rf_center_xy(chan,:) = fitout.fitresult{chan}(5:6);
        amp(chan,1) = fitresult{chan}(1);
    end
     %%% to get a z_predict for plotting, comment in below code and a provide a 'gaussfit_interp_spaces' parameter 
% % %     [x_predict y_predict] =...      %% make the finer x-y plane for predicting z values
% % %         meshgrid(1:1/(1+gaussfit_interp_spaces):stim_pars.columns, 1:1/(1+gaussfit_interp_spaces):stim_pars.rows); 
% % %     z_predict = par(7) + fitresult(1)*exp(...    %% generate a 3D surface from the 2D Gaussian fitted parameters
% % %          -(((x_predict-fitresult(5)).*cosd(fitresult(2))+...
% % %           (y_predict-fitresult(6)).*sind(fitresult(2)))./fitresult(3)).^2-...
% % %           ((-(x_predict-fitresult(5)).*sind(fitresult(2))+...
% % %           (y_predict-fitresult(6)).*cosd(fitresult(2)))./fitresult(4)).^2);
end

%% Function for fitting a von Mises distribution to a circular distribtion of spiking responses
%%%%%%%%%%%%%%% this sub-function is incomplete
function fit_out = orient_vmfit(counts_in, event_labels, pars)
    
    for chan = 1:length(counts_in)
       
        fit_out{chan} = 
    end
end






