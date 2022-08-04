function [counts] = analyze_rf_responses_updated(stim_chan, loc_resp_file, ...
    rec_chans, savedir, save_append, overwrite_ok, avg_method, stim_file)
%%% Last updated 7/10/2015 on diesel

% maybe present rf mapping stim in sphere-projected space; this analysis
% should still work (then update grid-relative to screen-relative
% conversion in rf_mapping.m)

%ANALYZE_RF_RESPONSES Determine stimulus preferences of recorded cells. 
% This function is to be run on the recording computer in conjunction with
% running rf_mapping on the stimulus computer. This function constitutes Step 4
% in the outline below, after response data has been acquired via Merec2. 
% To skip location-preference or orientation-preference mapping, set
% the 'loc_resp_file'/'orient_resp_file' arguments to non-strings, and press
% Enter when prompted for filenames. 
% 'rec_chans' is a vector of the electrode channels to analyze for RF preferences.
%       (Use matlab/merecmm, not dumbbox/merec2 labels)
%       Use shankmaps.m to get channels from a specific shank. 
%       Input 'all' to analyze all recorded channels except stim channels. 
% 'stim_chan' is the channel through which stimulus timing information 
%       was sent by rf_mapping (matlab/merecmm, not dumbbox/merec2 labels).
% 'save_append' is a string (the subject/session ID) to be appended to the
%       names of files to be created; if empty, will append the date
% 'avg_method' is one of the following and determines how multiple channels are analyzed:
%   'skip': do not do 2d Gaussian fitting
%   'fitall' (default): fit gaussians to all channels, then average
%   'meanall': average responses from all channels, then fit Gaussian
%   'meannormall': average channel-normalized responses from all channels, then fit Gaussian   
%   'onechan': only map rf preferences for the first channel in 'rec_chans'
% 'stim_file' (optional) is the stim log file created by rf_loc_map_stim.m
%%% You must be in the directory containing loc_resp_file to prevent errors with readheader.m. 


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

%% Parameters 
pulse.Vbase = 0;    %%% must match the value of this variable in rf_mapping
max_vsteps = 10;    %% maximum number of values of positive target voltages on the stim chan; usually 10
generic_dir = '/home/illya/Andrew/recordings/for_stimulus_computer'; %directory for files for stim comp to find; must match expected dir in rf_mapping
reserved_stimchans = [55 63]; % these channels will not be set as recording channels by default
skip_autosnip2_autosort = 0; % skip autosnip2 and autosort
keep_autosnip2_autosort_files = 1; % keep outputs of autosnip2 and autosort except for the .vlv and .env files
resp_window = 2; %% time in seconds in which to look for spikes after a location-mapping stimulus
gaussfit_interp_spaces = 5; % # of spaces to interpolate between grid points when plotting the z prediction from 2D Gaussian fitting

locmap_stim.nAngles = 5; %%% must match the value of this variable in rf_mapping
locmap_stim.rows = 9;%%% must match the value of this variable in rf_mapping
locmap_stim.columns = 9; %%% must match the value of this variable in rf_mapping

if ~exist('loc_resp_file','var') || isempty(loc_resp_file)
    loc_resp_file = uigetfile('*.merec','Select a RF location-mapping response file.');
end
opt_autosnip2.files = loc_resp_file;
if length(opt_autosnip2.files)>5 && strcmp(opt_autosnip2.files(end-5:end), '.merec')  
    opt_autosnip2.files = opt_autosnip2.files(1:end-6);  %% remove extension
end
merec_obj = merecmm(strcat(opt_autosnip2.files,'.merec'));

opt_autosnip2.resume = 0;
opt_autosnip2.feedback_channels = 63;
opt_autosnip2.do_filter = 0;
opt_autosnip2.compare_hz60 = 0;
opt_autosnip2.snip_range = [-1 2] * merec_obj.scanrate*1e-3; % ms in brackets; standard = [1 3]
opt_autosnip2.time4data_segment = 5;
opt_autosnip2.polarity = 0;
opt_autosnip2.do_in_time_order = 1;
opt_autosnip2.extension = '.ssnp';
opt_autosnip2.iterative_threshold = 0; 

% % %      below line not needed if not predicting z
% % % loc_gauss_interp = 5; % number of spaces to interpolate between
% contiguous grid spaces when predicting z from fmgaussfit

if exist('stim_file','var') && ~isempty(stimfile)
% % %     matobj = matfile(stimfile);
    matobj = load(stimfile); % diesel doesn't have matfile.m function
    locmap_stim = matobj.pars;
    pulse = matobj.pulse;
    load(stim_file,'stim_centers','manual_center');
end

%% Variable Checks
if ~exist('stim_chan','var') || isempty(stim_chan)
    stim_chan = input('Stimulus channel?');
end
if numel(stim_chan)>1
    error('Only one channel may be assigned as the stimulus channel.')
end

if ~exist('avg_method','var')  || isempty(avg_method)
    avg_method = 'fitall'; 
elseif ~any(strcmp(avg_method,{'skip' 'fitall' 'meanall' 'meannormall' 'onechan'})) % we don't recognize avg_method
    disp('Specified channel averaging method not recognized; defaulting to ''fitall''.')
    avg_method = 'fitall';
end

if ~exist('save_append','var') || isempty(save_append)
    save_append = input('ID to append to file? (Leave blank to append today''s date.)','s');
    if isempty(save_append)
        save_append = [date,'_rf_mapping']; % may not be identical to generic file date
    end
end



origdir = pwd; 
cd(fileparts(which([opt_autosnip2.files '.merec']))); % go to the dir of the response file

if ~exist('rec_chans','var') || isempty(rec_chans) || strcmp(rec_chans,'all')
    disp(sprintf(['Using all channels except '...
        '[%s].'], num2str(reserved_stimchans)));
    allchans = merec_obj.channels;
    rec_chans = allchans(arrayfun(@(x)(~any(x==reserved_stimchans)), allchans));
end
opt_autosnip2.snip_channels = rec_chans;

if strcmp(avg_method, 'onechan')
    if numel(rec_chans) > 1
    warning(['More than one recording channel was specified while using averaging method ''onechan''.',...
        'Only channel %g will be analyzed.'],rec_chans(1))
    end
    rec_chans = rec_chans(1);
end

% Prompt to check that specified parameters match what was presented from the stimulus computer.
input(sprintf(['\nIf the following parameters match the stimulus script, press Enter:\n',...
    'Location-mapping stimulus rows = %g \nLocation-mapping columns = %g\n',...
    'Angles (iterations) = %g\nStimulus channel = %g'],...
    locmap_stim.rows, locmap_stim.columns, locmap_stim.nAngles, stim_chan));

% Assign save directory. 
if ~exist('savedir','var') || isempty(savedir);  % if save directory not specified, save everything in a new folder in the working directory
    savedir = pwd;
end
if ~isdir(savedir) %% if savedir does not exist yet
    mkdir(savedir); %% make a directory in which to save analysis results 
elseif ~(exist('overwrite_ok', 'var') && overwrite_ok)   %% unless specified by 'overwrite_ok', warn if overwriting an existing directory
    go_on = input('Save directory already exists. Enter ''y'' to overwrite.','s');
    if ~(strcmp(go_on,'y') || strcmp(go_on,'Y') || strcmp(go_on,'yes'))
        disp('Will not overwrite directory; quitting analyze_rf_responses...')
        return
    end
    disp('Overwriting save directory.')
else
    disp('Overwriting save directory.')
end

%% Start analysis
disp(sprintf('Computing RF center location (post-stimulus response window = %gs)...',resp_window))
stimwave = merec_obj([stim_chan], [1:merec_obj.nscans]);
window_scans = resp_window * merec_obj.scanrate;    % convert window size from seconds to scans
event_bins = extract_event_timing(stimwave, pulse.Vbase, max_vsteps, window_scans);  %get stim timing info from stim channel

% Checks: is #events correct, is each row/column combo found once per
% iteration, and are iterations identical? 
%%%%%%%%%%%%%%%% maybe add a check that stim duration is the expected length (during event parsing section)    
if locmap_stim.nAngles*(locmap_stim.rows * locmap_stim.columns) ~= size(event_bins,2) % check that event_bins is the expected size
    error(['Event timing-parsing error: number of detected stim events (%g)',...
        ' does not equal of iterations (%g) x specified rows (%g) x specified columns (%g) = %g.'],...
        size(event_bins,2), locmap_stim.nAngles, locmap_stim.rows, locmap_stim.columns, locmap_stim.rows*locmap_stim.columns);
else  % if the number of events is correct
    first_iter = event_bins(:, 1 : (size(event_bins,2)/locmap_stim.nAngles)); % first iteration of stimuli
    for row = locmap_stim.rows
        for column = locmap_stim.columns
            if numel(find( 10*row + column == first_iter(3,:))) ~= 1 % if this row/column combination isn't represented exactly once
                error(['Event timing-parsing error: in each iteration, each row/column combination',...
                    'must be represented exactly once on the stim channel.'])
            end
        end
    end
    if event_bins(3,:) ~= repmat(first_iter(3,:), 1, locmap_stim.nAngles); % if each iteration isn't identical
        error('Event timing-parsing error: each iteration must be identical')
    end
end

% Get spike snippets from the electrode data.
opt_autosnip2.files = {[opt_autosnip2.files '.merec']}; % needs to be cell for autosnip2 to handle it correctly
if ~skip_autosnip2_autosort
    autosnip2(opt_autosnip2);
    autosort(snipfile2sortheader(strcat(opt_autosnip2.files{:}(1:end-6),'.ssnp')), savedir);
end
load([savedir '/overview.mat'])
snipdata = sorthead_from_raw(sorthead, opt_autosnip2.snip_channels);
%%% maybe show random examples of snips to visually check that snips are good

%% Count stimulus-evoked spike snippets after each stimulus event
% Find mean number of spike snippets in the response window
% following each stimulus type. All times should be measured in scans. 
%%% counts_mean is a vector of length = # unique event labels listing the 
%%%     number of spike following  the corresponding event within time 
%%%     'window', averaged over all iterations
events_per_iter = length(event_bins(1,:))/locmap_stim.nAngles;
response_time_total = (size(event_bins,2)*resp_window); % total time counted as response
baseline_time_total = (merec_obj.nscans / merec_obj.scanrate) - response_time_total; % total time counted as baseline
emptygrid = {NaN(locmap_stim.rows, locmap_stim.columns)}; 
chanNames = [cellfun(@(x)['ch_',num2str(x)],num2cell(rec_chans),'UniformOutput',false) 'avg'];
iterNames = [cellfun(@(x)['iter_',num2str(x)],num2cell(1:locmap_stim.nAngles),'UniformOutput',false) 'avg'];
counts = dataset({repmat(emptygrid,length(rec_chans)+1,locmap_stim.nAngles+1), iterNames{:}}, ...
    {NaN(length(chanNames),4),'response_hz','baseline_hz','response_total','baseline_total'}, ...
    'ObsNames', chanNames);
for chan = 1:length(snipdata.sniptimes)
    tsnips = snipdata.sniptimes{chan};
    countsmat = NaN(locmap_stim.rows, locmap_stim.columns, locmap_stim.nAngles); % initialize and clear
    for event = 1:length(event_bins(1,:))  %% for each event, get number of snips following that event on this channel
        row = floor(event_bins(3,event)/max_vsteps);
        column = mod(event_bins(3,event),max_vsteps);
        iteration = floor((event-1)/(locmap_stim.rows*locmap_stim.columns)) + 1; 
        counts{chanNames(chan),iterNames(iteration)}(row,column) =...
            length(find(tsnips > event_bins(1,event) & tsnips <= event_bins(1,event)+window_scans));
    end
    for iter = 1:locmap_stim.nAngles
        countsmat(:,:,iter) = counts{chanNames(chan),iterNames(iter)};  
    end
    counts.avg{chanNames(chan)} = mean(countsmat,3); % average for this channel over all iterations
    counts.response_total(chanNames(chan)) = sum(sum(sum(countsmat))); % total snips in response windows for this chan 
    counts.response_hz(chanNames(chan)) = counts.response_total(chanNames(chan)) / response_time_total; % avg hz during resp window
    counts.baseline_total(chanNames(chan)) = size(snipdata.snips{chan},2) - counts.response_total(chanNames(chan)); % total snips during baseline
    counts.baseline_hz(chanNames(chan)) = counts.baseline_total(chanNames(chan)) / baseline_time_total; % avg hz during baseline
end

% Get averages for a given trial across all channels and grand-average.
countsmat = NaN(locmap_stim.rows, locmap_stim.columns, length(rec_chans), locmap_stim.nAngles); % initialize and clear
for iter = 1:locmap_stim.nAngles
    for chan = 1:length(snipdata.sniptimes)
        countsmat(:,:,chan,iter) = counts{chanNames(chan),iterNames(iter)};
    end
    counts{'avg',iterNames(iter)} = mean(countsmat(:,:,:,iter),3); % average over all channels for this iteration
end
counts.avg{'avg'} = mean(mean(countsmat,3),4); % grand average across all channels and iterations
counts.response_total('avg') = nanmean(counts.response_total);
counts.response_hz('avg') = nanmean(counts.response_hz);
counts.baseline_total('avg') = nanmean(counts.baseline_total);
counts.baseline_hz('avg') = nanmean(counts.baseline_hz);

%% Fit a 2D Gaussian to the data to compute point of peak response (RF center).  
%%      will need to edit 'mean spike counts' to take the 'counts mean columnfrom 'counts'
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
%%%% weighting channels by amp may be problematic - may overly favor
%%%% high-rate cells low-rate cells with well-defined rfs; maybe just set
%%%% threshold to eliminate nonresponsive chans then weight all responsive
%%%% chans equally
    case 'skip'
        RF_CENTER = [NaN NaN];
    case 'onechan'
        [RF_CENTER amps gaussfit zpred] = loc_gaussfit(mean_spike_counts,...
            first_iter(3,:),locmap_stim,opt_autosnip2.snip_channels,...
            max_vsteps,gaussfit_interp_spaces); 
    case 'fitall'
        [xypeaks amps gaussfit zpred] = loc_gaussfit(mean_spike_counts,...
            first_iter(3,:),locmap_stim,opt_autosnip2.snip_channels,...
            max_vsteps,gaussfit_interp_spaces); 
        %%% still need to check for significant differences in rf center between chans
        RF_CENTER = bsxfun(@times,xypeaks,amps)/mean(amps); %weight by amps to minimize effect of unresponsive chans  
        RF_CENTER = mean(RF_CENTER); % take weighted mean of rf center from all analyzed channels
    case 'meanall' 
        counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
        [RF_CENTER amps gaussfit zpred] = loc_gaussfit(mean(counts_mat),...
            first_iter(3,:),locmap_stim,opt_autosnip2.snip_channels,...
            max_vsteps, gaussfit_interp_spaces);
    case 'meannormall'
        counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
        counts_mat =  bsxfun(@rdivide,counts_mat,mean(counts_mat,2)); % normalize to mean response on each channel
        [RF_CENTER amps gaussfit zpred] = loc_gaussfit(mean(counts_mat),...
            first_iter(3,:),locmap_stim,opt_autosnip2.snip_channels,...
            max_vsteps,gaussfit_interp_spaces);
end

save([savedir,'/rf_location_mapping_',save_append]); %% save all rf location-mapping variables into a permanent directory
if ~strcmp(avg_method,'skip')
%%% in generic file, maybe add: order of stimuli, approximate
%%% stim/stimchan timing, Vbase
    figure; surf(mean(zpred,3)); % display 2D Gaussian fit of receptive field 
    savethis = [];
    savethis(2,:) = clock; % for checking that this file was created in the same session as session when stimuli were presented
    savethis(1,1:2) = RF_CENTER;
    delete(strcat(generic_dir,'/loc_response_generic')); % remove old generic files
    save(strcat(generic_dir,'/loc_response_generic'),'savethis','-ascii'); %save summary file for stim comp. to find
end

% Delete everything except the results and generic files. 
delete(strcat(opt_autosnip2.files{:}(1:end-6),'.env'))
delete(strcat(opt_autosnip2.files{:}(1:end-6),'.vlv'))
if ~keep_autosnip2_autosort_files
    delete(strcat(opt_autosnip2.files{:}(1:end-6),'.ssnp'))
    delete(strcat(savedir,'/overview.mat'))
    delete(strcat('chan',str2double(opt_autosnip2.snip_channels),'/autosort_info.mat'))
    for ch = 1:length(opt_autosnip2.snip_channels) % delete chanX folders if they only contain autosort_info.mat
        thischan = opt_autosnip2.snip_channels(ch);
        if length(extractfield(dir(strcat(savedir,'/chan',num2str(thischan))),'name')) < 4; 
            rmdir(strcat(savedir,'/chan',num2str(thischan)),'s'); 
        end
    end
end

cd(origdir); % go back to the directory we were in when this function was called
disp('RF location analysis complete. Press Enter on the stimulus computer.')

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    int_wave = round(wave);  %% round stim wave to nearest integer
    transitions(1,:) = [1, find(~diff(int_wave)==0) + 1]; % event timings, starting with initial state (1)... +1 is to find post-transition time
    transitions(2,:) = int_wave(transitions(1,:)); %% get stim channel values after transitioning from one value to another
    
    % Delete all event transitions immediately following the preceding transition, 
    % in case it takes more than one timestep to complete a transition. 
    trans_diff = [diff(transitions(1,:)) 0]; %% get time length in time steps between events, append 0 for last timestep
    count = 0;
    for i = 1:size(transitions,2)
        if trans_diff(i)==1 %% if this transition is 1 timestep before the next transition
            count = count+1; 
        elseif count  %% if this transition is 1 step away from the preceding but not the following transition
            transitions(1,i-count+1:i) = -1;        %% mark problematic transitions with -1 for deletion
            transitions(2,i-count) = transitions(2,i); % use the final stimchan value in this set of adjacent transitions
            count = 0;
        end
    end
    transitions = transitions(:,transitions(1,:)>0);  %% eliminate all transitions one timestep after the preceding transition
    
    % If the initial pre-stim Vbase value is preceded by another stimchan
    % value, give a warning but delete the first value. The first value was
    % probably the stimchan value before stimulus presentation began. 
   if transitions(end,2)==Vbase & transitions(end,1)~=Vbase
       warning(['Ignoring stimchan value preceding first stim-off channel value'...
           ' (value = %g from scan 1 to scan %g).'],transitions(end,1),transitions(1,2));
       transitions = transitions(:,2:end); % delete the preceding nonzero column
   end
    
    % Check that stim events follow the prescribed order of of Vbase,negative,positive,Vbase....positive,Vbase
    if rem(length(transitions)-1,3) % if 'transitions'-1 is NOT divisible by 3
        error('Event timing-parsing error: number of events must equal 3 * nUniqueStim + 1')
    elseif any(transitions(2,1:3:end) ~= Vbase) ||...    % if Vbase values are inccorrect
           any(transitions(2,2:3:end) >= 0) || ...       % ... or negative values are incorrect
           any(transitions(2,3:3:end) <= 0)              % .... or positive values are incorrect
       error('Event timing-parsing error: events must follow this sequence: [Vbase, -, +, Vbase, -.... +, Vbase].')
    end
    
    % Compress labels: reverse sign of negative events and multiply by max_vsteps, then add the following positive event, 
    % then delete the positive events (positive events only encode information about the preceding negative event). 
    transitions(2,2:3:end) = -max_vsteps*transitions(2,2:3:end) + transitions(2,3:3:end); % negative event = tens place, positive event = ones place
    transitions(:,3:3:end) = []; % eliminate positive-event columns
    
    % Convert to event-block form as function output: first row is stim onset, 
    % second row is duration, third row is event label. 
    stim_info = [transitions(1,:); [diff(transitions(1,:)) 0]; transitions(2,:)];%add 2nd row: event durations (ignore last nonstim event duration)
    stim_info(:,1:2:end) = [];    %% eliminate stim-absent (Vbase) events
    
    % Check that response window is shorter than the time between stim onsets. (This will not the check final stim vs. end of recording.) 
    if any(diff(stim_info(1,:) < window))
        error('Event timing-parsing error: response window to analyze must be shorter than time between stimulus onsets.')
    end
end


%% Function for fitting 2D Gaussian to a grid of spiking responses
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
function [rf_center_xy amp fitout zpredict] = loc_gaussfit(counts_in, event_labels,...
    pars, channels, max_vsteps, gaussfit_interp_spaces)
    fitout = struct('channels', channels); 
    %%%%% do we need to switch x and y in the following 3 lines? %%%%%
    [meanspikes_grid_x meanspikes_grid_y] = meshgrid(1:pars.columns, 1:pars.rows); 
    [x_predict y_predict] =...      %% make the finer x-y plane for predicting z values
        meshgrid(1:1/(1+gaussfit_interp_spaces):pars.columns, 1:1/(1+gaussfit_interp_spaces):pars.rows); 
    zpredict = -inf(size(x_predict,1),size(x_predict,2),size(channels));
    for chan = 1:length(counts_in)
        mean_spike_counts{chan} = [event_labels; floor(event_labels/max_vsteps);...%event labels (row1); grid-rows (row2)
           rem(event_labels-1,max_vsteps)+1; counts_in{chan}]; %% grid-columns (row3); spike counts (row4)
        fitout.meanspikes_grid_z{chan} = NaN(pars.rows, pars.columns); % initialize
            for i = 1:size(mean_spike_counts{chan},2)         %% put the event-response spike data into grid form
                fitout.meanspikes_grid_z{chan}(mean_spike_counts{chan}(2,i), mean_spike_counts{chan}(3,i))...
                    = mean_spike_counts{chan}(4,i); % rows and cols should be same as in stim grid
            end
     %%%%% do we need to switch x and y in the following 2 lines? %%%%%%%%%%%
        [fitout.fitresult{chan}, fitout.zfit{chan}, fitout.fiterr{chan}, fitout.zerr{chan},...%fitresult=row vector
            fitout.resnorm(chan), fitout.rr{chan}] = ...  % fit the grid-spike data with a 2D Gaussian
            fmgaussfit(meanspikes_grid_x, meanspikes_grid_y, fitout.meanspikes_grid_z{chan});
        rf_center_xy(chan,:) = fitout.fitresult{chan}(5:6);
        amp(chan,1) = fitout.fitresult{chan}(1);
        % Generate a 3D surface from the 2D Gaussian fitted parameters for
        % plotting.
        zpredict(:,:,chan) = fitout.fitresult{chan}(7) + fitout.fitresult{chan}(1)*exp(...    
             -(((x_predict-fitout.fitresult{chan}(5)).*cosd(fitout.fitresult{chan}(2))+...
              (y_predict-fitout.fitresult{chan}(6)).*sind(fitout.fitresult{chan}(2)))./fitout.fitresult{chan}(3)).^2-...
              ((-(x_predict-fitout.fitresult{chan}(5)).*sind(fitout.fitresult{chan}(2))+...
              (y_predict-fitout.fitresult{chan}(6)).*cosd(fitout.fitresult{chan}(2)))./fitout.fitresult{chan}(4)).^2);
    end

end



