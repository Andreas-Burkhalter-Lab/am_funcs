function [counts] = rf_main(resp_file_rf, rec_chans,...
     savedir, save_prepend, overwrite_ok, avg_method, stim_file)
%RF_MAIN Determine stimulus preferences of recorded cells. 
% This function is to be run on the recording computer in conjunction with
% running rf_mapping on the stimulus computer. This function constitutes Step 4
% in the outline below, after response data has been acquired via Merec2. 
% To skip location-preference or orientation-preference mapping, set
% the 'resp_file_rf'/'orient_resp_file' arguments to non-strings, and press
% Enter when prompted for filenames. 
%%% INPUTS: 
% 1. 'resp_file_rf' is the .merec file containing electrode responses.
% 2. 'rec_chans' is a vector of the electrode channels to analyze for RF preferences.
%       (Use matlab/merecmm, not dumbbox/merec2 labels)
%       Input 'all' or [] to analyze all recorded channels except stim channels. 
% 3. 'savedir' is the directory in which to put created files.
% 4. 'save_prepend' is a string to be prepended to the
%       names of files to be created; if empty, will prepend the date
% 5. 'overwrite_ok' toggles overwriting existing files without user check.
% 6. 'avg_method' is one of the following and determines how multiple channels are analyzed:
%   'skip' (default): do not do 2d Gaussian fitting
%   'fitall': fit gaussians to all channels, then average
%   'meanall': average responses from all channels, then fit Gaussian
%   'meannormall': average channel-normalized responses from all channels, then fit Gaussian   
%   'onechan': only map rf preferences for the first channel in 'rec_chans'
% 7. 'stim_file' (optional) is the stim log file created by rf_loc_map_stim.m
%%% You must be in the directory containing resp_file_rf to prevent errors with readheader.m. 
%%%%% Last updated 1/30/2015 on vivid

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
%%%%%% add: if it's a 16ch recording, make ch55 a reserved stimchan.... or
%%%%%% change merec layout file to not include ch55 for 16ch recordings
stim_chan = 63; % the (merecmm) channel through which stimulus timing information was sent by rf_mapping   
reserved_stimchans = [55 63]; % these channels will not be set as recording channels by default
rerun_rf_getevents = 0; % run rf_getevents even if the .trialdata file already exists
rerun_autosnip2_autosort = 0; % rerun autosnip2 and autosort even if the files exist
keep_autosnip2_autosort_files = 1; % keep outputs of autosnip2 and autosort except for the .vlv and .env files
resp_window = 2;        %% time in seconds in which to look for spikes after a location-mapping stimulus
do_filter = 1;          % if on, use filtering option in autosnip2; significantly increases autosnip2 time
    highpass = 70;     % hz; lowpass will be set to 0.48*merec_obj.scanrate (no effect if ~do_filter)
gaussfit_interp_spaces = 0; % # of spaces to interpolate between grid points when plotting the z prediction from 2D Gaussian fitting
generic_dir = '/home/Andrew/recordings/for_stimulus_computer'; %directory for files for stim comp to find; must match expected dir in rf_mapping

manual_trial_check = 0; % manually check high-spiking trials and remove noisy trials    
    opt_trialCheck.zScoreThresh = .5; % above-average-spiking trials with this z-score or greater will be flagged
    opt_trialCheck.nChansThresh = 1; % trials with at least this many flagged channels will be presented for checking
    opt_trialCheck.maxRows = 4; %  max number of subplot rows to use for plotting traces - fill all rows before columns
    opt_trialCheck.maxColumns = 1; %  max number of subplot rows to use for plotting traces
    opt_trialCheck.chans = [15];      % channels to check; if empty, defaults to all chans in opt_autosnip2.snip_channels
    opt_trialCheck.keepRaw = 1;     % keep data on trials flagged as noisy for reference
    opt_trialCheck.plotTrialWaveforms = 1; % plot the full waveform of the trial alongside the snips    
    
stimpars_rf.nAngles = 8; %%% must match the value of this variable in rf_mapping
stimpars_rf.rows = 9;%%% must match the value of this variable in rf_mapping
stimpars_rf.columns = 9; %%% must match the value of this variable in rf_mapping
pulse.Vbase = 0;    %%% must match the value of this variable in rf_mapping
max_vsteps = 10;    %% maximum number of values of positive target voltages on the stim chan; usually 10

if ~exist('resp_file_rf','var') || isempty(resp_file_rf)
    resp_file_rf = uigetfile('*.merec','Select a RF location-mapping response file.');
end
opt_autosnip2.files = resp_file_rf;
if length(opt_autosnip2.files)>5 && strcmp(opt_autosnip2.files(end-5:end), '.merec')  
    opt_autosnip2.files = opt_autosnip2.files(1:end-6);  %% remove extension
end
merec_obj = merecmm(strcat(opt_autosnip2.files,'.merec'));

opt_autosnip2.thresh_factor = 1; % =2 in Vaiceliunaite et al. 2013: thresh  = twice 6*median
opt_autosnip2.resume = 0; % 0 = default
opt_autosnip2.feedback_channels = 63;
opt_autosnip2.do_filter = do_filter;
opt_autosnip2.compare_hz60 = 0; % 0 = default
opt_autosnip2.snip_range = round([-1 2] * merec_obj.scanrate*1e-3); % ms in brackets; default = [-1 3]
opt_autosnip2.time4data_segment = 5; % 5 = default
opt_autosnip2.polarity = 0;
opt_autosnip2.do_in_time_order = 1; % 1 =  = default
opt_autosnip2.extension = '.ssnp'; % '.ssnp' = default
opt_autosnip2.iterative_threshold = 0; 
opt_autosnip2.no_env = 1; 
opt_autosnip2.no_redo_vlv = 1;

if do_filter
    passband = [highpass 0.48*merec_obj.scanrate]; % highpass slightly below Nyquist
    opt_snippetfile = snipoptions(struct(... % make the conditioning filter
        'Fs',merec_obj.scanrate,'bandpass',passband)); 
else
    passband = [];
    opt_snippetfile = []; % use snippetfile defaults only
    highpass = 'not filtered';
end

% % %      below line not needed if not predicting z
% % % loc_gauss_interp = 5; % number of spaces to interpolate between
% contiguous grid spaces when predicting z from fmgaussfit

if exist('stim_file','var') && ~isempty(stimfile)
% % %     matobj = matfile(stimfile);
    load(stim_file,'stim_centers','manual_center');
    stimpars_rf = pars;
end

%% Variable Checks
if numel(stim_chan)>1
    error('Only one channel may be assigned as the stimulus channel.')
end

if ~exist('avg_method','var')  || isempty(avg_method)
    avg_method = 'skip'; 
elseif ~any(strcmp(avg_method,{'skip' 'fitall' 'meanall' 'meannormall' 'onechan'})) % we don't recognize avg_method
    disp('Specified channel averaging method not recognized; defaulting to ''fitall''.')
    avg_method = 'fitall';
end

if ~exist('save_prepend','var') || isempty(save_prepend)
    save_prepend = input('ID to prepend to file? (Leave blank to append response file name.) ','s');
    if isempty(save_prepend)
        save_prepend = ['analyzed',getfname(resp_file_rf)]; % may not be identical to generic file date
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
opt_autosnip2.snip_channels = rec_chans; % 'semiauto' = default (snip entire file, then pick channels)
nchans = length(rec_chans);

if strcmp(avg_method, 'onechan')
    if numel(rec_chans) > 1
    warning(['More than one recording channel was specified while using averaging method ''onechan''.',...
        'Only channel %g will be analyzed.'],rec_chans(1))
    end
    rec_chans = rec_chans(1);
end

% Prompt to check that specified parameters match what was presented from the stimulus computer.
input(sprintf(['\nIf the following parameters match the stimulus script, press Enter:\n',...
    '   Location-mapping stimulus rows = %g \n   Location-mapping columns = %g\n',...
    '   Angles (iterations) = %g\n   Stimulus channel = %g\n'],...
    stimpars_rf.rows, stimpars_rf.columns, stimpars_rf.nAngles, stim_chan));

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

if ~rerun_rf_getevents && exist([getfname(resp_file_rf) '.trialdata'],'file') % check to see if we already get event timing data
    fprintf('Found event timing data in ''%s''; will not rerun ss_getevents.\n',...
        [getfname(resp_file_rf) '.trialdata'])
    load([getfname(resp_file_rf) '.trialdata'],'trialdata_rf','event_bins','window_scans','nevents','first_iter','-mat'); % load event data
else
    rf_getevents; %% extract event timing
end

rf_getspikes; % get spike snippets from the electrode data

%% Fit a 2D Gaussian to the data to compute point of peak response (RF center).  
%%      will need to edit 'mean spike counts' to take the 'counts mean' column from 'counts'
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
        mean_spike_counts = {counts{chanNames{1},'avg'}};
        [RF_CENTER amps gaussfit zpred] = rf_gaussfit(mean_spike_counts,...
            first_iter(3,:),stimpars_rf,opt_autosnip2.snip_channels,...
            max_vsteps,gaussfit_interp_spaces); 
    case 'fitall'
        mean_spike_counts = cell(1,size(counts,1)-1);
        for i = 1:size(counts,1)-1 %% simplify this: make the subfunction take the dataset as argument rather than cell array
            mean_spike_counts{i} = counts{i,'avg'};
        end     
        [xypeaks amps gaussfit zpred] = rf_gaussfit(mean_spike_counts,...
            first_iter(3,:),stimpars_rf,opt_autosnip2.snip_channels,...
            max_vsteps,gaussfit_interp_spaces); 
        %%% still need to check for significant differences in rf center between chans
        RF_CENTER = bsxfun(@times,xypeaks,amps)/mean(amps); %weight by amps to minimize effect of unresponsive chans  
        RF_CENTER = mean(RF_CENTER); % take weighted mean of rf center from all analyzed channels
    case 'meanall' 
        counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
        [RF_CENTER amps gaussfit zpred] = rf_gaussfit(mean(counts_mat),...
            first_iter(3,:),stimpars_rf,opt_autosnip2.snip_channels,...
            max_vsteps, gaussfit_interp_spaces);
    case 'meannormall'
        counts_mat = cell2mat(mean_spike_counts'); % rows=channels, columns = event labels
        counts_mat =  bsxfun(@rdivide,counts_mat,mean(counts_mat,2)); % normalize to mean response on each channel
        [RF_CENTER amps gaussfit zpred] = rf_gaussfit(mean(counts_mat),...
            first_iter(3,:),stimpars_rf,opt_autosnip2.snip_channels,...
            max_vsteps,gaussfit_interp_spaces);
end

%% Save analysis results and clean up.
if ~manual_trial_check
    clear opt_trialcheck; % clear options structure so we know trial checking wasn't run
end
saveIfExists([savedir,filesep,save_prepend,'_rf_mapping_'],...%% save rf location-mapping variables into a permanent directory
    {'RF_CENTER','xypeaks','counts','rec_chans','resp_window','snipdata','tsnips','zpred',...
    'do_filter','passband','flagged_trials','trialnoise_rf','nFlaggedTrials','manual_trial_check',...
    'amps','avg_method','resp_file_rf','stim_chan','stimpars_rf','mean_spike_counts',...
    'pulse', 'trialdata_rf','gaussfit','generic_dir',...
    'window_scans','response_time_total','baseline_time_total','max_vsteps',...
    'merec_obj','opt_autosnip2','opt_snippetfile','sorthead'});
if ~strcmp(avg_method,'skip')
%%% in generic file, maybe add: order of stimuli, approximate stim/stimchan timing, Vbase
    figure; surf(mean(zpred,3)); % display 2D Gaussian fit of receptive field 
    zpredict = dataset(repmat({cell(size(zpred(:,:,1)))},size(zpred,3),1),...
        'ObsNames',counts.Properties.ObsNames(1:size(zpred,3)),'VarNames',{'zpred'});
    for indChan = 1:size(zpred,3) %%%% put zpred into dataset form, change name to 'zpredict'
        zpredict{indChan,'zpred'} = zpred(:,:,indChan);
    end
    clear zpred
    title('2D Gaussian fit of rec_chans responses')
    savethis = [];
    savethis(2,:) = clock; % for checking that this file was created in the same session as session when stimuli were presented
    savethis(1,1:2) = RF_CENTER;
    delete(strcat(generic_dir,'/loc_response_generic')); % remove old generic files
% % %     save(strcat(generic_dir,'/loc_response_generic'),'savethis','-ascii'); %save summary file for stim comp. to find
end

% display mean spikes 
figure; surf(counts{'avg','avg'}); 
title('Average grid response from rec_chans')
rotate3d on;

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


%% rf_gaussfit: subfunction for fitting 2D Gaussian to a grid of spiking responses
% 'counts_in' is a cell array with each cell containing the grid of mean 
%       spikes for one channel
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
function [rf_center_xy, amp, fitout, zpredict] = rf_gaussfit(counts_in, event_labels,...
    pars, channels, max_vsteps, gaussfit_interp_spaces)

    if any(any(isnan(cell2mat(counts_in))))
        error('Some grid locations average == NaN (all trials were eliminated during manual trial checking).')
    end
    fitout = struct('channels', channels); 
    %%%%% do we need to switch x and y in the following 3 lines? %%%%%
    [meanspikes_grid_x meanspikes_grid_y] = meshgrid(1:pars.columns, 1:pars.rows); 
    [x_predict y_predict] =...      %% make the finer x-y plane for predicting z values
        meshgrid(1:1/(1+gaussfit_interp_spaces):pars.columns, 1:1/(1+gaussfit_interp_spaces):pars.rows); 
    zpredict = -inf(size(x_predict,1),size(x_predict,2),length(channels));
    for chan = 1:length(counts_in)         
            %%%%% do we need to switch x and y in the following 2 lines? %%%%%%%%%%%
        [fitout.fitresult{chan}, fitout.zfit{chan}, fitout.fiterr{chan}, fitout.zerr{chan},...%fitresult=row vector
            fitout.resnorm(chan), fitout.rr{chan}] = ...  % fit the grid-spike data with a 2D Gaussian
            fmgaussfit(meanspikes_grid_x, meanspikes_grid_y, counts_in{chan});
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



