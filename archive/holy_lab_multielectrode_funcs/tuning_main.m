%TUNING_MAIN Create tuning curves for stimulus parameters.
%%% Inputs: 
%   1. resp_file is the .merec file of MEA responses
%   2. stim_file is the stimulus log file created on the stimulus computer
%   3. rec_chans specifies which channels contain electrode data to analyze
%       (default = all channels except 'reserved_stimchans')
%   4. datadir specifies the directory (preferably absolute) in which to save results 
%       and look for autosnip2, autosort, and cass results (default = current directory) 
%   5. save_prepend = label to prepend to saved analysis results file
%%% See run_experiment.m for relevant units.
%%% updated 1/7/16
% Intended for single-unit data analysis, not multiple electrode channels

%%% Last updated 1/7/16
function [] = tuning_main(resp_file, stim_file, rec_chans, datadir,save_prepend)

%% Parameters
stim_chan = 63;
resp_window = 2; % post=stim time (s) in which to count spikes as responding to the stim
rerun_autosnip2_autosort = 0;  
rerun_ss_getevents = 0; % run ss_getevents even if the .trialdata file already exists
manual_clustering = 0; % run cass_sort_apply and treat all manual clusters as separate channels
reserved_stimchans = [55 63]; % do not treat these channels as electrode channels by default...maybe autodetect 16vs64ch?
snip_range_ms = [-1 2]; % spike-relative snippet window in ms (default = [-1 3])
analyzeLimitedTrialRange = 0; % if true, analyze only a range of trials
    trialStartEndToAnalyze = [1 1500]; % has no effect if analyzeLimitedTrialRange = 0
do_filter = 0; % if on, use filtering option in autosnip2; significantly increases autosnip2 time
    highpass = 90; % hz; lowpass will be set to 0.48*merec_obj.scanrate (no effect if ~do_filter)
diamIndsForPrefAngle = 1; %only these diam indices (low = small diam) will be used for computing the pref angle
disp_trialBeingAnalyzed = 0; % output trial being analyzed on command line
n_param_vals_forplotting = 50; % interpolate values to have this many values when plotting curvefit   

manual_trial_check = 1; % manually check high-spiking trials and remove noisy trials    
    opt_trialCheck.zScoreThresh = 1.5; % above-average-spiking trials with this z-score or greater will be flagged
    opt_trialCheck.nChansThresh = 4; % trials with at least this many flagged channels will be presented for checking
    opt_trialCheck.maxRows = 4; %  max number of subplot rows to use for plotting traces - fill all rows before columns
    opt_trialCheck.maxColumns = 8; %  max number of subplot rows to use for plotting traces
    opt_trialCheck.plotTrialWaveforms = 1; % plot the full waveform of the trial alongside the snips
    opt_trialCheck.chans = [];      % channels to check; if empty, defaults to all chans in opt_autosnip2.snip_channels
    opt_trialCheck.keepRaw = 1;     % keep data on trials flagged as noisy for reference
    
opt_autosnip2.thresh_factor = 1; % =2 in Vaiceliunaite et al. 2013: thresh  = twice 6*median
opt_autosnip2.resume = 0;
opt_autosnip2.feedback_channels = 63;
opt_autosnip2.do_filter = do_filter; % 0 = default
opt_autosnip2.compare_hz60 = 0; % 0 = default
opt_autosnip2.time4data_segment = 5; % 5 = default
opt_autosnip2.polarity = -1;
opt_autosnip2.do_in_time_order = 1;
opt_autosnip2.extension = '.ssnp';
opt_autosnip2.iterative_threshold = 0; 
opt_autosnip2.no_env = 1; % 0 = default
opt_autosnip2.no_redo_vlv = 1; % 0 = default

%% Run analysis
analysisTotalTic = tic;
get_response_data; % snip and sort spikes if necessary
tuning_curvefit; % tabulate stimulus-response results and fit tuning curve
analysis_runtime = toc(analysisTotalTic); % does not reflect full runtime if getspikes/autsort/autosnip were already run

%% Save results
[junk file_to_save junk] = fileparts(stim_file);
if exist('datadir','var') && ~isempty(datadir)
    file_to_save = [datadir filesep save_prepend '_' file_to_save '.mat'];
end

%%%%% may not save large variables ('snipdata' if the recording was
%%%%% long)... use -v7.3 to save larger variables
saveIfExists(file_to_save,[{'stimpars','trialdata','stimrec','analyzeLimitedTrialRange','stimPresentationStartTime',...
    'datadir''do_filter','file_to_save','distToPlot','parValsToPlot'...
    'highpass','manual_clustering','manual_trial_check','merec_obj','opt_autosnip2','opt_snippetfile','options',...
    'par_sets','passband','pulsepars','rec_chans','resp_file','resp_window','analysis_runtime','savefile',...
    'snip_file','snip_range_ms','snipdata','sorthead','stim_chan',...
    'stim_file','trialStartEndToAnalyze','window_scans'} tuningVarsToSave]); 
fprintf('\nAnalysis saved into %s\n',file_to_save)

end