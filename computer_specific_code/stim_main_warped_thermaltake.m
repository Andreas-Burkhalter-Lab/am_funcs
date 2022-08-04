function [] = stim_main_warped_thermaltake()
%RUN_EXPERIMENT: Top-level function for calling all other scripts required
%for the surround-suppression experiment.
%%% Provide the 'rf_center' argument if it has already been mapped;
%%% otherwise rf_mapping will be run find the center.
% This script first uses rf mapping and tuningstim_main to determine receptive
% field and stimulus preferences, then tests surround suppression
% modulationg uses these preferences. Tuning preferences are successively
% saved into stimpref. Stimulus data from each tuning experiment are saved into
% individual log files. 
% Use generate_ss_stim.m before ssmod to generate warped surround stimuli. 

%%%%% last updated 19/2/28 on thermaltake


% Parameters for stimulus presentation to pass to surround_suppression_stim.m   
% Specify below the maximum and minimum and number of values of annulus sf and tf to 
% present and the maximum and minimum values. Values of tf and sf will be 
% presented ranging along a log scale between and including the maximum and minimum. 
%%% grats_in_workspace starts to cause missed flips at 2s dur, 60x 60deg gratings
% % trackmouse causes missed flips... instead run track_mouse_movements 
% % other Matlab window... would need to edit trackmouse_callback and
% % specify xResetDist and yResetDist to reanable this option

% prerender_file = 'C:\Users\Burkhalter Lab\Documents\matlab_functions\computer_specific\stimset_prerender_test_sizetuning1.mat'; % name of file where prerendered warped gratings are stored
prerender_file = 'C:\Users\Burkhalter Lab\Documents\matlab_functions\computer_specific\stimset_prerender.mat';
save_stim_file = 1;

%%% parameters for stimuli; gratings will be
%%% warped and prerendered and will not use defaultpars values
stimpars.SHOW_CENTER = 0;% if true, includes center grating on top of larger grating
stimpars.diam_minmax = [30 30]; %% range of annulus diameters in degrees ([1 64] in Gao fig 1)
stimpars.stim_center_yx = [500 540]; %%% % x and y coords of stim center;...... warning: this value might not function if prerend.prerenderedLocAndDiam == true
stimpars.n_diams = 1;   % number of values of annulus diameter to 
stimpars.repetitions = 2; %% 10 in Gao
stimpars.isi = 1; %% inter-stimulus interval in seconds (onset to offset; 0.5 in Vaiceliunaite et al. 2013)
stimpars.include_blank_trials = 1; % include trials with no stimuli (1 per repetition)
stimpars.allow_offscreen_stim = 1; % allow grating sizes which exceed the dimensions of the presentation screen
stimpars.grats_in_workspace = 0; % attempt to load all gratings into the workspace prior to stimulus presentation; may cause missed flips
    stimpars.maxFrameSeqLength = 200; % max # of grats to load at once; large vals interrupt mouseTimer, small vals lengthen isi
stimpars.makeTexDuringIsi = 1; % do Screen('MakeTexture') between trials rather than between flips; prevent missed flips, lengthen isi


% Variables for psychtoolbox and display setup
setupvars.send_pulse = 0; 
    setupvars.device_name = 'Dev2'; % device name for addAnalogOutputChannel
setupvars.overwrite_check = 1;  % check before overwriting stimulus log files (always turn on for experiments)
setupvars.do_sync_tests = 0; % sets 'SkipSyncTests' to 0; setting 'SkipSyncTests' to 2 skips all tests
setupvars.ptb_debug_level = 1; % 1 to 4 (4 = most thorough)
setupvars.ptb_verbosity = 0; % 0 = no output, 1 = critical errors, 2 = warnings, 3 = startup info
setupvars.ifi_check = 0; % check that sampled ifi doesn't deviate from nominal ifi by too much
    setupvars.allowableIfiDifference = 1e-5; % measured ifi can differ from nominal refresh interval by this many seconds
    setupvars.ifiSamples = 500; % get this many flip samples for ifi checking
setupvars.missedFlipsWarnProportion = 0.01; % display a warning if the missed flips:total IFIs ratio exceeds this on a trial    
    
%% Stim channel pulse parameters
pulsepars.Vstimon = 6; % voltage at which to hold the stim channel as stimulus is being presented (except for rf mapping)
pulsepars.Vstimoff = 0; % voltage at which to hold the stim channel when no stimulus is being presented

%% Run scripts
warped_runstim; 


