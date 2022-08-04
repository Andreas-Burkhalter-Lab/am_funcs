function [] = run_experiment(rf_center)
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
% For ssmod, angles, sf, tf, duration, and outer grating amplitude must all be specified 
%   by generate_ss_stim.m.
%%%%% last edited 1/27/16 on stim comp

% Specify which experiments to run.
do_rfMapping = 1; % do rf mapping.... rf mapping will be skipped if stim center is provided as input to run_experiment, even if this var==1
do_stimPreferences = 0; % get tuning preferences of a single unit for orientation, sf, tf, and size
do_ssModulation = 1; % run ss modulation stimuli

% Parameters for stimulus presentation to pass to surround_suppression_stim.m   
% Specify below the maximum and minimum and number of values of annulus sf and tf to 
% present and the maximum and minimum values. Values of tf and sf will be 
% presented ranging along a log scale between and including the maximum and minimum. 
%%% grats_in_workspace starts to cause missed flips at 2s dur, 60x 60deg gratings
% % trackmouse causes missed flips... instead run track_mouse_movements 
% % other Matlab window... would need to edit trackmouse_callback and
% % specify xResetDist and yResetDist to reanable this option
prerender_file = 'stimset_prerender_test_sizetuning.mat'; % name of file where prerendered warped gratings are stored

%%% Default stim parameters for stimulus preference tuning; as tuning
%%% properties for preferred rf center, orientation, sf, tf, and stim size
%%% are measured, these measured values will be used in place of the values below. 
defaultpars.diam = 18; %% diameter of stimuli in degrees; divide by 2 to get radius (approximate for off-center stim)... 5 in Gao
defaultpars.orient = 0; %% stim orientation 
defaultpars.isi = 0.5; %% time in seconds between end of one stimulus and beginning of next stimulus
defaultpars.tf = 2;   %% temporal frequency in hz; peak in Gao ~4hz; 1.5 in Vaic
defaultpars.sf = .1; % = grating frequency in cycles/degree if support width and height = rect_wh... ~0.06 peak in Gao; 0.02 in Vaic
defaultpars.dur_stim = 2; %% duration of presentation of each unique stimulus (2s in Gao)
defaultpars.modulateColor = [255 255 255 0];  %RGBA Change to grating color
defaultpars.backgroundColorOffset = [0 0 0 0];   
defaultpars.startphase = 0;   % starting phase for drifting gratings in degrees
defaultpars.amp = 1e38; % amplitude; grating becomes invisible if this value =>1e(38.6)
defaultpars.contrastPreMultiplicator = 1;  
defaultpars.background_color = [0 0 0];

%%% parameters for surround suppression modulation stimuli; these will be
%%% warped and prerendered and will not use defaultpars values
ssmod_stim.SHOW_CENTER = 1;% if true, includes center grating on top of larger grating
ssmod_stim.diam_minmax = [13 64]; %% range of annulus diameters in degrees ([1 64] in Gao fig 1)
ssmod_stim.n_diams = 1;   % number of values of annulus diameter to 
ssmod_stim.repetitions = 10; %% 10 in Gao
ssmod_stim.isi = 1; %% inter-stimulus interval in seconds (onset to offset; 0.5 in Vaiceliunaite et al. 2013)
ssmod_stim.grats_in_workspace = 1; % attempt to load all gratings into the workspace prior to stimulus presentation; may cause missed flips
    ssmod_stim.maxFrameSeqLength = 200; % max # of grats to load at once; large vals interrupt mouseTimer, small vals lengthen isi
ssmod_stim.makeTexDuringIsi = 1; % do Screen('MakeTexture') between trials rather than between flips; prevent missed flips, lengthen isi
ssmod_stim.trackmouse = 0; 
ssmod_stim.tuningParameter = 'ssmod';
ssmod_stim.filetag = 'ssmod'; %%% tag for log file

% RF Location-Mapping Stimuli Parameters
% RF mapping stimuli are currently NONwarped
% NOTE: h_spacing and v_spacing may be changed later if the grid does not fit on the screen. 
rf_stim.getRfCenterRemotely = 0; % if on, get rf center from diesel via plink; otherwise, have user enter results from diesel
rf_stim.repetitions = 8; %% iterations of the full grid; originally intended to change angle but does not
rf_stim.rows = 9;       %% number of rows of stimuli in the grid
rf_stim.columns = 9;    %% number of columns of stimuli in the grid
rf_stim.h_spacing = 10; %% horizontal distance in degrees between stim centers (not specified in Gao)
rf_stim.v_spacing = 10; %% vertical distance in degrees between stim centers (not specified in Gao)
rf_stim.tuningParameter = 'rf'; 
rf_stim.filetag = 'rf'; %%% tag for log file

% Orientation tuning testing parameters; currently nonwarped stimuli.
% Evenly spaced angles from 0 to 360-angle_increment will be presented.
orient_stim.n_orients = 2; %% number of evenly-spaced angles to use;
orient_stim.zeroorient =  0;  % reference orientation for creating orientation distribution
orient_stim.repetitions = 2; %% number of times to show each angle
orient_stim.tuningParameter = 'orient'; 

% Spatial frequency tuning testing parameters; nonwarped
sf_stim.n_sfs = 1;      % number of values of annulus spatial frequency to present
sf_stim.sf_minmax = [0.01 1];      %% range of  spatial frequencies in cycles per degree (~[0.01 1.6] in Gao)
sf_stim.repetitions = 1; %% number of times to show each sf
sf_stim.tuningParameter = 'sf';

% Temporal frequency tuning testing parameters; nonwarped
tf_stim.n_tfs = 1;      % number of values of annulus temporal frequency to present
tf_stim.tf_minmax = [0.1 13];       %% range of temporal frequencies in hz (~[0.1 13] in Gao fig 1)
tf_stim.repetitions = 1; %% number of times to show each tf
tf_stim.tuningParameter = 'tf';

% Stimulus size tuning parameters; nonwarped
size_stim.n_diams = 2; 
size_stim.diam_minmax = [3 64]; %% range of annulus diameters in degrees ([1 64] in Gao fig 1)
size_stim.repetitions = 2; %% number of times to show each size
size_stim.tuningParameter = 'diam';

% Variables for psychtoolbox and display setup
setupvars.send_pulse = 0; 
setupvars.overwrite_check = 1;  % check before overwriting stimulus log files (always turn on for experiments)
setupvars.mouseScanRate = 300; % hz
setupvars.mouseBufferFactor = 1.5; %%% multiply expected number of scans by this factor when pre-allocating mouse_movements values
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
pulsepars.Vstimoff_dur_post = 0.5; %% length of time in seconds during which stim chan is held at Vstimoff after stimulus ends.... should = isi
pulsepars.V1_dur = 0.5; %% duration in seconds of first voltage value in each pulse if there are two stim-on pulse values (rf mapping only)

%% Run scripts
commonvars = setup_surround_suppression(setupvars, prerender_file);
if ~exist('rf_center','var') %% get rf center by mapping if it's not already specified
    stimpref.rf_center = rfstim_main(rf_stim,pulsepars,defaultpars,commonvars);
else
    stimpref.rf_center = rf_center; clear rf_center; % can't input to run_experiment as a structure field
end
if do_stimPreferences % present stimuli for determining tuning for orientation, sf, tf, size, in that order
    [orient_stimrec, orient_stim, stimpref] = tuningstim_main(orient_stim,stimpref,pulsepars,defaultpars,commonvars);
    [sf_stimrec, sf_stim, stimpref] = tuningstim_main(sf_stim,stimpref,pulsepars,defaultpars,commonvars);
    [tf_stimrec, tf_stim, stimpref] = tuningstim_main(tf_stim,stimpref,pulsepars,defaultpars,commonvars);
    [size_stimrec, size_stim, stimpref] = tuningstim_main(size_stim,stimpref,pulsepars,defaultpars,commonvars);
end  
if do_ssModulation
    ssmodstim_main; 
end

% Cleanup after stim presentation
Screen('CloseAll')
Screen('Preference','SkipSyncTests',commonvars.origSyncLevel); % reset to original sync testing setting
Screen('Preference', 'VisualDebugLevel', commonvars.origDebugLevel); % reset to iriginal PTB debug level
Screen('Preference','Verbosity',commonvars.origVerbosity); % reset to original PTB warning mode
commandwindow;
end
