%% add uption to save stim data

% STIM_MAIN_GRID
%   present stimuli in a grid to find cells' RF location
%%% Last updated 2018-10-19 on thermaltake

function stim_main_grid_thermaltake()

%% Setup
stimpars.manual_grid_center = 0;
    stimpars.grid_center = [10 30]; % x and y coords of grid center; if blank, center will be at center of stim screen; overridden by stimpars.manual_grid_center = true
stimpars.keypress_between_trials = 0; % get keyboard input between every trial
stimpars.random_order = 0; % if true, show stim in random order; otherwise, show by rows then columns
stimpars.repetitions = 8; %%  number of full iterations of grid presentation
stimpars.rows = 4;       %% number of rows of stimuli in the grid
stimpars.columns = 6;    %% number of columns of stimuli in the grid
stimpars.h_spacing = 45; %% horizontal distance in degrees between stim centers (not specified in Gao)
stimpars.v_spacing = 45; %% vertical distance in degrees between stim centers (not specified in Gao)
stimpars.diam = 30; %% diameter of stimuli in degrees; divide by 2 to get radius (approximate for off-center stim)... 5 in Gao
stimpars.stimdur = 4; %% duration of presentation of each unique stimulus (2s in Gao)
stimpars.isi = 2; %% time in seconds between end of one stimulus and beginning of next stimulus
stimpars.tf = 2;   %% temporal frequency in hz; peak in Gao ~4hz; 1.5 in Vaic
stimpars.sf = .15; % = grating frequency in cycles/degree if support width and height = rect_wh... ~0.06 peak in Gao; 0.02 in Vaic
stimpars.orient = 0;
stimpars.zeroAngle = 0;            % reference orientation
stimpars.modulateColor = [255 255 255 0];  %RGB Change to grating color
stimpars.backgroundColorOffset = [0 0 0 0];   
stimpars.startphase = 0;   %Degrees
% stimpars.amp = 1e38; % amplitude; grating becomes invisible if this value =>1e(38.6)
stimpars.amp = 1;
stimpars.contrastPreMultiplicator = 1;  
stimpars.background_color = [0 0 0];

% Variables for psychtoolbox and display setup
setupvars.save_stim_file = 1;
setupvars.screen_height = 13;   %% height in inches of stimulus screen; width = 23.5 inches
setupvars.eye2screen_top_bottom = 15; %% distance in inches from eye to both top and bottom of stimulus screen
setupvars.stimScreen = max(Screen('Screens'));
setupvars.do_sync_tests = 0; % sets 'SkipSyncTests' to 0; setting 'SkipSyncTests' to 2 skips all tests
setupvars.ptb_debug_level = 1; % 1 to 4 (4 = most thorough)
setupvars.ptb_verbosity = 0; % 0 = no output, 1 = critical errors, 2 = warnings, 3 = startup info
setupvars.ifi_check = 0; % check that sampled ifi doesn't deviate from nominal ifi by too much
    setupvars.allowableIfiDifference = 1e-5; % measured ifi can differ from nominal refresh interval by this many seconds
    setupvars.ifiSamples = 500; % get this many flip samples for ifi checking
setupvars.missedFlipsWarnProportion = 0.01; % display a warning if the missed flips:total IFIs ratio exceeds this on a trial   

setupvars.send_pulse = 0; 
    setupvars.device_name = 'Dev2'; % device name for addAnalogOutputChannel
    pulsepars.Vstimon = 6; % voltage at which to hold the stim channel a stimulus is being presented
    pulsepars.Vstimoff = 0; % voltage at which to hold the stim channel no stimulus is being presented

%  present RF-mapping stimuli.
grid_runstim(stimpars, pulsepars, setupvars); 
