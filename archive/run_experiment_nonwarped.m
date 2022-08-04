function [] = run_experiment_nonwarpedq(rf_center,override_stimpars)
%Top-level function for calling all other scripts required
%for the experiment.
%%%   override_stimpars allows to change stim params from defaults
%%%%% last edited 12/21/17 on msi

%% switch to gray screen when stimoff?

% % % %%% This version uses non-warped gratings; updated versions use warped
% % % %%% gratings and require the surround gratings to be generated ahead of time 
% % % %%% with generate_ss_stim.m. 

%%% Provide the 'rf_center' argument if it has already been mapped;
%%% otherwise rf_mapping will be run to find the center.

screenstats.screen_height_inches = 13;   %% height in inches of stimulus screen; width = 23.5 inches
screenstats.eye2screen_top_bottom_inches = 9; %% distance in inches from eye to both top and bottom of stimulus screen
save_stim_file = 1;

%% Parameters for stimulus presentation to pass to surround_suppression_stim.m
% Turn on stimpars.SHOW_CENTER to include center grating with parameters
% described by sf_fixed, tf_fixed, and diam_inner. If stimpars.SHOW_CENTER
% is off, only the 'outer' grating will be shown (including central region)
% and diam_inner will have no effect. 
% To only test size tuning, set one of n_sfs or n_tfs to 1 and the other to 0, 
% and set min, max, and fixed values equal to the desired fixed values. 
% Specify below the maximum and minimum and number of values of annulus sf 
% and tf to present. Values of tf and sf will be presented ranging along a 
% log scale between and including the max and min. 
% diam_minmax and diam_inner may be shrunk 
stimpars.SHOW_CENTER = 0;% if true, includes center grating on top of larger grating
stimpars.brightfraction = 0.1; % percentage of maximum brightness to use (divide pix values by this value)
stimpars.sf_fixed = 0.2; %% annulus sf when varying tf and constant inner sf in c/deg (0.02 in in Vaiceliunaite et al. 2013)
stimpars.tf_fixed = 1.5; %% annulus tf when varying sf and constant inner tf in hz (1.5 in Vaiceliunaite et al. 2013)
stimpars.diam_inner = 0; %% diameter of inner grating in degrees
stimpars.orient_offset = 0; %% add this value in degrees to the inner angle to get outer angle; in Webb and Self, 0  suppresses more than 90 (270) do
stimpars.sf_minmax = [0.01 1.6];       %% range of annulus spatial frequencies in cycles per degree (~[0.01 1.6] in Gao)
stimpars.tf_minmax = [0.1 13];       %% range of annulus temporal frequencies in hz (~[0.1 13] in Gao fig 1)
stimpars.diam_minmax = [150 150]; %% range of annulus diameters in degrees ([1 64] in Gao fig 1)
stimpars.angle_range = [0 360]; % range of linearly-space angles to present, not including@ the upper limit
stimpars.n_sfs = 0;      % number of values of annulus spatial frequency to present
stimpars.n_tfs = 1;      % number of values of annulus temporal frequency to present
stimpars.n_diams = 1;   % number of values of annulus diameter to present 
stimpars.nAngles = 1; %% number of times to present the stim set, using one of this number of evenly-spaced orientations 

stimpars.repetitions = 2; %% 10 in Gao
stimpars.stimdur = 4; %% duration of each stimulus presentation in seconds
stimpars.isi = 2; %% inter-stimulus interval in seconds (onset to offset; 0.5 in Vaiceliunaite et al. 2013)
stimpars.do_stimcheck = 0; % if yes, check that stim fit on screen
stimpars.do_sync_tests = 0; 
stimpars.show_splash_screen = 0; % use strict debug level and show white PTB splash screen
stimpars.send_pulse = 0; % send daq pulse

% RF Location-Mapping Stimuli Parameters
% RF mapping stimuli are currently NONwarped
% NOTE: h_spacing and v_spacing may be changed later if the grid does not fit on the screen. 
rf_stim.locmapfile = strcat('locmapping_',date); %% file in which to save rf mapping stimulus information
rf_stim.nAngles = 12; %% CURRENTLY DOES NOT CHANGE Angle; number of evenly-spaced angles to use; this value is also the number of full iterations
rf_stim.rows = 7;       %% number of rows of stimuli in the grid
rf_stim.columns = 7;    %% number of columns of stimuli in the grid
rf_stim.h_spacing = 10; %% horizontal distance in degrees between stim centers (not specified in Gao)
rf_stim.v_spacing = 10; %% vertical distance in degrees between stim centers (not specified in Gao)
rf_stim.diam = 18; %% diameter of stimuli in degrees; divide by 2 to get radius (approximate for off-center stim)... 5 in Gao
rf_stim.dur_stim = 2; %% duration of presentation of each unique stimulus (2s in Gao)
rf_stim.dur_intertrial = 1; %% time in seconds between end of one stimulus and beginning of next stimulus
rf_stim.tfreq = 2;   %% temporal frequency in hz; peak in Gao ~4hz; 1.5 in Vaic
rf_stim.sfreq = .15; % = grating frequency in cycles/degree if support width and height = rect_wh... ~0.06 peak in Gao; 0.02 in Vaic
rf_stim.zeroAngle = 0;            % reference orientation
rf_stim.modulateColor = [255 255 255 0];  %RGB Change to grating color
rf_stim.backgroundColorOffset = [0 0 0 0];   
rf_stim.startphase = 0;   %Degrees
rf_stim.amp = 1e38; % amplitude; grating becomes invisible if this value =>1e(38.6)
rf_stim.contrastPreMultiplicator = 1;  
rf_stim.background_color = [0 0 0];

pulsepars.Vstimon = 6; % voltage at which to hold the stim channel a stimulus is being presented
pulsepars.Vstimoff = 0; % voltage at which to hold the stim channel no stimulus is being presented

%% Run scripts
if save_stim_file
    savefile = strcat('stimdata_',strrep(datestr(now),':',';')); 
else
    checksave = input('Not saving stim file. Press Enter to continue.','s');
    savefile = '';
end
if exist('override_stimpars','var')
    stimpars = copyAllFields(stimpars,override_stimpars);
end
if ~exist('rf_center','var') %% get rf center by mapping if it's not already specified
    [rf_center] = rf_mapping(rf_stim);
end
surround_suppression_stim_nonwarpedq(stimpars, rf_center, pulsepars, savefile, screenstats); 

%%% other possible stim to present: fine-grained determination of tuning
%%% for orientation, sf, tf, direction, speed

% Screen('CloseAll')
commandwindow;
end