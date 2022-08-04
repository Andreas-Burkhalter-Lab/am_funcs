function [] = stim_main_nonwarped_thermaltake(override_stimpars)
%Top-level function for calling all other scripts required
%for the experiment.
%%%   override_stimpars allows to change stim params from defaults
%%%%% last edited 18/10/27 on msi

% % % %%% This version uses non-warped gratings; updated versions use warped
% % % %%% gratings and require the surround gratings to be generated ahead of time 
% % % %%% with generate_ss_stim.m. 

screenstats.screen_height_inches = 13;   %% height in inches of stimulus screen; width = 23.5 inches
screenstats.eye2screen_top_bottom_inches = 11; %% distance in inches from eye to both top and bottom of stimulus screen
stimpars.save_stim_file = 1;
    savefile = strcat('stimdata_',strrep(datestr(now),':',';')); 

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
stimpars.orient_fixed = 0; % orientation to use when not testing orientation tuning
stimpars.diam_inner = 0; %% diameter of inner grating in degrees
stimpars.orient_offset = 0; %% add this value in degrees to the inner angle to get outer angle; in Webb and Self, 0  suppresses more than 90 (270) do
stimpars.sf_minmax = [0.01 1.6];       %% range of annulus spatial frequencies in cycles per degree (~[0.01 1.6] in Gao)
stimpars.tf_minmax = [0.1 13];       %% range of annulus temporal frequencies in hz (~[0.1 13] in Gao fig 1)
stimpars.diam_minmax = [150 150]; %% range of annulus diameters in degrees ([1 64] in Gao fig 1)
stimpars.orient_range = [0 360]; % range of linearly-space angles to present, not including@ the upper limit
stimpars.n_sfs = 2;      % number of values of annulus spatial frequency to present; if zero, don't test this param
stimpars.n_tfs = 2;      % number of values of annulus temporal frequency to present; if zero, don't test this param
stimpars.n_orients = 2; %% number of values of annulus orientations to present; if zero, don't test this param
stimpars.n_diams = 1;   % number of values of annulus diameter to present
stimpars.include_blank_trials = 1; % include trials with no stimuli (1 per repetition)

stimpars.stim_center = [500 500]; %  [500 500] = approx screen center
stimpars.repetitions = 2; %% 10 in Gao
stimpars.stimdur = .2; %% duration of each stimulus presentation in seconds
stimpars.isi = .2; %% inter-stimulus interval in seconds (onset to offset; 0.5 in Vaiceliunaite et al. 2013)
stimpars.do_stimcheck = 0; % if yes, check that stim fit on screen
stimpars.do_sync_tests = 0; 
stimpars.show_splash_screen = 0; % use strict debug level and show white PTB splash screen

stimpars.send_pulse = 0; % send daq pulse
    pulsepars.Vstimon = 6; % voltage at which to hold the stim channel a stimulus is being presented
    pulsepars.Vstimoff = 0; % voltage at which to hold the stim channel no stimulus is being presented

%% Run scripts
grating_stim(stimpars, pulsepars, savefile, screenstats); 

