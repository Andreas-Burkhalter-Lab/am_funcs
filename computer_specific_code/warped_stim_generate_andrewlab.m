function [  ] = warped_stim_generate_andrewlab(  )
%GENERATE_stimpars Create warped grating stimuli and save them
%as .mat files. 
%   Call this function to generate stimulus files before presenting warped
%   grating stimuli. Gratings are drawn to the full screen; to specify location and size,
%   create an aperture in the stimulus presentation script to overlay onto the
%   grating. Full-screen 120-frame presentations take ~700ms to
%   load on the stimulus computer, so use isi of 1s+. 
% last updated 18/8/26 on thermaltake

tic
%% Paremeters
%%% inner grating sf and tf are fixed; inner and outer amp is fixed;
%%% gratings always iso-oriented
%%% maybe add orient_offset option
screen_height = 13;   %% height in inches of stimulus screen; width = 23.5 inches
eye2screen_top_bottom = 14; %% distance in inches from eye to both top and bottom of stimulus screen

prerend.prerenderedLocAndDiam = 0; % prerender only the portion of the grating filling the required location and diameter
    prerend.stim_center = [854  676]; % only save gratings centered here if prerenderedLocAndDiam = 1
    prerend.max_diam = 13; % only save gratings within the area of this diameter if prerenderedLocAndDiam = 1
    prerend.diam_inner = 13; %% diameter of inner grating in degrees... only has an effect if prerenderedLocAndDiam = 1
    prerend.amp_inner = 1; %% amplitude of inner sine gratings (1=max)... only has an effect if prerenderedLocAndDiam = 1
prerend.brightfraction = .1;
prerend.sf_fixed = 0.1; % cpd
prerend.tf_fixed = 2; % hz
prerend.orient_fixed = 0; % deg
prerend.sf_minmax = [0.01 1.6];      %% range of annulus spatial frequencies in cycles per degree (~[0.01 1.6] in Gao)
prerend.tf_minmax = [0.1 13];       %% range of annulus temporal frequencies in hz (~[0.1 13] in Gao fig 1)
prerend.outer_amp = 1;
prerend.orient_range = [0 360];  %%%% degrees; positive is clockwise
prerend.n_sfs = 1;      % number of values of annulus spatial frequency to present
prerend.n_tfs = 1;      % number of values of annulus temporal frequency to present
prerend.n_orients = 1; 
prerend.stimdur = 4; % seconds
prerend.tfsCommonMultiples = 0; % make all tfs factors of the next-smallest tf to reduce # of prerendered gratings
prerend.tfsIfiMultiples = 1; % set all tfs >= than 1/stimdur to integer values of flips/cycle to reduce # of gratings
% prerend.stimScreen = max(Screen('Screens'));  % create stimuli for highest-index screen
prerend.stimScreen = 2;
stimset_filename = 'stimset_prerender'; % save gratings and parameter values into a .mat file of this name


% PTB checks
prerend.prerend_ifi_checked = 0; % check that sampled ifi doesn't deviate from nominal ifi by too much
    allowableIfiDifference = 1e-5; % measured ifi can differ from nominal refresh interval by this many seconds
    ifiSamples = 500; % get this many flip samples for ifi checking
do_sync_tests = 0; % sets 'SkipSyncTests' to 0 during execution of this script; setting 'SkipSyncTests' to 2 skips all tests
ptb_debug_level = 1; % 1 to 4 (4 = most thorough)
ptb_verbosity = 0; % 0 = no output, 1 = critical errors, 2 = warnings, 3 = startup info

warped_stim_generate_makestim();


