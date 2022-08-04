function [ commonvars ] = setup_surround_suppression( setupvars, prerender_file )
%SCREEN_SETUP Basic setup for run_experiment prior to rf_mapping and
%surround_suppresion_stim
%   Open a screen, get ifi, start daq session, make mouse timer.
%%%% last updated on stim comp 9/23/15
commonvars = setupvars; % copy setup preferences into commonvars 

%% fix next line = we might not have a prerender file yet
if exist(prerender_file,'file') && ~isempty(prerender_file)
    load(prerender_file,'screen_height','eye2screen_top_bottom'); % get the same screen specifications used when prerendering
else % use defaults
    screen_height = 13;   %% height in inches of stimulus screen; width = 23.5 inches
    eye2screen_top_bottom = 9; %% distance in inches from eye to both top and bottom of stimulus screen
end
commonvars.screen_height = screen_height;
commonvars.eye2screen_top_bottom = eye2screen_top_bottom;
global daq_sess % causes problems if we clear and recreate the session (hardware stays reserved for first session)

% Start pyschtoolbox and Data Acquisition Toolbox sessions
% Note: all 3 screens (stim, stim-duplicate, and control) have max resolution 
%    of 1920x1080, where aspect ratio of screen resolution = aspect ratio of 
%    physical screen size (square pixels).  
% input('Press Enter when screen is away from eye or turned off.');    %% avoid exposing eye to psychtoolbox splash
commonvars.origSyncLevel = Screen('Preference','SkipSyncTests');
Screen('Preference', 'SkipSyncTests', double(~commonvars.do_sync_tests));
commonvars.origDebugLevel = Screen('Preference', 'VisualDebugLevel', commonvars.ptb_debug_level);
commonvars.origVerbosity = Screen('Preference', 'Verbosity', commonvars.ptb_verbosity);

Screen('CloseAll');
stimScreen = max(Screen('Screens'));
[winPtr, winrect] = Screen('OpenWindow',stimScreen,BlackIndex(stimScreen));  
commonvars.screenstats = Screen('Resolution',stimScreen);

% Get interflip interval.
nominal_refresh = 1/Screen('NominalFrameRate',winPtr);
if commonvars.ifi_check
    [commonvars.ifi ifi_samp std ] = Screen('GetFlipInterval', winPtr, commonvars.ifiSamples);
    if abs(commonvars.ifi - nominal_refresh) > commonvars.allowableIfiDifference % if the measured ifi differs too much from the nominal ifi
        error(['Measured interflip interval more than 1e-5 seconds different from nomimal screen refresh interval:',...
            'IFI - Nominal Refresh = %g'],(commonvars.ifi - nominal_refresh));
    end
else
    commonvars.ifi = nominal_refresh;
    warning('Not computing interflip interval, using ifi = nominal refresh = %g instead.', nominal_refresh);
end

% Set up daq session.
if setupvars.send_pulse
    if isempty(daq_sess) || ~any(size(daq_sess.Channels)>0) % if we somehow cleared the daq session (by stopping run_experiment.m)
        daq_sess = daq.createSession('ni');
        addAnalogOutputChannel(daq_sess,'Dev1','ao0','Voltage'); 
    end
end
    
%% Mouse Tracking
%%% assign the callback function in the stim script
commandScreenStats = Screen('Resolution',1); % stats from the command screen  - where the measurement mouse cursor will be
commonvars.mouseScreenCenter = [commandScreenStats.width/2 commandScreenStats.height/2]; % center of command screen
mouseTimer = timer; 
mouseTimer.ExecutionMode = 'fixedSpacing';
mouseTimer.Period = round(1000/commonvars.mouseScanRate)/1000;
    
commonvars.winPtr = winPtr;
commonvars.winrect = winrect;
commonvars.mouseTimer = mouseTimer;
end

