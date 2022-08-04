function [computed_rf_center preferred_angle] = rf_mapping_gabors(varargin)
% Old incomplete version of this function that uses Gabors rather than
% gratings as stimuli. 

%%% Last updated 5/1/2015
%RF_MAPPING Map receptive field
%   First get approximate RF location by manually moving stimulus while
%   listening for audio responses. Once an approximate center is found
%   manually, automatically present stimuli in a grid to more precisely
%   find the RF location around the approximate center by analyzing spiking
%   responses. Preferred orientation and other parameters can be
%   subsequently mapped via the same procedures.
%
%   This process currently requires interacting with recording computer in
%   addition to stimulus computer as indicated in Outline. 
%
%   Begin running before bringing screen near eye to avoid exposure to the
%   Psychtoolbox splash screen. 

%% Outline
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
% 6. repeat steps 2-5 for other next stimulus parameters, replacing 'location' with 
%   appropriate parameter name

%% Setup
% Measurements from the physical setup
screen_height = 23.5;   %% height in inches of stimulus screen
eye2screen_top_bottom = 9; %% distance in inches from eye to both top and bottom of stimulus screen

ifi_check = 0; %% compare interflip interval to screen refresh rate; must always be on for real experiments

% get rid of warnings and splash screen... comment out for real experiments
old_visualdebug = Screen('Preference', 'VisualDebugLevel');
Screen('Preference', 'VisualDebugLevel', 0);    

% RF Location-Mapping Stimuli Parameters
%%%%%% maybe give option to show each stimulus more than once (shouldn't
%%%%%% require extra code on the analysis side if pulse code is consistent 
% (Stim channel noise (max-min) = ~0.3V over 5s... we can safely accomodate ~100 stimuli using 2-pulse IDs.)
locmap_stim.iterations = 1; %% number of times to show the entire grid of stimuli
locmap_stim.rows = 4;       %% number of rows of stimuli in the grid
locmap_stim.columns = 4;    %% number of columns of stimuli in the grid
locmap_stim.horz_spacing = 200; %% should input these values as degrees, convert with deg2pixels %% horizontal spacing between center of stimuli in the grid
locmap_stim.vert_spacing = 200; %% should input these values as degrees, convert with deg2pixels%% vertical spacing between center of stimuli in the grid
locmap_stim.v_diam = 90; %% vertical diameter of stimuli in degrees - Gao et al. used 2 degrees
locmap_stim.h_diam = 50; %% horizontal diameter of stimuli in degrees - Gao et al. used 2 degrees
locmap_stim.dur_stim = 2; %% duration of presentation of each unique stimulus
locmap_stim.dur_intertrial = 2; %% time in seconds between end of ones stimulus and beginning of next stimulus
locmap_stim.width = 16000;      %% controls Gabor shape
locmap_stim.height = 14000;      %% controls Gabor shape
locmap_stim.speed = 24.7219;
locmap_stim.Angle = 170;            %Degrees
locmap_stim.modulateColor = [255 255 255 0];  %RGB Change to grating color
locmap_stim.backgroundColorOffset = [0 0 0 0];   
locmap_stim.phase = 0;   %Degrees
locmap_stim.freq = .0059;     %Spatial Frequency - adjusts speed so that temporal frequency does not changelocmap_stim.
locmap_stim.contrast = 1e38; % Gabor becomes invisible if this value =>1e(38.6)
locmap_stim.sc = 50.0;         %Spatial Constant of Gaussian hull function - the sigma value
locmap_stim.nonSymmetric = 0;    % 0 makes ellipsoid, 1 makes a symmetric circle
locmap_stim.aspectratio = 2;         % Ignored if nonSymmetric = 1
locmap_stim.disableNorm = 0;     % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.
locmap_stim.radius = 0.1;         %  These are optional arguments
locmap_stim.contrastPremultiplicator = 1;  
locmap_stim.background_color = [0 0 0];

% Orientation Preference-Mapping Stimuli Parameters
orientmap_stim.iterations = 1;      % number of times to show the entire set of orientations
orientmap_stim.nAngles = 16;    %% number of orientations 
orientmap_stim.angle_range = [0  360*(1-1/orientmap_stim.nAngles)];  %% upper and lower limit on angles to test.... maybe add a manual angle-mapping function like with rf center hand-mapping, then limit this range accordingly
orientmap_stim.v_diam = 2; %% vertical diameter of stimuli in degrees - Gao et al. used 2 degrees
orientmap_stim.h_diam = 2; %% horizontal diameter of stimuli in degrees - Gao et al. used 2 degrees
orientmap_stim.dur_stim = 2; %% duration of presentation of each unique stimulus
orientmap_stim.dur_intertrial = 2; %% time in seconds between end of one stimulus/iteration and beginning of next stimulus/iteration
orientmap_stim.width = 16000;      %% controls Gabor shape
orientmap_stim.height = 14000;      %% controls Gabor shape
orientmap_stim.speed = 24.7219;
orientmap_stim.modulateColor = [255 255 255 0];  %RGB Change to grating color
orientmap_stim.backgroundColorOffset = [0 0 0 0];   
orientmap_stim.phase = 0;   %Degrees
orientmap_stim.freq = .0059;     %Spatial Frequency - adjusts speed so that temporal frequency does not changelocmap_stim.
orientmap_stim.contrast = 1e38; % Gabor becomes invisible if this value =>1e(38.6)
orientmap_stim.sc = 50.0;         %Spatial Constant of Gaussian hull function - the sigma value
orientmap_stim.nonSymmetric = 0;    % 0 makes ellipsoid, 1 makes a symmetric circle
orientmap_stim.aspectratio = 2;         % Ignored if nonSymmetric = 1
orientmap_stim.disableNorm = 0;     % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.
orientmap_stim.radius = 0.1;         %  These are optional arguments
orientmap_stim.contrastPremultiplicator = 1;  
orientmap_stim.background_color = [0 0 0];

% Stim channel parameters
% These values are used for both location mapping and orientation preference mapping.  
pulse.Vbase = 0;  %% stim channel value when no stimulus is present
pulse.Vbase_dur_post = 0.05; %% length of time in seconds during which stim chan is held at Vbase after stimulus ends
pulse.V1_dur = 0.05; %% duration in seconds of first voltage value in each pulse 

% Make sure that we aren't trying to hold stim chan at Vbase longer than
% the intertrial interval.
if pulse.Vbase_dur_post > locmap_stim.dur_intertrial || pulse.Vbase_dur_post > orientmap_stim.dur_intertrial
    error('Stim channel post-stimulus duration (%gs) exceeds',...
        'stimulus intertrial interval (%gs)',Vbase_dur_post,pars.dur_intertrial)
end

% Make sure we aren't trying to send more than 10 rows or columns or 100
% orientations (stim coding not yet set up for this many.)
if locmap_stim.rows > 10 || locmap_stim.columns > 10
    error('Stim-channel coding not yet set up to handle more than 10 rows or columns.')
elseif orientmap_stim.nAngles > 100
    error('Stim-channel coding not yet set up to handle more than 100 orientations.')
end

    
                      
% Start pyschtoolbox and Data Acquisition Toolbox sessions
% Note: all 3 screens (stim, stim-duplicate, and control) have max resolution 
%    of 1920x1080, where aspect ratio of screen resolution = aspect ratio of 
%    physical screen size.  
input('Press Enter when screen is away from eye or turned off.');    %% avoid exposing eye to psychtoolbox splash
Screen('CloseAll');
stimScreen = max(Screen('Screens'));
[windowPtr, rect] = Screen('OpenWindow',stimScreen,BlackIndex(stimScreen));
nominal_refresh = 1/Screen('NominalFrameRate',windowPtr);
screenstats = Screen('Resolution',stimScreen);

% Convert degrees to pixels. We are assuming that aspect ratio of screen
% resolution = aspect ratio of physical screen size, and therefore using
% vert_res for horz_res and screen_height for screen_width. 
locmap_stim.rect_wh = [deg2pixels(locmap_stim.v_diam,screen_height,screenstats.height,eye2screen_top_bottom),...
                       deg2pixels(locmap_stim.h_diam,screen_height,screenstats.height,eye2screen_top_bottom)]; 
orientmap_stim.rect_wh = [deg2pixels(orientmap_stim.v_diam,screen_height,screenstats.height,eye2screen_top_bottom),...
                          deg2pixels(orientmap_stim.h_diam,screen_height,screenstats.height,eye2screen_top_bottom)]; 

if ifi_check
    [ifi ifi_samp std ] = Screen('GetFlipInterval', windowPtr, 500);
    if abs(ifi - nominal_refresh) > 1e-5
        error(['Interflip interval more than 1e-5 seconds different from nomimal screen refresh rate:'...v 
            sprintf('IFI - Nominal Refresh = %g',(ifi - nominal_refresh))]);
    end
else
    ifi = nominal_refresh;
    warning('Not computing interflip interval, using ifi = nominal refresh = %g instead.', nominal_refresh);
end

daq_sess = daq.createSession('ni');
addAnalogOutputChannel(daq_sess,'Dev1','ao0','Voltage');  %%% use digital/counter output channel instead?

commandwindow
input(sprintf(['\nAssuming that top and bottom of screen are both %g inches from eye,\n',...     %% check parameters
    '   vertical screen resolution is %g,\n   and left and right edges are equidistant from eye.\n',...
    'Press Enter key when screen is in position.'], eye2screen_top_bottom, screenstats.height));    %% check that screen is in position

%% RF Location Mapping
% LocMap Step 1: obtain approximate RF location by hand. 
%   Saved location will be the center of the stimulus shape upon keypress. 
%%%%% should include check/warning if the grid to be presented will be partly off-screen   
disp(sprintf('\nBeginning RF-location mapping...'))
disp(sprintf(['\nFind approximate receptive field location manually by ear.\n',...
    '   Adjust the screen to place the RF center at the screen center.\n',...
    'Press any key when bar is centered on receptive field.']))
manual_rf_center = follow_cursor;

% LocMap Step 2: begin Merec2 recording.
input('Begin Merec2 recording on diesel, then press Enter on stimulus computer.');

% LocMap Step 3: present RF-mapping stimuli.
disp(sprintf('\nNow presenting RF location-mapping stimuli...'));
loc_stimorder = rf_loc_map_stim(windowPtr, rect, ifi, daq_sess, locmap_stim, pulse, screenstats, manual_rf_center); 

% LocMap Step 4: remotely run analysis program on diesel to extract spiking
% responses. 
disp('End the recording on diesel, then run location-mapping-response-analysis program on diesel.');
input('Press Enter when diesel completes RF location analysis.')

% LocMap Step 5: Import spike response information to determine RF center 
[junk sysout] = system(['plink -pw abcd1234 illya@n5.wustl.edu ssh -X -q -t -t illya@diesel head -50 /home/illya/Andrew/recordings/for_stimulus_computer/loc_response_generic 2>&1 2> deletethis'])
sysout = str2num(sysout); % values are impoarted as char
datecheck = clock;
datecheck = datecheck(1:3); % only check year, month, and day
if ~isequal(sysout(2,1:3),datecheck) %% check that generic output file was created today
    warning(['Error: information imported from recording computer was not',...
        'created today and is probably not from this session.']) % exception: performing experiment around midnight
end
%%%%% need to check that that coordinates below are in the format we expect
computed_rf_center = sysout(1,1:2);   %% ? need to check: first is x coordinate, second is y?
%%%% -get this value as (x,y) for CenterRectOnPOint for destinationRect in rf_orient_map_stim


%% Orientation Preference Mapping
% OrientMap Step 2: begin Merec2 recording.
disp('Beginning orientation preference mapping...')
input('Begin Merec2 recording on diesel, then press Enter on stimulus computer.');

% OrientMap Step 3: present orientation preference mapping stimuli.
disp('Now running orientation mapping stimuli...');
orient_stimorder = rf_orient_map_stim(windowPtr, rect, ifi, daq_sess, orientmap_stim, pulse, computed_rf_center); %% last arg will be generated in locmap step 5

% OrientMap Step 4: remotely run analysis program on diesel to extract spiking
% responses. 
disp('End the recording on diesel, then run orientation-preference-response-analysis program on diesel.')
input('Press Enter when diesel completes orientation preference analysis.')

% OrientMap Step 5: Import spike response information to determine preferred
% orientation. 
[junk sysout] = system(['plink -pw abcd1234 illya@n5.wustl.edu ssh -X -q -t -t illya@diesel head -50 /home/illya/Andrew/recordings/for_stimulus_computer/orient_response_generic 2>&1 2> deletethis'])
sysout = str2num(sysout); % values are impoarted as char
datecheck = clock;
datecheck = datecheck(1:3); % only check year, month, and day
if ~isequal(sysout(2,1:3),datecheck)  %% check that generic output file was created today
    warning(['Error: information imported from recording computer was not',...
        'created today and is probably not from this session.']) % exception: performing experiment around midnight
end
%%%% need to check that the angle below is in the same format as we expect
preferred_angle = sysout(1,1);   %% take imported values of the preffered orientation

%%% maybe call this function from the larger program which presents the
%%% experimental stimuli of interest (size tuning vs. sf and tf, etc)

%%%% maybe prompt at some point for other parameters of the experiment to
%%%% save: location of screen, animal ID.... save date and time... other
%%%% notes about animal condition or observations of recording

if exist('old_visualdebug','var')
    Screen('Preference','VisualDebugLevel',old_visualdebug);
end

end

%% deg2pixels: function to convert degrees subtended to pixels on the screen
% All measurements and output are taken in one dimension of interest.
% We are assuming that the stimulus is symmetric and centered on the screen center.
% 
    % deg is desired length of stimulus on screen in degrees
    % screen_length is the length of the screen in the chosen length unit
    % res is the number of pixels on the screen along this dimension (resolution) 
    % eye2screen_edge is the measured distance from the eye to the closest edge of the screen
    %       in the chosen length unit along this dimension
function pixels = deg2pixels(deg,screen_length,res,eye2screen_edge)
    eye2screen_center = sqrt( (screen_length/2)^2 - eye2screen_edge^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
    stim_diam = 2*eye2screen_center*tan(deg2rad(deg)/2); % because (stim_diam/2) / eye2screen_center = tan(deg/2)
    pixels = stim_diam * res/screen_length;  %% convert length to pixels
end