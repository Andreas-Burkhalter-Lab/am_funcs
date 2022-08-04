function [computed_rf_center preferred_angle] = rf_mapping(varargin)
%%% Last updated 5/16/2015 in Holy Lab
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
global screen_height eye2screen_top_bottom winPtr winrect screenstats ifi daq_sess

ifi_check = 0; %% compare interflip interval to screen refresh rate; must always be on for real experiments

% get rid of warnings and splash screen... comment out for real experiments
old_visualdebug = Screen('Preference', 'VisualDebugLevel');
Screen('Preference', 'VisualDebugLevel', 0);    

% RF Location-Mapping Stimuli Parameters
%%%%%% maybe give option to show each stimulus more than once (shouldn't
%%%%%% require extra code on the analysis side if pulse code is consistent 
% (Stim channel noise (max-min) = ~0.3V over 5s... we can safely accomodate ~100 stimuli using 2-pulse IDs.)
locmapfile = 'deletethis'; %strcat('locmapping_',date,'_sub1'); %% maybe customize more for each session/subject
locmap_stim.iterations = 2; %% number of times to show the entire grid of stimuli
locmap_stim.rows = 3;       %% number of rows of stimuli in the grid
locmap_stim.columns = 3;    %% number of columns of stimuli in the grid
locmap_stim.h_spacing = 10; %% horizontal distance in degrees between stim centers
locmap_stim.v_spacing = 10; %% vertical distance in degrees between stim centers
locmap_stim.diam = 2; %% diameter of stimuli in degrees; divide by 2 to get radius (approximate for off-center stim)... 5 in Gao
locmap_stim.dur_stim = 3; %% duration of presentation of each unique stimulus
locmap_stim.dur_intertrial = 3; %% time in seconds between end of ones stimulus and beginning of next stimulus
locmap_stim.tfreq = 13;   %% temporal frequency in hz; peak in Gao ~4hz
locmap_stim.sfreq = .06; % = grating frequency in cycles/degree if support width and height = rect_wh... ~0.06 peak in Gao
locmap_stim.Angle = 170;            %Degrees
locmap_stim.modulateColor = [255 255 255 0];  %RGB Change to grating color
locmap_stim.backgroundColorOffset = [0 0 0 0];   
locmap_stim.startphase = 0;   %Degrees
locmap_stim.amp = 1e38; % amplitude; grating becomes invisible if this value =>1e(38.6)
locmap_stim.contrastPreMultiplicator = 1;  
locmap_stim.background_color = [0 0 0];

% Orientation Preference-Mapping Stimuli Parameters
orientmap_stim.iterations = 1;      % number of times to show the entire set of orientations
orientmap_stim.nAngles = 16;    %% number of orientations 
orientmap_stim.angle_range = [0  360*(1-1/orientmap_stim.nAngles)];  %% upper and lower limit on angles to test.... maybe add a manual angle-mapping function like with rf center hand-mapping, then limit this range accordingly
orientmap_stim.diam = 2; %% diameter of stimuli in degrees (divide by 2 to get radius) - Gao et al. used 5 degrees
orientmap_stim.dur_stim = 2; %% duration of presentation of each unique stimulus
orientmap_stim.dur_intertrial = 2; %% time in seconds between end of one stimulus/iteration and beginning of next stimulus/iteration
orientmap_stim.modulateColor = [255 255 255 0];  %RGB Change to grating color
orientmap_stim.backgroundColorOffset = [0 0 0 0];   
orientmap_stim.startphase = 0;   %Degrees
orientmap_stim.tfreq = 4;  %% temporal frequency in hz; peak in Gao ~4hz
orientmap_stim.sfreq = .06;  % = grating frequency in cycles/degree if support width and height = rect_wh... ~0.06 peak in Gao
orientmap_stim.amp =  1e38; % amplitude; grating becomes invisible if this value =>1e(38.6)
orientmap_stim.contrastPreMultiplicator = 1;  
orientmap_stim.background_color = [0 0 0];

% Stim channel parameters
% These values are used for both location mapping and orientation preference mapping.  
pulse.Vbase = 0;  %% stim channel value when no stimulus is present

%%%%%%%%%Vbase_dur_post is really the isi; it should be declared as part of locmap_stim 
pulse.Vbase_dur_post = 1; %% length of time in seconds during which stim chan is held at Vbase after stimulus ends
pulse.V1_dur = 0.5; %% duration in seconds of first voltage value in each pulse 

% Check whether or not the locmap log file already exists. 
if exist(locmapfile,'file') | exist(strcat(locmapfile,'.mat'),'file')
    overwrite = input(sprintf('Mapping stimulus log file ''%s'' already exists. Enter ''y'' to overwrite.', locmapfile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
end

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

% Make sure grating spatial period is between 0 and 180 degrees. 
if 1/locmap_stim.sfreq<0 | 1/locmap_stim.sfreq>180 | ...
   1/orientmap_stim.sfreq<0 | 1/orientmap_stim.sfreq>180
    error('Grating spatial period must be between 0 and 180 degrees (spatial frequency greater than 0.0056).')
end
                      
% Start pyschtoolbox and Data Acquisition Toolbox sessions
% Note: all 3 screens (stim, stim-duplicate, and control) have max resolution 
%    of 1920x1080, where aspect ratio of screen resolution = aspect ratio of 
%    physical screen size.  
input('Press Enter when screen is away from eye or turned off.');    %% avoid exposing eye to psychtoolbox splash
Screen('CloseAll');
stimScreen = max(Screen('Screens'));
[winPtr, winrect] = Screen('OpenWindow',stimScreen,BlackIndex(stimScreen));  
nominal_refresh = 1/Screen('NominalFrameRate',winPtr);
screenstats = Screen('Resolution',stimScreen);

if ifi_check
    [ifi ifi_samp std ] = Screen('GetFlipInterval', winPtr, 500);
    if abs(ifi - nominal_refresh) > 1e-5
        error(['Interflip interval more than 1e-5 seconds different from nomimal screen refresh rate:',...
            'IFI - Nominal Refresh = %g'],(ifi - nominal_refresh));
    end
else
    ifi = nominal_refresh;
    warning('Not computing interflip interval, using ifi = nominal refresh = %g instead.', nominal_refresh);
end

% Make sure specified temporal frequency is positive and doesn't exceed 0.5/IFI. 
if locmap_stim.tfreq<0 | locmap_stim.tfreq>0.5/ifi |...
        orientmap_stim.tfreq<0 | orientmap_stim.tfreq>0.5/ifi
    error('Stimulus temporal frequency must be between 0 and 0.5/IFI (=%g) or direction will appear reversed.',0.5/ifi)
end

daq_sess = daq.createSession('ni');
addAnalogOutputChannel(daq_sess,'Dev1','ao0','Voltage');  %%% use digital/counter output channel instead?

% Convert units (degs to pix, hz to degs/flip, cycles/deg to pix/cycle). We are assuming that 
% aspect ratio of screen resolution = aspect ratio of physical screen size, and therefore 
% assuming that the same vert and horz diameter can be used for a circular stimulus. 
%%% Note: deg2pixels is only valid for lengths with one end at the screen center.  
locmap_stim.diam_pix = deg2pixels(locmap_stim.diam); 
locmap_stim.h_spacing_pix = deg2pixels(locmap_stim.h_spacing);
locmap_stim.v_spacing_pix = deg2pixels(locmap_stim.v_spacing);
locmap_stim.sfreq_pix = 1 / deg2pixels(1/locmap_stim.sfreq); % spatial period in pixels
locmap_stim.tfreq_degperflip = 360 * ifi * locmap_stim.tfreq; %%%% deg along a sine wave oscillation, not deg in visual space
locmap_stim.edge_refresh = (locmap_stim.tfreq/locmap_stim.sfreq_pix); % velocity in pixels per second... maybe warn if too low

orientmap_stim.diam_pix = deg2pixels(orientmap_stim.diam); 
orientmap_stim.sfreq_pix = 1 / deg2pixels(1/orientmap_stim.sfreq);% spatial period in pixels
orientmap_stim.tfreq_degperflip =360 * ifi * orientmap_stim.tfreq;
orientmap_stim.edge_refresh = (orientmap_stim.tfreq/orientmap_stim.sfreq_pix); % velocity in pixels per second... maybe warn if too low

% Make sure that spatial period is at least 2 pixels large (should actually be much larger).  
if 1/locmap_stim.sfreq_pix<2*screen_height/screenstats.height | 1/orientmap_stim.sfreq_pix<2*screen_height/screenstats.height
    error('Grating spatial period must be at least 2 pixels wide.')
end

% Warn if spatial period exceeds stimulus diameter.
%%% not really an issue if amp is reasonably small... 
if 1/locmap_stim.sfreq_pix > 2*locmap_stim.diam_pix
    warning('Grating spatial period exceeds 2x stimulus diameter in RF center location mapping.')
end
if 1/orientmap_stim.sfreq_pix > orientmap_stim.diam_pix
    warning('Grating spatial period exceeds stimulus diameter in orientation preference mapping.')
end
    
commandwindow;
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
%%%% need to move daq session initialization here to make sure that stim starts at Vbase or else analyze_rf_responses will get confused 
input('Begin Merec2 recording on diesel, then press Enter on stimulus computer.');

% LocMap Step 3: present RF-mapping stimuli.
disp(sprintf('\nNow presenting RF location-mapping stimuli...'));
loc_stimorder = rf_loc_map_stim(locmap_stim, pulse, manual_rf_center, locmapfile); 

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

%% Orientation Preference Mapping
% OrientMap Step 2: begin Merec2 recording.
disp('Beginning orientation preference mapping...')
input('Begin Merec2 recording on diesel, then press Enter on stimulus computer.');

% OrientMap Step 3: present orientation preference mapping stimuli.
disp('Now running orientation mapping stimuli...');
orient_stimorder = rf_orient_map_stim(winPtr, winrect, ifi, daq_sess, orientmap_stim, pulse, computed_rf_center); %% last arg will be generated in locmap step 5

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
    % deg is desired length of stimulus on screen in degrees (can be a vector)   
    % screen_length is the length of the screen in the chosen length unit
    % res is the number of pixels on the screen along this dimension (resolution) 
    % eye2screen_edge is the measured distance from the eye to the closest edge of the screen
    %       in the chosen length unit along this dimension
function pixels = deg2pixels(deg)
    global screen_height eye2screen_top_bottom screenstats
    eye2screen_center = sqrt( eye2screen_top_bottom^2 - (screen_height/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
    stim_diam = 2*eye2screen_center*tan(deg2rad(deg)/2); % because (stim_diam/2) / eye2screen_center = tan(deg/2)
    pixels = stim_diam * screenstats.height/screen_height;  %% convert length to pixels
end