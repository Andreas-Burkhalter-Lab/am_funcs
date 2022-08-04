function [rf_center] = rfstim_main(rf_stim, pulsepars, defaultpars, commonvars)
%%% Last updated 1/15/16 on stim comp
%%% rf_stim and defaultpars structs are provided by run_experiment with parameters
%%%    determining the characteristics of rf mapping stimulus presentation

%RFSTIM_MAIN Map receptive field
%   First get approximate RF location by manually moving stimulus while
%   listening for audio responses. Once an approximate center is found
%   manually, automatically present stimuli in a grid to more precisely
%   find the RF location around the approximate center by analyzing spiking
%   responses. Preferred orientation and other parameters can be
%   subsequently mapped via the same procedures.

%% Setup
checkRfSideLength = 50; % side length of square displayed when checking computed rf center
checkRfColor = 100;  % color of square displayed when checking computed rf center
singleVStimOn = 0; % multiple VStimOn values will be used to indicate grid location
winrect = commonvars.winrect;
rf_stim = copyAllFields(rf_stim,defaultpars); % use default parameters for all stimulus variables
global daq_sess

% Make sure we aren't trying to send more than 9 rows or columns 
% (stim coding not yet set up for this many).
if rf_stim.rows > 9 || rf_stim.columns > 9
    error('Stim-channel coding not yet set up to handle more than 10 rows or columns.')
end

% Convert units (degs to pix, hz to degs/flip, cycles/deg to pix/cycle). We are assuming that 
% aspect ratio of screen resolution = aspect ratio of physical screen size, and therefore 
% assuming that the same vert and horz diameter can be used for a circular stimulus. 
%%% Note: deg2pixels is only valid for lengths with one end at the screen center.  
rf_stim.diam_pix = deg2pixels(rf_stim.diam, commonvars); 
rf_stim.h_spacing_pix = deg2pixels(rf_stim.h_spacing, commonvars);
rf_stim.v_spacing_pix = deg2pixels(rf_stim.v_spacing, commonvars);
rf_stim.sf_pix = 1 / deg2pixels(1/rf_stim.sf, commonvars); % spatial period in pixels
rf_stim.tf_degperflip = 360 * commonvars.ifi * rf_stim.tf; %%%% deg along a sine wave oscillation, not deg in visual space
rf_stim.edge_refresh = (rf_stim.tf/rf_stim.sf_pix); % velocity in pixels per second... maybe warn if too low

% run stim checks; not using stimrec
[~,pulsepars] = tuningstim_preparestim(rf_stim,pulsepars,commonvars,daq_sess,singleVStimOn);

% Check that grid is not too large to fit on the screen.
grid_width = rf_stim.h_spacing_pix*(rf_stim.columns - 1) + rf_stim.diam_pix; % pixels
grid_height = rf_stim.v_spacing_pix*(rf_stim.rows - 1) + rf_stim.diam_pix; % pixels
reduce_h_spacing = []; reduce_v_spacing = []; % initialize
if grid_width > winrect(3) - winrect(1) % if grid is too wide
    max_h_spacing = floor([winrect(3)-winrect(1)-rf_stim.diam_pix]/rf_stim.columns); % pixels
    reduce_h_spacing = input(sprintf(['Stimulus grid is too wide to fit on screen.\n'...
        'Enter ''y'' to reduce horizontal center-to-center stimulus spacing from %g pixels to %g pixels.\n'],...
        rf_stim.h_spacing_pix, max_h_spacing));
    if strcmp(reduce_h_spacing,'y')
        rf_stim.h_spacing_pix = max_h_spacing;
    else
        error('Not reducing horizontal center-to-center stimulus spacing. Quitting rf_mapping...')
    end
end
if grid_height > winrect(4) - winrect(2); % if grid is too tall
    max_v_spacing = floor([winrect(4)-winrect(2)-rf_stim.diam_pix]/rf_stim.rows); % pixels
    reduce_v_spacing = input(sprintf(['Stimulus grid is too tall to fit on screen.\n'...
        '   Enter ''y'' to reduce vertical center-to-center stimulus spacing from %g pixels to %g pixels.\n'],...
        rf_stim.v_spacing_pix, max_v_spacing));
    if strcmp(reduce_v_spacing,'y')
        rf_stim.v_spacing_pix = max_v_spacing;
    else
        error('Not reducing vertical center-to-center stimulus spacing. Quitting rf_mapping...')
    end
        
end
        
commandwindow;
disp(sprintf(['\nAssuming that top and bottom of screen are both %g inches from eye,\n',...     %% check parameters
    '   vertical screen resolution is %g,\n   and left and right edges are equidistant from eye.\n'],...
    commonvars.eye2screen_top_bottom, commonvars.screenstats.height));   

%% RF Location Mapping
% obtain approximate RF location by hand. 
%   Saved location will be the center of the stimulus shape upon keypress.  
% % % % % disp(sprintf('\nBeginning RF-location mapping...'))
if strcmp(reduce_h_spacing,'y') & strcmp(reduce_v_spacing,'y') % if both hrz and vert spacing had to be reduced
    disp('Stimulus grid will be centered at screen center; skipping manual grid centering.')
else
    disp(sprintf(['\nFind approximate receptive field location manually by ear.\n',...
        '   Adjust the screen to place the RF center at the screen center.\n',...
        '   Use <-- and --> to rotate.\n',...
        '   Use up and down arrows to change bar length.\n',...
        'Press Shift when bar is centered on receptive field.']))
    manual_rf_center = follow_cursor;
    if strcmp(reduce_h_spacing,'y')
        disp('\nGrid will be horizontally centered at screen horizontal center.');
        manual_rf_center(1) = round(winrect(3) - winrect(1));
    elseif strcmp(reduce_v_spacing,'y')
        disp('\n Grid will be vertically centered at screen vertical center.');
        manual_rf_center(2) = round(winrect(4) - winrect(2));
    end
end

%  present RF-mapping stimuli.
rf_stim.warping = 'nonwarped'; % indicate nonwarped stimuli
[loc_stimorder stimcenters] = rfstim_stimloop(rf_stim, pulsepars, commonvars, manual_rf_center); 

%  remotely run analysis program on diesel to extract spiking responses. 
disp('End the recording on diesel, then run location-mapping-response-analysis program on diesel.');

if rf_stim.getRfCenterRemotely %% not currently functional
    input('Press Enter when diesel completes RF location analysis.')
    % Import spike response information to determine RF center 
    [junk sysout] = system(['plink -pw abcd1234 illya@n5.wustl.edu ssh -X -q -t -t illya@diesel head -50 /home/illya/Andrew/recordings/for_stimulus_computer/loc_response_generic 2>&1 2> deletethis'])
    sysout = str2num(sysout); % values are imported as char
    datecheck = clock;
    datecheck = datecheck(1:3); % only check year, month, and day
    if ~isequal(sysout(2,1:3),datecheck) %% check that generic output file was created today
        warning(['Error: information imported from recording computer was not',...
            'created today and is probably not from this session.']) % exception: performing experiment around midnight
    end
    %%%%% need to check that that coordinates below are in the format we expect
    rf_center_gridrelative = sysout(1,1:2);   %% ? need to check: first is x coordinate, second is y?
    grid_width = stimcenters{1,end}(1)-stimcenters{1,1}(1); % leftmost stim center to rightmost stim center in pixels
    grid_height = stimcenters{end,1}(2)-stimcenters{1,1}(2); % topmost stim centfrer to bottommost stim center in pixels
    center2center_x = grid_width/(rf_stim.columns-1); % hrz distance in pixels from one stim center to the next
    center2center_y = grid_height/(rf_stim.rows-1); % vert distance in pixels from one stim center to the next
    rf_center = stimcenters{1,1} + rf_center_gridrelative.*[center2center_x center2center_y]; % rf center in screen coordinates in pixels
elseif ~rf_stim.getRfCenterRemotely % have user enter xypeaks from diesel
    xypeaks(1) = input('Enter xypeaks(1) calculated on recording computer (grid column index):\n'); 
    xypeaks(2) = input('Enter xypeaks(2) calculated on recording computer (grid row index):\n');
    rf_center = grid2screen(stimcenters,xypeaks(1),xypeaks(2)); % convert grid index to screen pixel coordinates
    fprintf('Computed rf_center is x,y = [%g %g].\n',rf_center(1),rf_center(2))
    displayRfCenter = input('Enter ''y'' to display rf_center location on screen. Otherwise press Enter.','s');
    if strcmp(displayRfCenter,'y')
        rect = CenterRectOnPoint([0 0 checkRfSideLength checkRfSideLength],rf_center(1),rf_center(2));
        Screen('FillRect',commonvars.winPtr,checkRfColor,rect); % display gray square centered on computere rf center
        Screen('Flip',commonvars.winPtr);
        input('Press Enter to continue to next stim set.')
        commandwindow;
        Screen('FillRect',commonvars.winPtr,BlackIndex(commonvars.winPtr)); % remove the rf center-checking squqare
        Screen('Flip',commonvars.winPtr);
    end
end

clear global 