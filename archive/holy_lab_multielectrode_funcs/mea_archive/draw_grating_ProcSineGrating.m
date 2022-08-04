function [  ] = draw_grating(  )
%DRAW_GRATING Draw simple drifting grating for testing.
% Last modified 5/5/15

dur = 5; %% seconds
backgroundColorOffset = [0 0 0 0]; 
diam = 2; %% diameter of stimuli in degrees; divide by 2 to get radius (approximate for off-center stim)... 2 in Gao
contrastPreMultiplicator = 1; %% defaults to 1
angle = 100;
startphase = 135; 
modulateColor = [255 255 255 0];  %RGB Change to grating color
sfreq = .24; % = grating freq in cyc/deg if support width and height = rect_WH; Gao range ~0.01:1.6, peak ~0.06 cyc/deg
tfreq = .4; % temporal frequency in hz; Gao range ~0.1:13, peak ~4hz
amplitude = 1e38; 

rect_center = [500,500];

ifi_check = 0; %% compare interflip interval to screen refresh rate; must always be on for real experiments
Screen('Preference', 'SkipSyncTests', 1);
debugoriginal = Screen('Preference','VisualDebugLevel');
Screen('Preference', 'VisualDebugLevel', 0);    % set splash screen color to black... comment out for experiments

screen_height = 13;   %% height in inches of stimulus screen; width = 23.5 inches
eye2screen_top_bottom = 9; %% distance in inches from eye to both top and bottom of stimulus screen

if isempty(Screen('Windows'))
    myScreen = max(Screen('Screens'));
    [ win , winrect ] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
else
    win = Screen('Windows');
    win = win(1);
    winrect = Screen('Rect', win);
    Screen('FillRect',win, BlackIndex(win), winrect);
end
nominal_refresh = 1/Screen('NominalFrameRate',win);
screenstats = Screen('Resolution',win);

if ifi_check
    [ifi ifi_samp std ] = Screen('GetFlipInterval', windowPtr, 500);
    if abs(ifi - nominal_refresh) > 1e-5
        error(['Interflip interval more than 1e-5 seconds different from nomimal screen refresh rate:',...
            'IFI - Nominal Refresh = %g'],(ifi - nominal_refresh));
    end
else
    ifi = nominal_refresh;
end

% Make sure grating spatial period is between 0 and 180 degrees. 
if 1/sfreq<0 | 1/sfreq>180
    error('Grating spatial period must be between 0 and 180 degrees (spatial frequency greater than 0.0056).')
end

% Make sure specified temporal frequency is positive and doesn't exceed 0.5/IFI. 
if tfreq<0 | tfreq>0.5/ifi
    error('Stimulus temporal frequency must be between 0 and 0.5/IFI (=%g) or direction will appear reversed.',0.5/ifi)
end

% Convert units (degs to pix, hz to degs/flip, cycles/deg to pix/cycle). We are assuming that 
% aspect ratio of screen resolution = aspect ratio of physical screen size, and therefore 
% assuming that the same vert and horz diameter can be used for a circular stimulus. 
%%% Note: this method isn't quite correct for for any distances not in
%%% center of screen. 
diam_pix = deg2pixels(diam, screen_height, screenstats.height, eye2screen_top_bottom); 
sfreq_pix = 1 / deg2pixels(1/sfreq, screen_height, screenstats.height, eye2screen_top_bottom);% spatial period in pixels
tfreq_degperflip = 360 * ifi * tfreq;
edge_refresh = (tfreq/sfreq_pix) % velocity in pix/s... maybe warn if too low


% Warn if spatial period exceeds stimulus diameter.
if 1/sfreq_pix > 2*diam_pix
    warning('Grating spatial period exceeds 2x stimulus diameter in RF center location mapping.')
end

nframes = dur/ifi;
rect_WH = round([diam_pix diam_pix]); % larg enough not to cut off the aperture

%%%% Make support_width and support_height = rect_WH so that freq = period of grating in pixels.   
support_width = rect_WH(1);
support_height = rect_WH(2);

% rect_center = winrect(3:4)/2;
dstRect = CenterRectOnPoint([0 0 rect_WH], rect_center(1), rect_center(2));

[gratingid, gratingrect] = CreateProceduralSineGrating(win, support_width, support_height,...
    backgroundColorOffset, diam_pix/2, contrastPreMultiplicator); % radius = diam/2
Screen('DrawTexture', win, gratingid, [], dstRect, angle, [], [], modulateColor,...
    [], kPsychDontDoRotation, [startphase sfreq_pix amplitude 0]);
Screen('Flip',win);
    
for frame = 2:nframes
    Screen('DrawTexture', win, gratingid, [], dstRect, angle, [], [], modulateColor,...
        [], kPsychDontDoRotation, [startphase+frame*tfreq_degperflip sfreq_pix amplitude 0]);
    Screen('Flip',win);
end

% Screen('CloseAll')
Screen('MATLABToFront');
Screen('Preference', 'SkipSyncTests', 0);
if exist('debugoriginal','var')
    Screen('Preference', 'VisualDebugLevel', debugoriginal);
end

end

%% deg2pixels: function to convert degrees subtended to pixels on the screen
% All measurements and output are taken in one dimension of interest (length or height).
% We are assuming that the stimulus is symmetric and centered on the screen center.
    % deg is desired length of stimulus on screen in degrees
    % screen_length is the length of the screen in the chosen length unit
    % res is the number of pixels on the screen along this dimension (resolution) 
    % eye2screen_edge is the measured distance from the eye to the closest edge of the screen
    %       in the chosen length unit along this dimension
function pixels = deg2pixels(deg,screen_length,res,eye2screen_edge)
    eye2screen_center = sqrt( eye2screen_edge^2 - (screen_length/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
    stim_diam = 2*eye2screen_center*tan(deg2rad(deg)/2); % because (stim_diam/2) / eye2screen_center = tan(deg/2)
    pixels = stim_diam * res/screen_length;  %% convert length to pixels
end

