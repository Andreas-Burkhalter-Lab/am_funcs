function [] = draw_grating()
%DRAW_GRATING Draw two concentric black/white drifting (nonwarped) sine gratings of 
% different parameters separated by a ring.
%%% last updated 8/16/17

repetitions = 2;

dur = 4; %% seconds
rect_center = [950 540];
% rect_center = follow_cursor;
% Screen('CloseAll')

%%% Webb fig 2 - spacing didn't 'necessarily' exist
spacing_pix = [0 0]; % [hrz vert] spacing in pixels between inner and outer grating 

% First column = inner grating, second column = outer grating
%%% maybe add parameters for color, amplitude (contrast), start phase
Axes = [0 0; 150 100]; %% diameter of [hrz vert] ellipse axes in degrees; (approximate for off-center stim)... 2 each in Gao... screen height ~92.5deg
Angle = [100; 20];
sfreq = [0.01; 0.01]; % cycles per degree; Gao range ~0.01:1.6, peak ~0.06 cyc/deg
tfreq = [.12; 0.1]; % hz

ifi_check = 0; %% compare interflip interval to screen refresh rate; must always be on for real experiments
Screen('Preference', 'SkipSyncTests', 1);
debugoriginal = Screen('Preference','VisualDebugLevel');
Screen('Preference', 'VisualDebugLevel', 0);    % set splash screen color to black... comment out for experiments

screen_height = 13;   %% height in inches of stimulus screen;
eye2screen_top_bottom = 9; %% distance in inches from eye to both top and bottom of stimulus screen

%% Start psychtoolbox and check parameters.
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
gray = ceil((WhiteIndex(win) - BlackIndex(win))/2);





winmain = Screen(1,'OpenWindow',BlackIndex(1));



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
% aspect ratio of screen resolution = aspect ratio of physical screen size (square pixels), 
% and therefore that we can use the same the same conversion for the horz and vert axes.  
%%% Note: this method isn't quite correct for for any distances not in
%%% center of screen. 
Axes_pix = deg2pixels(Axes, screen_height, screenstats.height, eye2screen_top_bottom); 
pix_per_cycle = deg2pixels(1./sfreq, screen_height, screenstats.height, eye2screen_top_bottom);
sfreq_pix = 1 ./ pix_per_cycle; % spatial period in pixels
tfreq_degperflip = 360 * ifi * tfreq; %% do we need this?
shiftperframe = tfreq .* pix_per_cycle * ifi; % pix/frame = cycles/s * pix/cycle * s/frame 
edge_refresh = (tfreq/sfreq_pix); % velocity in pix/s... maybe warn if too low
nframes = dur/ifi; % convert duration from seconds to flips

% Warn if spatial period exceeds largest axis.   
if 1/sfreq_pix > 2*max(max(Axes_pix))
    warning('Grating spatial period exceeds 2x largest stimulus axis in RF center location mapping.')
end




for patch = 1:2
    [x ~] = meshgrid(1:Axes_pix(patch,1:2) + pix_per_cycle(patch), 1);
    grating = 0.5*WhiteIndex(win) + 0.5*WhiteIndex(win)*cos(2*pi*sfreq(patch)*x);
    rect_aperture(patch,:) = CenterRectOnPoint([0 0 Axes_pix(patch,1:2)], rect_center(1), rect_center(2)); 
    gratingtex(patch) = Screen('MakeTexture', win, grating);
    diag = sqrt(Axes_pix(patch,1)^2 + Axes_pix(patch,2)^2); % make grating at least this big to cover the aperture
    rect_clear(patch,:) = CenterRectOnPoint([0 0 sqrt(2)*[diag diag]], rect_center(1), rect_center(2));
    rect_grating(patch,:) = CenterRectOnPoint([0 0 diag diag], rect_center(1), rect_center(2));
end 

rect_spacing = CenterRectOnPoint([0 0 Axes_pix(1,:)+spacing_pix], rect_center(1), rect_center(2));

for rept = 1:repetitions
%% Present Stimuli
for frame = 1:nframes
    xoffset = mod(frame*shiftperframe, pix_per_cycle);
    srcRect = [xoffset [0; 0] xoffset+Axes_pix(:,1) Axes_pix(:,2)];
    
    %%% Outer Grating   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Disable alpha-blending, restrict following drawing to alpha channel:    
    Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]);   

        Screen('Blendfunction', winmain, GL_ONE, GL_ZERO, [0 0 0 1]);   
    
    % Clear 'dstRect' region of framebuffers alpha channel to zero:
    Screen('FillRect', win, [0 0 0 0], rect_clear(2,:)); % maybe take smaller of winrect and this rect

        Screen('FillRect', winmain, [0 0 0 0], rect_clear(2,:)); % maybe take smaller of winrect and this rect
    
    % Fill oval-shaped 'dstRect' region with an alpha value of 255:
    Screen('FillOval', win, [0 0 0 255], rect_aperture(2,:));
    
        Screen('FillOval', winmain, [0 0 0 255], rect_aperture(2,:));
    
    % Enable DeSTination alpha blending and reenable drawing to all
    % color channels. Following drawing commands will only draw there
    % the alpha value in the framebuffer is greater than zero, ie., in
    % our case, inside the circular 'dstRect' aperture where alpha has
    % been set to 255 by our 'FillOval' command:
    Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
    
        Screen('Blendfunction', winmain, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
    
    Screen('DrawTexture', win, gratingtex(2), srcRect(2,:), rect_grating(2,:), Angle(2));
    
        Screen('DrawTexture', winmain, gratingtex(2), srcRect(2,:), rect_grating(2,:), Angle(2));
    
    %%%% Spacer  %%%%%%%%%%%%%%%%%%%%%%%%%%
    Screen('FillOval', win, BlackIndex(win), rect_spacing); % add the spacer
    
    Screen('FillOval', winmain, BlackIndex(win), rect_spacing); % add the spacer
    
    
    %%% Inner Grating %%%%%%%%%%%%%%%%%%%%%%
    % Disable alpha-blending, restrict following drawing to alpha channel:    
    Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]);
    
    % Clear 'dstRect' region of framebuffers alpha channel to zero:
    Screen('FillRect', win, [0 0 0 0], rect_clear(1,:)); % maybe take smaller of winrect and this rect
    
        Screen('FillRect', winmain, [0 0 0 0], rect_clear(1,:)); % maybe take smaller of winrect and this rect
    
    % Fill oval-shaped 'dstRect' region with an alpha value of 255:
    Screen('FillOval', win, [0 0 0 255], rect_aperture(1,:));
    
        Screen('FillOval', winmain, [0 0 0 255], rect_aperture(1,:));

    % Enable DeSTination alpha blending and reenable drawing to all
    % color channels. Following drawing commands will only draw there
    % the alpha value in the framebuffer is greater than zero, ie., in
    % our case, inside the circular 'dstRect' aperture where alpha has
    % been set to 255 by our 'FillOval' command:
    Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
    
        Screen('Blendfunction', winmain, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
        
    Screen('DrawTexture', win, gratingtex(1), srcRect(1,:), rect_grating(1,:), Angle(1));   % draw the grating
    
        Screen('DrawTexture', winmain, gratingtex(1), srcRect(1,:), rect_grating(1,:), Angle(1));   % draw the grating
    
    vbl = Screen('Flip',win);
    
        vbl = Screen('Flip',winmain);
end

pause(1)

end

Screen('Close',win);

Screen('Close',winmain);

commandwindow;

end

%% deg2pixels: function to convert degrees subtended to pixels on the screen
% All measurements and output are taken in one dimension of interest (length or height).
% We are assuming that the stimulus is symmetric and centered on the screen center.
    % deg is desired length of stimulus on screen in degrees (can be a vector)   
    % screen_length is the length of the screen in the chosen length unit
    % res is the number of pixels on the screen along this dimension (resolution) 
    % eye2screen_edge is the measured distance from the eye to the closest edge of the screen
    %       in the chosen length unit along this dimension
function pixels = deg2pixels(deg,screen_length,res,eye2screen_edge)
    eye2screen_center = sqrt( eye2screen_edge^2 - (screen_length/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
    stim_diam = 2*eye2screen_center*tan(deg2rad(deg)/2); % because (stim_diam/2) / eye2screen_center = tan(deg/2)
    pixels = stim_diam * res/screen_length;  %% convert length to pixels
end