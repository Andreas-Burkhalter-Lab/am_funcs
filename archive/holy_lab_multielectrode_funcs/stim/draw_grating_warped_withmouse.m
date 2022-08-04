function [] = draw_grating_warped_withmouse()
%DRAW_GRATING_WARPED Draw surround suppression stimuli warped to simulate a
%spherical surface on the flat screen. 
% Present gratings warped to resemble a sinusoidal grating on a sphere
% projected onto the 2D screen, as per Marshel et al. 2011
%%% First create outer-grating stim files with generate_ss_stim.m; these files will be
%%% loaded by this script to create stimuli. Duration, Angle, SF, TF, and outer amplitude must
%%% all be specified beforehand in generate_ss_stim. Diameter and location of the outer
%%% grating are controlled by this script by creating the outer aperture. Inner gratings
%%% can be created within this script without overloading memory because there
%%% are only as many as the number of orientations. This script will present the
%%% stimuli contained in stimfile in a randomized order. 
%%%%% last updated 9/15/15 on stim comp

%% selection of dstRect and srcRect is probably wrong; see surround_suppression_stim for correct method

clear
%% Setup
trackmouse = 0; % if on, record mouse movements during stimulus presentation
isi = 0.5; % interstimulus interval in seconds
stim_center = [1300 200]; % projective center for drawing aperture; not the 2D stim center in pixels
diam_inner = 13; % diameter of inner grating in degrees
diam_outer = [80];% 50 70 100 %stim diameter of outer grating in degrees
amp_inner = 1; % inner grating amplitude (determines contrast)
prerender_file = 'stimset_prerender.mat'; % file containing prerendered outer gratings to load 
file_savename = 'ss_stimrec'; % filename for saving stim parameters info

% Create variables for making gratings. 
Screen('CloseAll'); %% clear all previously drawn textures to free up memory
myScreen = max(Screen('Screens'));
[win , winrect] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
screenstats = Screen('Resolution',win);
scrnperpx = round(screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
scrnperpy = round(screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
[xcentermesh ycentermesh] = meshgrid(1:screenstats.width,1:screenstats.height);
perp2stimcenterx = scrnperpx - stim_center(1); % x distance from perpendicular to stimulus center
perp2stimcentery = scrnperpy - stim_center(2); % y distance from perpendicular to stimulus center
perp2meshx = scrnperpx-xcentermesh;
perp2meshy = scrnperpy-ycentermesh;
[xstraight ystraight] = meshgrid(-screenstats.width/2 : screenstats.width/2, -screenstats.height/2 : screenstats.height/2);

% Check that all stimuli fit on screen and shrink if necessary. 
%    edgeVals are values from the specified screen edge listing the aperture
%        radius (in radians) required from stim_center to contain each pixel
matobj = matfile(prerender_file, 'Writable', false); % memory map data from the prerender .m file 
eye2screen_center_pix = matobj.eye2screen_center_pix; % distance from eye to center of screen in pixels
edgeVals.left = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(:,1) + perp2stimcentery*perp2meshy(:,1)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(:,1).^2+perp2meshy(:,1).^2)) );
edgeVals.right = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(:,end) + perp2stimcentery*perp2meshy(:,end)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(:,end).^2+perp2meshy(:,end).^2)) );
edgeVals.top = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(1,:) + perp2stimcentery*perp2meshy(1,:)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(1,:).^2+perp2meshy(1,:).^2)) );
edgeVals.bottom = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(end,:) + perp2stimcentery*perp2meshy(end,:)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(end,:).^2+perp2meshy(end,:).^2)) );
radMaxAllowable = min([edgeVals.left; edgeVals.right; edgeVals.top'; edgeVals.bottom']); % largest radius that can fit on screen

if deg2rad(diam_inner/2) > radMaxAllowable % if inner diam does not fit on screen
    error(sprintf(['Inner grating does not fit on screen. Diameters must be \n    ',...
        num2str(2*rad2deg(radMaxAllowable)) ' degrees or smaller to fit on screen.'])) % convert to diameter radians
end
if deg2rad(max(diam_outer)/2) > radMaxAllowable
    shrink_ok = input(sprintf(['Largest gratings (%g deg) are cut off by screen. ',... % give option to automatically shrink
        'Enter ''y'' to reduce grating diameters larger than',...
        '\n      %g degrees to %g degrees.\n'],...
        max(diam_outer), 2*rad2deg(radMaxAllowable), 2*rad2deg(radMaxAllowable)),'s');
    if strcmp(shrink_ok,'y')
        [junk cutOffGratings] = find(diam_outer > 2*rad2deg(radMaxAllowable));
        diam_outer(cutOffGratings) = 2*rad2deg(radMaxAllowable); % shrink the cut-off gratings to max size that will fit
    else
        error('Will not shrink cut-off gratings.')
    end
end
    
% Construct aperture for the inner grating. 
apt_fullscreen = deg2rad(diam_inner/2) >... % get the area covered by diameters larger than diam_inner/2 to mask out
    acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
[apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
inner_apt_rect = [min(aptx) min(apty) max(aptx) max(apty)];
aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of aperture into alpha; RGB=0
inner_apertex = Screen('MakeTexture', win, aperture);
                
% Construct inner gratings of different angles  for each frame. Use the
% angles, sf_fixed, and tf_fixed already specified by stimset.mat. 
% (Could prerender these along with outer gratings.)
diamtable = table([kron(diam_outer, ones(1,height(matobj.tex_outer)))]','VariableNames',{'diam'}); % outer grating diameters
ingratingtex_table  = table(repmat({NaN(1,matobj.nframes)},height(diamtable),1),'VariableNames',{'inner_gratingtex'}); % inner grating texture pointers
outer_apertex_table = table(NaN(size(diamtable)),'VariableNames',{'outer_apertex'}); % outer aperture texture pointers
outer_aperrect_table = table(cell(size(diamtable)),'VariableNames',{'outer_apt_rect'}); % outer aperture rect values
mouseScanStimOn_table = table(NaN(size(diamtable)),'VariableNames',{'mouseScanStimOn'}); % mouse scan at stim start per trial
mouseScanStimOff_table = table(NaN(size(diamtable)),'VariableNames',{'mouseScanStimOff'}); % mouse scan at stim end per trial
stimrec = [diamtable, repmat(matobj.tex_outer,length(diam_outer),1), outer_apertex_table,... % stimulus record
    outer_aperrect_table, ingratingtex_table, mouseScanStimOn_table, mouseScanStimOff_table];

for Angle = 1:length(matobj.theta); % Make the inner gratings. Gratings move along the yrotate dimension.
    spatial_period = 1/matobj.sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
    spatial_period_rad = deg2rad(spatial_period); % convert to radians
    sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
    tfrad = 2*pi*matobj.tf_fixed; % convert from cycles per second to radians per second
    thetarad = deg2rad(matobj.theta(1,Angle));
    xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
    yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
    eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % angles in radians; pi/2 - acos is faster than asin
    thisAngle_ptrs = inf(1,matobj.nframes);
    for frame = 1:matobj.nframes 
        grating = 0.5*WhiteIndex(win) + 0.5*WhiteIndex(win)*amp_inner*...
            cos(2*pi*sfrad*(eye2screen) - frame*matobj.ifi*tfrad); % from Marshel et al. 2011
        thisAngle_ptrs(frame) = Screen('MakeTexture', win, grating);
    end
    match = stimrec.Angle == matobj.theta(1,Angle);
    stimrec.inner_gratingtex(match) = {thisAngle_ptrs};
end
            
% Construct apertures for the outer grating of sizes listed in diam_outer. 
% (Could prerender these along with outer gratings.)
 for diam = 1:length(diam_outer)
    clear aperture
    
    % R side of inequality lists aperture radius (radians) required from stim_center to contain each pixel
    apt_fullscreen = deg2rad(diam_outer(diam)/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of apt into alpha; RGB=0
    match = stimrec.diam == diam_outer(diam);
    stimrec.outer_apt_rect(match) = {[min(aptx) min(apty) max(aptx) max(apty)]};
    stimrec.outer_apertex(match) = Screen('MakeTexture', win, aperture);
 end
 
 stimrec = stimrec(randperm(height(stimrec)),:); %% shuffle the trial order

%% Mouse Tracking
% need to figure out which screen measurement mouse is on and Hide this
% mouse
global mouse_movements

%%% multiply expected number of scans per trial by sampleBufferFactor when 
%%%   pre-allocating mouse_movements values
sampleBufferFactor = 1.5; 
mouseScanRate = 300; % rate of recording mouse position in hz

totaldur = height(stimrec) * (matobj.dur + isi); % total approximate stim presentation duration in seconds
mouse_movements = NaN(round(sampleBufferFactor * totaldur * mouseScanRate),2);
trial = 0; 
mouseScanIndex = 0; % counter for mouse scans over the course of this experiment
[mousePre(1) mousePre(2)] = GetMouse; % needs two output arguments to give x and y
mouseTimer = timer;
mouseTimer.ExecutionMode = 'fixedSpacing';
mouseTimer.Period = 1/mouseScanRate;
mouseTimer.TimerFcn = @(h,~)trackmouse_callback(h);
mouseTimer.TasksToExecute = round(size(mouse_movements,1)); % stop at this point in case error before stop(mouseTimer)
    
if trackmouse
    start(mouseTimer); % begin measuring mouse movements
end
    
%% Present stimuli
%%%%% use while loop to correctly start after isi rather than using
%%%%% pause... record a timestamp after the final flip around when we send
%%%%% the stimoff daq signal
for trial = 1:height(stimrec)    
    % some flips still miss ifi - check on stim comp, see which steps are
    % most time-consuming
    %%% check isi limit on stim comp
    %%%%%% maybe also making gratings into uint8 to save memory
    
    clear outer_gratings
    outgratrow = stimrec.grating_cell_row(trial); % get the appropriate row from which to load outer gratings
    outer_gratings = matobj.tex_outer_grating(outgratrow,:); % loading these will take time; maybe check that it's not longer than isi
    outer_rect_apt = stimrec.outer_apt_rect{trial}; 
    
    stimrec.mouseScanStimOn(trial) = mouseScanIndex; % mouse scan corresponding to beginning of stim presentation
    for frame = 1:matobj.nframes
        % Outer Grating
        gratingtex = Screen('MakeTexture',win,outer_gratings{frame}); % will all of these texs overload memory? could speed by limiting Make to outer_rect_apt
        Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', win, [0 0 0 0], outer_rect_apt); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', win, stimrec.outer_apertex(trial), [], outer_rect_apt); % set alpha within aperture to 1
        Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        Screen('DrawTexture', win, gratingtex, outer_rect_apt, outer_rect_apt) %%    % draw the outer grating  
        
        % Inner Grating
        Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', win, [0 0 0 0], inner_apt_rect); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', win, inner_apertex, [], inner_apt_rect); % set alpha within aperture to 1
        Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        Screen('DrawTexture', win, stimrec.inner_gratingtex{trial}(frame), inner_apt_rect, inner_apt_rect) % draw the inner grating  
                
        vbl(frame) = Screen('Flip',win); % vbl for checking for missed flips
    end
    Screen('FillRect',win,BlackIndex(win)); Screen('Flip',win);
    stimrec.mouseScanStimOff(trial) = mouseScanIndex; % mouse scan corresponding to end of stim presentation
    
    %%%% better method for isi would be checking toc or using timer before
    %%%% beginning next stim
%     pause(isi);
end

stop(mouseTimer);
stimrec.grating_cell_row = []; % clear unnecessary variables from stimrec
stimrec.outer_apertex = [];
stimrec.outer_apt_rect = [];
stimrec.inner_gratingtex = [];
save(file_savename,'mouseScanRate', 'mouseTimer', 'mouse_movements', 'stimrec', 'isi',...
    'stim_center', 'diam_inner', 'diam_outer', 'amp_inner', 'prerender_file', 'vbl');
clear global mouse_movements
commandwindow;

%% reorganize mouse_movements into a dataset

%% Timer callback function that measures mouse displacement between time intervals
% Uses and modifies variables within the scope of draw_grating_warped. 
function trackmouse_callback(h)
    xResetDist = 300; % recenter the mouse if it gets farther than this many pixels from horizontal center
    yResetDist = 200; % recenter the mouse if it gets farther than this many pixels from vertical center

    if trial == 0; % first trial hasn't start yet
        return % don't start scanning until the first trial starts
    end
    
    mouseScanIndex = mouseScanIndex + 1;
    [xNow yNow] = GetMouse;
    
    % If the mouse moved since the last scan, record this scan number (col 1)
    % within the trial (col 1) and the x (col 2) and y (col 3) displacements. 
    if any([xNow yNow] ~= mousePre); % if there was a mouse displacement
        mouse_movements(mouseScanIndex,1:2) = [xNow yNow]-mousePre;
    end
    
    % Recenter mouse position if necessary and set current mouse position
    % as mousePre for the next scan. 
    if abs(xNow - mouseScreenCenter(1)) > xResetDist || abs(yNow - mouseScreenCenter(2)) > yResetDist
        SetMouse(mouseScreenCenter(1),mouseScreenCenter(2));
        mousePre = mouseScreenCenter;
    else
        mousePre = [xNow yNow];
    end
    
end

end