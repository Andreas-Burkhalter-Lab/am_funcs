function [] = draw_grating_warped()
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
%%%%% last updated 8/27/15 on msi

%% not yet working - mouse_movements only gets values for first trial and only stores 1 value
%%  (should be 3 - scan number, x displ, y displ
%% also started crashing system

clear
tic
%% Setup
isi = 0.5; % interstimulus interval in seconds
stim_center = [800 400]; % projective center for drawing aperture; not the 2D stim center in pixels
diam_inner = 47; % diameter of inner grating in degrees
diam_outer = [65];% 50 70 100 %stim diameter of outer grating in degrees
amp_inner = 1; % inner grating amplitude (determines contrast)
prerender_file = 'C:\Users\AM\Documents\Lab\matlab_functions\stimset_prerender.mat'; % file containing prerendered outer gratings 
file_savename = 'ss_stimpars'; % filename for saving stim parameters info

% Create variables for making gratings. 
Screen('CloseAll'); %% clear all previously drawn textures to free up memory
myScreen = max(Screen('Screens'));
[win , winrect] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
screenstats = Screen('Resolution',win);
scrnperpx = round(screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
scrnperpy = round(screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
[xcentermesh ycentermesh] = meshgrid(1:screenstats.width,1:screenstats.height);
perp2stimcenterx = scrnperpx - stim_center(1);
perp2stimcentery = scrnperpy - stim_center(2);
perp2meshx = scrnperpx-xcentermesh;
perp2meshy = scrnperpy-ycentermesh;
[xstraight ystraight] = meshgrid(-screenstats.width/2 : screenstats.width/2, -screenstats.height/2 : screenstats.height/2);

% Construct aperture for the inner grating. 
matobj = matfile(prerender_file, 'Writable', false); 
z0_pix = matobj.z0_pix;
apt_fullscreen = deg2rad(diam_inner/2) > acos( (z0_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
    sqrt( (z0_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(z0_pix^2+perp2meshx.^2+perp2meshy.^2)) );
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
stimpars = [diamtable, repmat(matobj.tex_outer,length(diam_outer),1), outer_apertex_table, outer_aperrect_table, ingratingtex_table];

for Angle = 1:length(matobj.theta); % Make the inner gratings. Gratings move along the yrotate dimension.
    sfrad = deg2rad(matobj.sf_fixed);
    thetarad = deg2rad(matobj.theta(1,Angle));
    xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
    yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
    eye2screen = pi/2 - acos(yrotate ./ sqrt(z0_pix^2 + xrotate.^2 + yrotate.^2)); % angles in radians; pi/2 - acos is faster than asin
    thisAngle_ptrs = inf(1,matobj.nframes);
    for frame = 1:matobj.nframes 
        grating = 0.5*WhiteIndex(win) + 0.5*WhiteIndex(win)*amp_inner*...
            cos(2*pi*sfrad*(eye2screen) - frame*matobj.ifi*matobj.tf_fixed); % from Marshel et al. 2011
        thisAngle_ptrs(frame) = Screen('MakeTexture', win, grating);
    end
    match = stimpars.Angle == matobj.theta(1,Angle);
    stimpars.inner_gratingtex(match) = {thisAngle_ptrs};
 end
            
% Construct apertures for the outer grating of sizes listed in diam_outer. 
% (Could prerender these along with outer gratings.)
 for diam = 1:length(diam_outer)
    clear aperture
    apt_fullscreen = deg2rad(diam_outer(diam)/2) > acos( (z0_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (z0_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(z0_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of apt into alpha; RGB=0
    match = stimpars.diam == diam_outer(diam);
    stimpars.outer_apt_rect(match) = {[min(aptx) min(apty) max(aptx) max(apty)]};
    stimpars.outer_apertex(match) = Screen('MakeTexture', win, aperture);
 end
 
 stimpars = stimpars(randperm(height(stimpars)),:); %% shuffle the trial order
 save(file_savename);

 toc
%% Present stimuli
%%%%% use while loop to correctly start after isi rather than using
%%%%% pause... record a timestamp after the final flip around when we send
%%%%% the stimoff daq signal

%%%% Mouse Tracking


% need to figure out which screen measurement mouse is on and Hide this
% mouse
% % % % % global mouse_movements nScansThisTrial trial mouseScreenCenter mousePre
mouseScreenCenter = [(winrect(3)-winrect(1))/2 (winrect(4)-winrect(2))/2]; 
mouse_movements = cell(height(stimpars),1);
mouseScansPerTrial = inf(height(stimpars),1); % number of mouse scans measured per trial
trial = 0; 
nScansThisTrial = 0; % counter for number of scans taken this trial
mousePre = GetMouse;

mouseScanRate = 300; % rate of recording mouse position in hz
mouseTimer = timer;
mouseTimer.ExecutionMode = 'fixedSpacing';
mouseTimer.Period = 1/mouseScanRate;
mouseTimer.TimerFcn = @(h,~)track_mouse_movement(h);
mouseTimer.TasksToExecute = inf; % keep getting mouse position until stop(mouseTimer)

start(mouseTimer); % begin measuring mouse movements
for trial = 1:height(stimpars)    
    % some flips still miss ifi - check on stim comp, see which steps are
    % most time-consuming
    %%% check isi limit on stim comp
    %%%%%% maybe also making gratings into uint8 to save memory
    
    clear outer_gratings
    outgratrow = stimpars.grating_cell_row(trial); % get the appropriate row from which to load outer gratings
    outer_gratings = matobj.tex_outer_grating(outgratrow,:); % loading these will take time; maybe check that it's not longer than isi
    outer_rect_apt = stimpars.outer_apt_rect{trial}; 
    
    for frame = 1:100 %%%%%%%%% matobj.nframes
        % Outer Grating
        gratingtex = Screen('MakeTexture',win,outer_gratings{frame}); % will all of these texs overload memory? could speed by limiting Make to outer_rect_apt
        Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', win, [0 0 0 0], outer_rect_apt); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', win, stimpars.outer_apertex(trial), [], outer_rect_apt); % set alpha within aperture to 1
        Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        Screen('DrawTexture', win, gratingtex, outer_rect_apt, outer_rect_apt) %%    % draw the outer grating  
        
        % Inner Grating
        Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', win, [0 0 0 0], inner_apt_rect); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', win, inner_apertex, [], inner_apt_rect); % set alpha within aperture to 1
        Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        Screen('DrawTexture', win, stimpars.inner_gratingtex{trial}(frame), inner_apt_rect, inner_apt_rect) % draw the inner grating  
                
        vbl(frame) = Screen('Flip',win);
    end
    Screen('FillRect',win,BlackIndex(win)); Screen('Flip',win);
    mouseScansPerTrial(trial) = nScansThisTrial;
    nScansThisTrial = 0; % reset for next trial
end

stop(mouseTimer);
save del
% % % % % % % clear global mouse_movements nScansThisTrial trial mouseScreenCenter mousePre

%% reorganize mouse_movements into a dataset

%% Subfunction for the timer that calls track_mouse_movement
% Uses and modifies variables within the scope of draw_grating_warped. 
function track_mouse_movement(h)
    xResetDist = 300; % recenter the mouse if it gets farther than this many pixels from horizontal center
    yResetDist = 200; % recenter the mouse if it gets farther than this many pixels from vertical center

    if trial == 0; % first trial hasn't start yet
        return % don't start scanning until the first trial starts
    end
    
    nScansThisTrial = nScansThisTrial + 1;
    [xNow yNow] = GetMouse;
    
    % If the mouse moved since the last scan, record this scan number (col 1)
    % within the trial (col 1) and the x (col 2) and y (col 3) displacements. 
    if any([xNow yNow] ~= mousePre); % if there was a mouse displacement
        mouse_movements{trial} = [mouse_movements{trial},...
            nScansThisTrial mouse_movements{trial} [xNow yNow]-mousePre]; % record the displacement
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