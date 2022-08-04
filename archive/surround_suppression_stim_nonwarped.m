function [] = surround_suppression_stim_nonwarpedq(pars, stim_center, pulse, savefile, screenstats)
%SURROUNG_SUPPRESSION_STIM Experimental stimuli for testing the
%relationship between surround suppression and temporal and spatial frequency. 
%   Call after rf_mapping from run_experiment to provide arguments and globals. 
%%% pars = struct of stimulus parameters ('stimpars' from run_experiment.m)
%%% stim_center = screen location on which to center the stimuli; use the 
%%%     output from rf_mapping (default is screen center)
% last updated 12/26/17 on msi

% need to run check that stimuli will not be cut off by the edge of the
% screen

%% Setup
    
commandline_output = 1;

if exist(savefile,'file') | exist(strcat(savefile,'.mat'),'file')
    overwrite = input(sprintf('File ''%s'' already exists. Enter ''y'' to overwrite.\n', savefile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
end

debugoriginal = Screen('Preference','VisualDebugLevel');
if pars.show_splash_screen
    Screen('Preference', 'VisualDebugLevel', 4); % more thorough but show white splash screen  
elseif ~pars.show_splash_screen
    Screen('Preference','VisualDebugLevel',0);   % turn off white splash screen and error warning screen
end
    
if ~pars.do_sync_tests
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference','SkypSyncTests',0);
end

if isempty(Screen('Windows'))
%     myScreen = max(Screen('Screens'));
    myScreen = 2;
    [ screenstats.winPtr , winrect ] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
else
    screenstats.winPtr = Screen('Windows');
    screenstats.winPtr = screenstats.winPtr(1);
end
screenstats = copyAllFields(screenstats, Screen('Resolution',screenstats.winPtr));
screenstats.ifi = 1/screenstats.hz;

% % % % % % % % % Create stim chan signals.
if pars.send_pulse
    daq_sess = daq.createSession('ni');
    chan = addAnalogOutputChannel(daq_sess,'Dev2','ao0', 'Voltage');
    pulse.stimon_dur_scans = pars.stimdur * daq_sess.Rate;    %% get stim-on duration in scans for the stim channel pulse
    pulse.stimoff_dur_scans = pars.isi * daq_sess.Rate; %% get stim-off duration in scans for the stim channel pulse
    pulse.stimon_signal = pulse.Vstimon * ones(pulse.stimon_dur_scans,1); % signal to send at stim onset
    pulse.stimoff_signal = pulse.Vstimoff * ones(pulse.stimoff_dur_scans,1); % signal to send at stim offsetk to Vbase
else
    daq_sess = struct;
end
    
% Check that all stimuli fit on screen and shrink if necessary. 
if pars.do_stimcheck
    stimcheck.diamTried_pix = deg2pixels(max(pars.diam_minmax),screenstats); % max diam specified in run_experiment in pixels
    stimcheck.edge_left = stim_center(1) - stimcheck.diamTried_pix/2; % x coordinate of the left edge of the largest grating
    stimcheck.edge_right = stim_center(1) + stimcheck.diamTried_pix/2;
    stimcheck.edge_top = stim_center(2) - stimcheck.diamTried_pix/2; % y coordinate of the top edge of the largest grating
    stimcheck.edge_bottom = stim_center(2) + stimcheck.diamTried_pix/2;
    stimcheck.pix_cut_off = max([... % positive values indicate that part of at least one grating is off of the screen
        -stimcheck.edge_left + winrect(1),...
        stimcheck.edge_right - winrect(3),...
        -stimcheck.edge_top + winrect(2),...
        stimcheck.edge_bottom - winrect(4)]);
    if stimcheck.pix_cut_off > 0 % if we need to shrink the gratings to fit on the stimulus screen
        stimcheck.diamMaxAllowable_pix = floor(stimcheck.diamTried_pix - 2*stimcheck.pix_cut_off); % max diam that will fit on screen
        stimcheck.diamMaxAllowable = pixels2deg(stimcheck.diamMaxAllowable_pix,screenstats); % largest diam that will fit in degrees
        if min(pars.diam_minmax) >= stimcheck.diamMaxAllowable % if smallest grating is cut off by screen
            error(['Smallest grating diameter (%g pixels, %g deg) is cut off by screen.',...
                '\nGratings must have diameter %g pixels (%g deg) or less to fit on screen.'],...
                deg2pixels(min(pars.diam_minmax),screenstats), min(pars.diam_minmax),...
                stimcheck.diamMaxAllowable_pix, stimcheck.diamMaxAllowable)
        end
        stimcheck.shrink_ok = input(sprintf(['Largest gratings are cut off by screen. ',...
            'Enter ''y'' to reduce maximum grating diameter from',...
            '\n      %g pixels (%g deg) to %g pixels (%g deg).\n'],...
            stimcheck.diamTried_pix, max(pars.diam_minmax), stimcheck.diamMaxAllowable_pix, stimcheck.diamMaxAllowable),'s');
        if strcmp(stimcheck.shrink_ok,'y')
            pars.diam_minmax(2) = stimcheck.diamMaxAllowable;
        else
            error('Will not shrink cut-off gratings.')
        end
    end
end
    
% Create the distributions of annulus sf, tf, sizes, and angles.
% Sf, tf, and size are on log-scales including highest value and angle is on
% a linear scale not including highest value. 
%%% maybe Vaiceliunaite-type distribution of diams is better than exponential; seems
%%% to capture steeper right-side slope which Gao doesn't see
log_sf_minmax = log(pars.sf_minmax);
log_tf_minmax = log(pars.tf_minmax);
log_diam_minmax = log(pars.diam_minmax);
log_sf_range = linspace(log_sf_minmax(1),log_sf_minmax(2),pars.n_sfs);
log_tf_range = linspace(log_tf_minmax(1),log_tf_minmax(2),pars.n_tfs);
log_diam_range = linspace(log_diam_minmax(1), log_diam_minmax(2), pars.n_diams);   
sf_vals = exp(log_sf_range);
tf_vals = exp(log_tf_range);
diam_vals = exp(log_diam_range);
angle_inc = (pars.angle_range(2) - pars.angle_range(1)) / pars.nAngles;
angle_vals = pars.angle_range(1) : angle_inc : pars.angle_range(2)-angle_inc; % even distribution not including the upper limit

% Clear unused variables
% analyze_ss_responses will expect these fields to be removed if the
% respective parameter was not tested.
sfs_plus_tfs = pars.n_tfs + pars.n_sfs; % get sum before clearing
if pars.n_tfs <= 0 % if we aren't testing a range of temporal frequencies 
    pars = rmfield(pars, {'n_tfs' 'sf_fixed'});
end
if pars.n_sfs <= 0 % if we aren't testing a range of spatial frequencies 
    pars = rmfield(pars, {'n_sfs' 'tf_fixed'});
end

% Make sure specified temporal frequency is positive and doesn't exceed 0.5/IFI. 
if any(tf_vals<0) | any(tf_vals>0.5/screenstats.ifi)
    error('Stimulus temporal frequency must be between 0 and 0.5/IFI (=%gc/d) or direction may appear reversed.',0.5/screenstats.ifi)
end

if pars.SHOW_CENTER && any(diam_vals < pars.diam_inner)
    warning('At least one annulus diameter is smaller than the center-grating diameter.')
end

% Make sure grating spatial period is between 0 and 180 degrees. 
if any(1./sf_vals<0 | 1./sf_vals>180)
    error('Grating spatial period must be between 0 and 180 degrees (spatial frequency greater than 0.0056).')
end

% Convert units (degs to pix, cycles/deg to cycles/pixel, secs to frames). We are assuming that 
% aspect ratio of screen resolution = aspect ratio of physical screen size (square pixels), 
% and therefore that we can use the same the same conversion for the horz and vert axes.  
%%% Note: deg2pixels is only valid for lines bisected by the screen center. 
diam_pix_vals = deg2pixels(diam_vals,screenstats); % convert degrees to pixels
sfreq_pix_vals = 1./deg2pixels(1./sf_vals,screenstats); % convert degrees to pixels
pars.nframes = pars.stimdur/screenstats.ifi; % convert duration from seconds to flips
if isfield(pars,'sf_fixed')
    pars.sfreq_pix_fixed = 1/deg2pixels(1/pars.sf_fixed,screenstats); % convert degrees to pixels 
end
if pars.SHOW_CENTER
    pars.diam_pix_inner = deg2pixels(pars.diam_inner,screenstats); % convert degrees to pixels
end

% Warn if spatial period exceeds largest axis.   
if any(1./sfreq_pix_vals > 2*max(max(diam_pix_vals)))
    warning('Outer grating spatial period exceeds 2x largest ellipse axis.')
elseif pars.SHOW_CENTER && pars.sfreq_pix_fixed > 2*max(max(pars.diam_pix_inner))
    warning('Inner grating spatial period exceeds 2x largest ellipse axis')
end

% Create struct with element for every trial, then shuffle trials to randomize order of parameter combinations.  
ncombos = pars.repetitions*pars.nAngles*pars.n_diams*sfs_plus_tfs; % number of stimulus parameter combinations, including repetitions
par_sets = dataset({NaN(ncombos,6),'Angle','diam','sfreq','tfreq','diam_pix','sfreq_pix'},'ObsNames',cellstr(num2str((1:ncombos)'))); 

randorder = randperm(ncombos); % randomized order in which present the stimulus parameter sets
combo = 0; % counter for assigning to successive elements of the par_sets structure   
for rep = 1:pars.repetitions
    for ind_angle = 1:pars.nAngles
        for ind_diam = 1:pars.n_diams
            if isfield(pars, 'n_sfs')
                for ind_sfreq = 1:pars.n_sfs  % vary sfreq, keep tfreq fixed
                    combo = combo + 1; % assign to next element (trial) of parameters structure  
                    par_sets.Angle(randorder(combo)) = angle_vals(ind_angle);
                    par_sets.diam(randorder(combo)) = diam_vals(ind_diam);    %% for reference during later analysis            
                    par_sets.sfreq(randorder(combo)) = sf_vals(ind_sfreq);    %% for reference during later analysis
                    par_sets.tfreq(randorder(combo)) = pars.tf_fixed; % keep tf constant
                    par_sets.diam_pix(randorder(combo)) = diam_pix_vals(ind_diam);
                    par_sets.sfreq_pix(randorder(combo)) = sfreq_pix_vals(ind_sfreq); % varying sf
                end
            end
            if isfield(pars, 'n_tfs')
                for ind_tfreq = 1:pars.n_tfs; % vary tfreq, keep sfreq fixed
                    combo = combo + 1;   % assign to next element (trial) of parameters structure  
                    par_sets.Angle(randorder(combo)) = angle_vals(ind_angle);
                    par_sets.diam(randorder(combo)) = diam_vals(ind_diam);    %% for reference during later analysis              
                    par_sets.sfreq(randorder(combo)) = pars.sf_fixed;    %% for reference during later analysis   
                    par_sets.tfreq(randorder(combo)) = tf_vals(ind_tfreq); % varying tf  
                    par_sets.diam_pix(randorder(combo)) = diam_pix_vals(ind_diam);
                    par_sets.sfreq_pix(randorder(combo)) = pars.sfreq_pix_fixed; % keep sf constant
                end
            end
        end
    end
end
    
clockcheck = clock;
clear randorder rep ind_angle ind_sfreq indtfreq combo

pars.sf_vals = sf_vals;
pars.tf_vals = tf_vals;
pars.warping = 'nonwarped'; % indicate that nonwarped stimuli were used
pars.screenstats = screenstats;
pars.pulse = pulse;
pars.stim_center = stim_center;
pars.starttime = datestr(now);
pars.computer = getenv('computername');
par_sets = dataset2table(par_sets);
if ~isempty(savefile)
    save(savefile,'pars','par_sets');
end

%% Present stimuli for each combination of specified parameters.
disp(sprintf(['Now presenting surround-suppression stimuli.\n',...
    'Approximate duration = %g minutes'],ncombos*(pars.stimdur+pars.isi)/60));
if pars.send_pulse
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(pulse.stimoff_signal);    %%%%%%%%%%%%%%%%%% maybe replace with outputSingleScan and/or digital output to make faster
    daq_sess.startForeground(); % delay period at Vbase before starting iteration
end
Screen('FillRect', screenstats.winPtr, BlackIndex(screenstats.winPtr)); % set a black background before presenting stim 
Screen('Flip', screenstats.winPtr);
for trial = 1:height(par_sets)
    if commandline_output
        par_sets(trial,:)
    end
    present_stim(par_sets(trial,:), pars, stim_center, pulse, screenstats, daq_sess); % present stimuli with a different combination of the varying parameters
end

% flip to black screen
Screen('FillRect', screenstats.winPtr, BlackIndex(screenstats.winPtr)); 
Screen('Flip',screenstats.winPtr);
pause(3*pars.isi);

if exist('debugoriginal','var') %  reset PTB visual debug level if it was changed
    Screen('Preference', 'VisualDebugLevel', debugoriginal);    % set splash screen color to black
end

disp('Surround-suppression stimulus presentation complete; end recording.');

end

%% Function for presenting grating with annulus
%%%% Call this function each time parameters of the stimulus change.
% varpars = struct containing an element listing the values of varying
%       parameters for each trial
% fixpars = struct of parameters which don't change per trial ('pars' from
%       the larger function surround_suppression_stim.m)
function [] = present_stim(varpars, fixpars, stim_center, pulse, screenstats, daq_sess)

% Set parameters. First row = inner grating, second row = outer grating.
if fixpars.SHOW_CENTER
    patches_to_show = [2 1]; % show separate center and surround 
    axes_pix = [fixpars.diam_pix_inner fixpars.diam_pix_inner; varpars.diam_pix varpars.diam_pix]; % hrzt axis = vert axis

else
    patches_to_show = 2; % show center and surround as one grating
    axes_pix = [NaN NaN; varpars.diam_pix varpars.diam_pix]; % hrzt axis = vert axis
end

Angle = [varpars.Angle; varpars.Angle+fixpars.orient_offset]; % rotate outer grating orientation relative to inner grating orientation

% Get spatial and temporal parameters. 
if isfield(fixpars,'sfreq_pix_fixed') % if n_tfs > 0
    sfreq_pix(1,1) = fixpars.sfreq_pix_fixed;
else
    sfreq_pix(1,1) = NaN;
end
sfreq_pix(2,1) = varpars.sfreq_pix; % outer diam sfreq
pix_per_cycle =  1 ./ sfreq_pix; % get spatial period from spatial frequency
if isfield(fixpars,'tf_fixed') % if n_sfs > 0
    shiftperframe(1,1) = fixpars.tf_fixed * pix_per_cycle(1,1) * screenstats.ifi; % pix/frame = cycles/s * pix/cycle * s/frame
else
    shiftperframe(1,1) = NaN;
end
shiftperframe(2,1) = varpars.tfreq * pix_per_cycle(2,1) * screenstats.ifi; % pix/frame = cycles/s * pix/cycle * s/frame 

% Create the inner and outer gratings.
for patch = patches_to_show
    [x ~] = meshgrid(1:axes_pix(patch,1:2) + pix_per_cycle(patch), 1);
    grating = fixpars.brightfraction *...
        [0.5*WhiteIndex(screenstats.winPtr) + 0.5*WhiteIndex(screenstats.winPtr)*cos(2*pi*sfreq_pix(patch)*x)];
    rect_aperture(patch,:) = CenterRectOnPoint([0 0 axes_pix(patch,1:2)], stim_center(1), stim_center(2)); 
    gratingtex(patch) = Screen('MakeTexture', screenstats.winPtr, grating);
    diag = sqrt(axes_pix(patch,1)^2 + axes_pix(patch,2)^2); % make grating at least this big to cover the aperture
    rect_clear(patch,:) = CenterRectOnPoint([0 0 sqrt(2)*[diag diag]], stim_center(1), stim_center(2));
    rect_grating(patch,:) = CenterRectOnPoint([0 0 diag diag], stim_center(1), stim_center(2));
end 
 
        % Send first flip and daq pulse as close in time as possible. See
        % 'frame' loop below for comments on psychtoolbox functions. 
        % Send daq pulse right after, rather than before the first flip, or else
        % the pulse will be sent up to ifi seconds too early. 
if fixpars.send_pulse
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(pulse.stimon_signal);
end
for patch = patches_to_show
    Screen('Blendfunction', screenstats.winPtr, GL_ONE, GL_ZERO, [0 0 0 1]);   
    Screen('FillRect', screenstats.winPtr, [0 0 0 0], rect_clear(patch,:)); % maybe take smaller of winrect and this rect
    Screen('FillOval', screenstats.winPtr, [0 0 0 255], rect_aperture(patch,:));
    Screen('Blendfunction', screenstats.winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
    Screen('DrawTexture', screenstats.winPtr, gratingtex(patch), [0 0 axes_pix(patch,:)], rect_grating(patch,:), Angle(patch));
end
Screen('Flip',screenstats.winPtr);
if fixpars.send_pulse
    daq_sess.startBackground(); % following first flip, set stimchan to stimon level to indicate stimulus  
end
    %%%%% Present Stimuli %%%%%%%
    % adapted from DriftDemo5
for frame = 2:fixpars.nframes
    xoffset = mod(frame*shiftperframe, pix_per_cycle); % = how many pixels we are from the starting phase
    srcRect = [xoffset [0; 0] xoffset+axes_pix(:,1) axes_pix(:,2)]; 
    
    for patch = patches_to_show
        %%% Outer Grating   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Disable alpha-blending, restrict following drawing to alpha channel:    
        Screen('Blendfunction', screenstats.winPtr, GL_ONE, GL_ZERO, [0 0 0 1]);   

        % Clear 'dstRect' region of framebuffers alpha channel to zero:
        Screen('FillRect', screenstats.winPtr, [0 0 0 0], rect_clear(patch,:)); % maybe take smaller of winrect and this rect

        % Fill oval-shaped 'dstRect' region with an alpha value of 255:
        Screen('FillOval', screenstats.winPtr, [0 0 0 255], rect_aperture(patch,:));

        % Enable DeSTination alpha blending and reenable drawing to all
        % color channels. Following drawing commands will only draw there
        % the alpha value if the framebuffer is greater than zero, ie., in
        % our case, inside the circular 'rect_aperture' area where alpha has
        % been set to 255 by our 'FillOval' command:
        Screen('Blendfunction', screenstats.winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);
        
        Screen('DrawTexture', screenstats.winPtr, gratingtex(patch), srcRect(patch,:), rect_grating(patch,:), Angle(patch));
    end
    vbl = Screen('Flip',screenstats.winPtr); % could make vbl a vector, compare predicted to actual flip times
end

% Flip back to black screen when stimulus presentation is done. 
Screen('FillRect', screenstats.winPtr, BlackIndex(screenstats.winPtr)); 
if fixpars.send_pulse
    stop(daq_sess); % end stimon pulse if it is still running
    daq_sess.queueOutputData(pulse.stimoff_signal);
end
Screen('Flip',screenstats.winPtr);
if fixpars.send_pulse
    daq_sess.startBackground(); % set stimchan to stimoff level to indicate stimulus offset    
end
pause(fixpars.isi); %% wait isi seconds before presenting next stim... should use tic toc instead to be more accurate
end



%% deg2pixels: function to convert degrees subtended to pixels on the screen
% All measurements and output are taken in one dimension of interest.
% Note: deg2pixels is only valid for lengths bisected by the screen center.  
% We are assuming that the stimulus is symmetric and centered on the screen center.
    % deg is desired length of stimulus on screen in degrees (can be a vector)   
function pixels = deg2pixels(deg,screenstats)
    eye2screen_center_inches = sqrt( screenstats.eye2screen_top_bottom_inches^2 - (screenstats.screen_height_inches/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
    stim_diam = 2*eye2screen_center_inches*tan(deg2rad(deg)/2); % because (stim_diam/2) / eye2screen_center = tan(deg/2)
    pixels = stim_diam * screenstats.height/screenstats.screen_height_inches;  %% convert length to pixels
end

%% pixels2deg: inverse of deg2pixels, with the same caveats
function deg = pixels2deg(pixels,screenstats)
   eye2screen_center = sqrt( screenstats.eye2screen_top_bottom_inches^2 - (screenstats.screen_height_inches/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
   stim_diam = pixels * screenstats.screen_height_inches/screenstats.height; % convert pixels to length
   deg = rad2deg( 2*atan(stim_diam /(2*eye2screen_center))); % because stim_diam/2 = eye2screen_center*tan(deg/2)
end


%%%%% maybe map 'average' preferred sf and tf of all units rather than
%%%%% starting with fixed value; Vaiceliunaite stared with fixed value   
%%%%%%% possibility: analyze mean tf and sf tuning and hwhm; compare to
%%%%%%% wq's and my collected data to determine likelihood that shank is in
%%%%%%% patch/nonpatch; if not in desired location, penetrate in new area
%%%%%%% (far enough away from first penetratation)