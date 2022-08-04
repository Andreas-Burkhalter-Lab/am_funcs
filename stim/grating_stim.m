function [] = grating_stim(stimpars, pulse, savefile, screenstats)
%%% arg1 = struct of stimulus parameters ('stimpars' from run_experiment.m)
%%% arg2 = parameters of synchronization pulse
%%% arg3 = name of file to save stim data into
%%% arg4 = info about screen setup
%%%%% last edited 18/12/04 on thermaltake












%% to add: change stimpars struct fields h spacing, v spacing, grid center, diam if these values are changed during setup





%% Setup
if ~stimpars.save_stim_file
    checksave = input('Not saving stim file. Press Enter to continue.','s');
    savefile = '';
end    

commandline_output = 0; % show params of each trial as it is presented

if exist(savefile,'file') | exist(strcat(savefile,'.mat'),'file')
    overwrite = input(sprintf('File ''%s'' already exists. Enter ''y'' to overwrite.\n', savefile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
end

debugoriginal = Screen('Preference','VisualDebugLevel');
if stimpars.show_splash_screen
    Screen('Preference', 'VisualDebugLevel', 4); % more thorough but show white splash screen  
elseif ~stimpars.show_splash_screen
    Screen('Preference','VisualDebugLevel',0);   % turn off white splash screen and error warning screen
end
    
if ~stimpars.do_sync_tests
    Screen('Preference', 'SkipSyncTests', 1);
else
    Screen('Preference','SkypSyncTests',0);
end

stimpars.computer = getenv('computername');
if isempty(Screen('Windows'))

    if strcmp(stimpars.computer, 'DESKTOP-H739EDS') % if on AB lab thermaltake
    %     myScreen = 1;
          myScreen = 2;
    else % on any other computer
    %     myScreen = 0;
        myScreen = max(Screen('Screens'));
    end
        [ screenstats.winPtr , winrect ] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
    else
        screenstats.winPtr = Screen('Windows');
        screenstats.winPtr = screenstats.winPtr(1);
end
screenstats = copyAllFields(screenstats, Screen('Resolution',screenstats.winPtr));
screenstats.ifi = 1/screenstats.hz;

% % % % % % % % % Create stim chan signals.
if stimpars.send_pulse
    if stimpars.isi <= 0 
        error('ISI must be greater than zero to send stim pulse.')
    end
    daq_sess = daq.createSession('ni');
    try 
        chan = addAnalogOutputChannel(daq_sess,'Dev2','ao0', 'Voltage'); % if the nidaq card is Dev2
    catch
        chan = addAnalogOutputChannel(daq_sess,'Dev1','ao0', 'Voltage');
    end
    pulse.stimon_dur_scans = stimpars.stimdur * daq_sess.Rate;    %% get stim-on duration in scans for the stim channel pulse
    pulse.stimoff_dur_scans = stimpars.isi * daq_sess.Rate; %% get stim-off duration in scans for the stim channel pulse
    pulse.stimon_signal = pulse.Vstimon * ones(pulse.stimon_dur_scans,1); % signal to send at stim onset
    pulse.stimoff_signal = pulse.Vstimoff * ones(pulse.stimoff_dur_scans,1); % signal to send at stim offsetk to Vbase
else
    daq_sess = struct;
end
    
% Check that all stimuli fit on screen and shrink if necessary. 
if stimpars.do_stimcheck
    stimcheck.diamTried_pix = deg2pixels(max(stimpars.diam_minmax),screenstats); % max diam specified in run_experiment in pixels
    stimcheck.edge_left = stimpars.stim_center(1) - stimcheck.diamTried_pix/2; % x coordinate of the left edge of the largest grating
    stimcheck.edge_right = stimpars.stim_center(1) + stimcheck.diamTried_pix/2;
    stimcheck.edge_top = stimpars.stim_center(2) - stimcheck.diamTried_pix/2; % y coordinate of the top edge of the largest grating
    stimcheck.edge_bottom = stimpars.stim_center(2) + stimcheck.diamTried_pix/2;
    stimcheck.pix_cut_off = max([... % positive values indicate that part of at least one grating is off of the screen
        -stimcheck.edge_left + winrect(1),...
        stimcheck.edge_right - winrect(3),...
        -stimcheck.edge_top + winrect(2),...
        stimcheck.edge_bottom - winrect(4)]);
    if stimcheck.pix_cut_off > 0 % if we need to shrink the gratings to fit on the stimulus screen
        stimcheck.diamMaxAllowable_pix = floor(stimcheck.diamTried_pix - 2*stimcheck.pix_cut_off); % max diam that will fit on screen
        stimcheck.diamMaxAllowable = pixels2deg(stimcheck.diamMaxAllowable_pix,screenstats); % largest diam that will fit in degrees
        if min(stimpars.diam_minmax) >= stimcheck.diamMaxAllowable % if smallest grating is cut off by screen
            error(['Smallest grating diameter (%g pixels, %g deg) is cut off by screen.',...
                '\nGratings must have diameter %g pixels (%g deg) or less to fit on screen.'],...
                deg2pixels(min(stimpars.diam_minmax),screenstats), min(stimpars.diam_minmax),...
                stimcheck.diamMaxAllowable_pix, stimcheck.diamMaxAllowable)
        end
        stimcheck.shrink_ok = input(sprintf(['Largest gratings are cut off by screen. ',...
            'Enter ''y'' to reduce maximum grating diameter from',...
            '\n      %g pixels (%g deg) to %g pixels (%g deg).\n'],...
            stimcheck.diamTried_pix, max(stimpars.diam_minmax), stimcheck.diamMaxAllowable_pix, stimcheck.diamMaxAllowable),'s');
        if strcmp(stimcheck.shrink_ok,'y')
            stimpars.diam_minmax(2) = stimcheck.diamMaxAllowable;
        else
            error('Will not shrink cut-off gratings.')
        end
    end
end
    
% Create the distributions of annulus sf, tf, sizes, and orients.
% Sf, tf, and size are on log-scales including highest value and orient is on
% a linear scale not including highest value. 
%%% maybe Vaiceliunaite-type distribution of diams is better than exponential; seems
%%% to capture steeper right-side slope which Gao doesn't see
log_sf_minmax = log(stimpars.sf_minmax);
log_tf_minmax = log(stimpars.tf_minmax);
log_diam_minmax = log(stimpars.diam_minmax);
log_sf_range = linspace(log_sf_minmax(1),log_sf_minmax(2),stimpars.n_sfs);
log_tf_range = linspace(log_tf_minmax(1),log_tf_minmax(2),stimpars.n_tfs);
log_diam_range = linspace(log_diam_minmax(1), log_diam_minmax(2), stimpars.n_diams);   
sf_vals = exp(log_sf_range);
tf_vals = exp(log_tf_range);
diam_vals = exp(log_diam_range);
orient_inc = (stimpars.orient_range(2) - stimpars.orient_range(1)) / stimpars.n_orients;
orient_vals = stimpars.orient_range(1) : orient_inc : stimpars.orient_range(2)-orient_inc; % even distribution not including the upper limit

% Clear unused variables
sfs_plus_tfs_plus_orients = stimpars.n_tfs + stimpars.n_sfs + stimpars.n_orients; % get sum before clearing
if stimpars.n_tfs <= 0 % if we aren't testing a range of temporal frequencies 
    stimpars = rmfield(stimpars, {'n_tfs'});
end
if stimpars.n_sfs <= 0 % if we aren't testing a range of spatial frequencies 
    stimpars = rmfield(stimpars, {'n_sfs'});
end
if stimpars.n_orients <= 0 % if we aren't testing a range of orientations 
    stimpars = rmfield(stimpars, {'n_orients'});
end

% Make sure specified temporal frequency is positive and doesn't exceed 0.5/IFI. 
if any(tf_vals<0) | any(tf_vals>0.5/screenstats.ifi)
    error('Stimulus temporal frequency must be between 0 and 0.5/IFI (=%gc/d) or direction may appear reversed.',0.5/screenstats.ifi)
end

if stimpars.SHOW_CENTER && any(diam_vals < stimpars.diam_inner)
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
stimpars.nframes = stimpars.stimdur/screenstats.ifi; % convert duration from seconds to flips
if isfield(stimpars,'sf_fixed')
    stimpars.sfreq_pix_fixed = 1/deg2pixels(1/stimpars.sf_fixed,screenstats); % convert degrees to pixels 
end
if stimpars.SHOW_CENTER
    stimpars.diam_pix_inner = deg2pixels(stimpars.diam_inner,screenstats); % convert degrees to pixels
end

% Warn if spatial period exceeds largest axis.   
if any(1./sfreq_pix_vals > 2*max(max(diam_pix_vals)))
    warning('Outer grating spatial period exceeds 2x largest ellipse axis.')
elseif stimpars.SHOW_CENTER && stimpars.sfreq_pix_fixed > 2*max(max(stimpars.diam_pix_inner))
    warning('Inner grating spatial period exceeds 2x largest ellipse axis')
end

% Create struct with element for every trial, then shuffle trials to randomize order of parameter combinations.  
ncombos = stimpars.repetitions * [stimpars.n_diams*sfs_plus_tfs_plus_orients  + double(stimpars.include_blank_trials)]; % number of stimulus parameter combinations plus blank trials, including repetitions
stimpar_sets = dataset({NaN(ncombos,6),'orient','sfreq','tfreq','diam','sfreq_pix','diam_pix'}); 
stimpar_sets.isblank = false(ncombos,1);

randorder = randperm(ncombos); % randomized order in which present the stimulus parameter sets
combo = 0; % counter for assigning to successive elements of the par_sets structure   
for rep = 1:stimpars.repetitions
    for ind_diam = 1:stimpars.n_diams
        if isfield(stimpars, 'n_sfs')
            for ind_sfreq = 1:stimpars.n_sfs  % vary sfreq
                combo = combo + 1; % assign to next element (trial) of parameters structure  
                stimpar_sets.orient(randorder(combo)) = stimpars.orient_fixed;
                stimpar_sets.diam(randorder(combo)) = diam_vals(ind_diam);    %% for reference during later analysis            
                stimpar_sets.sfreq(randorder(combo)) = sf_vals(ind_sfreq);    %% for reference during later analysis
                stimpar_sets.tfreq(randorder(combo)) = stimpars.tf_fixed; % keep tf constant
                stimpar_sets.diam_pix(randorder(combo)) = diam_pix_vals(ind_diam);
                stimpar_sets.sfreq_pix(randorder(combo)) = sfreq_pix_vals(ind_sfreq); % varying sf
            end
        end
        if isfield(stimpars, 'n_tfs')
            for ind_tfreq = 1:stimpars.n_tfs % vary tfreq
                combo = combo + 1;   % assign to next element (trial) of parameters structure  
                stimpar_sets.orient(randorder(combo)) = stimpars.orient_fixed;
                stimpar_sets.diam(randorder(combo)) = diam_vals(ind_diam);    %% for reference during later analysis              
                stimpar_sets.sfreq(randorder(combo)) = stimpars.sf_fixed;    %% for reference during later analysis   
                stimpar_sets.tfreq(randorder(combo)) = tf_vals(ind_tfreq); % varying tf  
                stimpar_sets.diam_pix(randorder(combo)) = diam_pix_vals(ind_diam);
                stimpar_sets.sfreq_pix(randorder(combo)) = stimpars.sfreq_pix_fixed; % keep sf constant
            end
        end
        if isfield(stimpars,'n_orients')
            for ind_orient = 1:stimpars.n_orients
                combo = combo + 1;
                stimpar_sets.orient(randorder(combo)) = orient_vals(ind_orient);
                stimpar_sets.diam(randorder(combo)) = diam_vals(ind_diam);
                stimpar_sets.sfreq(randorder(combo)) = stimpars.sf_fixed;
                stimpar_sets.tfreq(randorder(combo)) = stimpars.tf_fixed;
                stimpar_sets.diam_pix(randorder(combo)) = diam_pix_vals(ind_diam);
                stimpar_sets.sfreq_pix(randorder(combo)) = stimpars.sfreq_pix_fixed;
            end
        end
    end
end
%%% label blanks if present    
if stimpars.include_blank_trials
    for rep = 1:stimpars.repetitions
    combo = combo+1;
    stimpar_sets.isblank(randorder(combo)) = true;
    stimpar_sets.orient(randorder(combo)) = NaN;
    stimpar_sets.diam(randorder(combo)) = NaN;
    stimpar_sets.sfreq(randorder(combo)) = NaN;
    stimpar_sets.tfreq(randorder(combo)) = NaN;
    stimpar_sets.diam_pix(randorder(combo)) = NaN;
    stimpar_sets.sfreq_pix(randorder(combo)) = NaN;
    end
end

clockcheck = clock;
clear randorder rep ind_orient ind_sfreq indtfreq combo

stimpars.sf_vals = sf_vals;
stimpars.tf_vals = tf_vals;
stimpars.orient_vals = orient_vals;
stimpars.diam_vals = diam_vals;
stimpars.warping = 'nonwarped'; % indicate that nonwarped stimuli were used
stimpars.experiment_type = 'sf_tf_orient_diam'; 
stimpars.screenstats = screenstats;
stimpars.pulse = pulse;
stimpars.starttime = datestr(now);
stimpar_sets = dataset2table(stimpar_sets);
if stimpars.save_stim_file
    save(savefile,'stimpars','stimpar_sets');
end

%% Present stimuli for each combination of specified parameters.
disp(sprintf(['Now presenting grating stimuli.\n',...
    'Approximate duration = %g minutes'],ncombos*(stimpars.stimdur+stimpars.isi)/60));
if stimpars.send_pulse
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(pulse.stimoff_signal);    %%%%%%%%%%%%%%%%%% maybe replace with outputSingleScan and/or digital output to make faster
    daq_sess.startForeground(); % delay period at Vbase before starting iteration
end
Screen('FillRect', screenstats.winPtr, BlackIndex(screenstats.winPtr)); % set a black background before presenting stim 
Screen('Flip', screenstats.winPtr);
ntrials = height(stimpar_sets);
for trial = 1:ntrials
    if commandline_output
        [table(sprintf('%g/%g',trial,ntrials),'VariableNames',{'trial'}), stimpar_sets(trial,:)]
    end
    present_stim(stimpar_sets(trial,:), stimpars, stimpars.stim_center, pulse, screenstats, daq_sess); % present stimuli with a different combination of the varying parameters
end

% flip to black screen
Screen('FillRect', screenstats.winPtr, BlackIndex(screenstats.winPtr)); 
Screen('Flip',screenstats.winPtr);
pause(3*stimpars.isi); % short pause before starting stim

if exist('debugoriginal','var') %  reset PTB visual debug level if it was changed
    Screen('Preference', 'VisualDebugLevel', debugoriginal);    % set splash screen color to black
end

disp(['Stimulus presentation complete; stim data save into ' savefile '.']);
% Screen('CloseAll')
commandwindow;

end



%% Function for presenting grating with annulus
%%%% Call this function each time parameters of the stimulus change.
% varpars = struct containing an element listing the values of varying
%       parameters for each trial
% fixpars = struct of parameters which don't change per trial ('pars' from
%       the larger function surround_suppression_stim.m)
function [] = present_stim(varpars, fixpars, stim_center, pulse, screenstats, daq_sess)

if ~varpars.isblank % if not a blank trial
    % Set parameters. First row = inner grating, second row = outer grating.
    if fixpars.SHOW_CENTER
        patches_to_show = [2 1]; % show separate center and surround 
        axes_pix = [fixpars.diam_pix_inner fixpars.diam_pix_inner; varpars.diam_pix varpars.diam_pix]; % hrzt axis = vert axis

    else
        patches_to_show = 2; % show center and surround as one grating
        axes_pix = [NaN NaN; varpars.diam_pix varpars.diam_pix]; % hrzt axis = vert axis
    end

    this_orient = [varpars.orient; varpars.orient+fixpars.orient_offset]; % rotate outer grating orientation relative to inner grating orientation

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
        Screen('DrawTexture', screenstats.winPtr, gratingtex(patch), [0 0 axes_pix(patch,:)], rect_grating(patch,:), this_orient(patch));
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

            Screen('DrawTexture', screenstats.winPtr, gratingtex(patch), srcRect(patch,:), rect_grating(patch,:), this_orient(patch));
        end
        vbl = Screen('Flip',screenstats.winPtr); % could make vbl a vector, compare predicted to actual flip times
    end
    %
elseif varpars.isblank % if blank trial
    if fixpars.send_pulse
        stop(daq_sess); %% end any signals that are already being sent to the stim channel
        daq_sess.queueOutputData(pulse.stimon_signal);
        daq_sess.startBackground(); % indicate 'start' of blank trial
    end
    pause(fixpars.stimdur); % pause instead of showing stimulus
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
pause(fixpars.isi);
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
