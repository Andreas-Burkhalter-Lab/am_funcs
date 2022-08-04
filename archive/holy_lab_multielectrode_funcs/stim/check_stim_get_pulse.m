%CHECK_STIM_GET_PULSE Check stimulus paramaters and get stim channel pulse
%signals
%%% singleVStimOn == 0 if the VStimOn value is being used to convery
%%% information about the stimulus other than timing (rf mapping only),
%%% otherwise singleVStimOn == 1
% To be called by rf mapping, stimpref tuning, and ssmod scripts
%%% last updated 1/28/16
function pulsepars_out = check_stim_get_pulse(stimpars,pulsepars_in,daq_sess,singleVStimOn,commonvars)

%% Stim channel signal setup
pulsepars_out = pulsepars_in;
if singleVStimOn % everything except rf mapping
    stimon_dur_scans = stimpars.dur_stim * daq_sess.Rate;    %% get stim-on duration in scans for the stim channel pulsepars
    stimoff_dur_scans = max([stimpars.isi * daq_sess.Rate, 1]); %% get stim-off duration in scans for the stim channel pulsepars
    pulsepars_out.stimon_signal = pulsepars_in.Vstimon * ones(stimon_dur_scans,1); % signal to send at stim onset
    pulsepars_out.stimoff_signal = pulsepars_in.Vstimoff * ones(stimoff_dur_scans,1); % signal to send at stim offset
    pulsepars_out = rmfield(pulsepars_out,'V1_dur');
elseif ~singleVStimOn
    pulsepars_out = rmfield(pulsepars_in,'Vstimon'); % this vstimon value is not used in rf mapping
end
    
%% Parameter checks
% Make sure that we aren't trying to hold stim chan at Vbase longer than
% the intertrial interval.
if pulsepars_out.Vstimoff_dur_post > stimpars.isi 
    error('Stim channel post-stimulus duration (%gs) exceeds',...
        'stimulus intertrial interval (%gs)',Vbase_dur_post,pars.isi)
end

% Make sure grating spatial period is between 0 and 180 degrees. 
if (isfield(stimpars,'sf') && [1/stimpars.sf<0 || 1/stimpars.sf>180]) ||... % 1 sf
     (isfield(stimpars,'sf_vals') && [1/max(stimpars.sf_vals)<0 | 1/min(stimpars.sf_vals)>180]) % multiple sfs
    error('Grating spatial period must be between 0 and 180 degrees (spatial frequency greater than 0.0056).')
end

% Make sure that spatial period is at least 2 pixels large (should be much larger).  
if (isfield(stimpars,'sf') && [deg2pixels(1/stimpars.sf,commonvars) < ...
        2*commonvars.screen_height/commonvars.screenstats.height]) ||... % 1 sf
      (isfield(stimpars,'sf_vals') && [min(deg2pixels(1./stimpars.sf_vals,commonvars)) < ...
      2*commonvars.screen_height/commonvars.screenstats.height]) % multiple sfs
    error('Grating spatial period must be at least 2 pixels wide.')
end

% Make sure specified temporal frequency is positive and doesn't exceed 0.5/IFI. 
if isfield(stimpars,'tf') && [stimpars.tf<0 || stimpars.tf>0.5/commonvars.ifi] ||... % 1 tf
      (isfield(stimpars,'tf_vals') &&...% multiple tfs
      ( min(stimpars.tf_vals) <=0  || max(stimpars.tf_vals) > 0.5/commonvars.ifi )) 
    error('Temporal frequencies must be between 0 and 0.5/IFI (=%g) or direction may appear reversed.',0.5/ifi)
end

% Warn if spatial period exceeds stimulus diameter.
%%% should not be an issue if amp is reasonably small.
if any(strcmp(stimpars.tuningParameter,{'diam','ssmod'})) % multiple diams
    mindiam = min(stimpars.diam_vals);
    maxdiam = max(stimpars.diam_vals); % for later check
else % 1 diam
    mindiam = stimpars.diam;
    maxdiam = stimpars.diam; % for later check
end
if any(strcmp(stimpars.tuningParameter,{'sf','ssmod'})) % multiple sfs
    minsf = min(stimpars.sf_vals);
    max_spatialperiod = 1/minsf;
else % 1 sf
    max_spatialperiod = 1/stimpars.sf;
end
if max_spatialperiod > 2*mindiam
    warning('Grating spatial period exceeds 2x stimulus diameter in RF center location mapping.')
end

% Check that all stimuli in rf grid fit on screen and shrink if necessary. 
if isfield(stimpars,'rf_center') % if rf mapping has already been done
    stimcheck.diamTried_pix = deg2pixels(maxdiam,commonvars); % max diam in pixels
    stimcheck.edge_left = stimpars.rf_center(1) - stimcheck.diamTried_pix/2; % x coordinate of the left edge of the largest grating
    stimcheck.edge_right = stimpars.rf_center(1) + stimcheck.diamTried_pix/2;
    stimcheck.edge_top = stimpars.rf_center(2) - stimcheck.diamTried_pix/2; % y coordinate of the top edge of the largest grating
    stimcheck.edge_bottom = stimpars.rf_center(2) + stimcheck.diamTried_pix/2;
    stimcheck.pix_cut_off = max([... % positive values indicate that part of at least one grating is off of the screen
        -stimcheck.edge_left + commonvars.winrect(1),...
        stimcheck.edge_right - commonvars.winrect(3),...
        -stimcheck.edge_top + commonvars.winrect(2),...
        stimcheck.edge_bottom - commonvars.winrect(4)]);
    if stimcheck.pix_cut_off > 0 % if we need to shrink the gratings to fit on the stimulus screen
        stimcheck.diamMaxAllowable_pix = floor(stimcheck.diamTried_pix - 2*stimcheck.pix_cut_off); % max diam that will fit on screen
        stimcheck.diamMaxAllowable = pixels2deg(stimcheck.diamMaxAllowable_pix); % largest diam that will fit in degrees
        if min(pars.diam_minmax) >= stimcheck.diamMaxAllowable % if smallest grating is cut off by screen
            error(['Smallest grating diameter (%g pixels, %g deg) is cut off by screen.',...
                '\nGratings must have diameter %g pixels (%g deg) or less to fit on screen.'],...
                deg2pixels(min(pars.diam_minmax)), min(pars.diam_minmax),...
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