%SSMODSTIM_PREPARESTIM: Create and/or load stimulus gratings for ssstim_main.
% Access the prerendered gratings, load variables, and put prerendered parameters into 'ssmod_stim'.
% Note: with large gratings files, it appears faster to load data (aside
% from gratings) rather than accessing it via matfile. 
%%% last updated 1/27/16 on msi

fprintf(['Loading prerendered parameters and memory-mapping ''' prerender_file '''...'])
matfileTic = tic; % for timing how long matfile() takes
load(prerender_file,'perp2meshx','perp2meshy',...
    'xstraight','ystraight','eye2screen_center','eye2screen_center_pix','scrnperpx','scrnperpy',...
    'sf_vals','tf_vals','orient_vals','prerend','outGratTrialTable');
matobj = matfile(prerender_file);
fprintf([' took ' num2str(toc(matfileTic)) 's.\n']);

% Copy parameters from the prerended gratings file; copy without overwrite
% fields from defaultssmod_stim.
prerend_fields = fieldnames(prerend);
for f = 1:length(prerend_fields) % analyze_ss_responses looks for a struct named 'ssmod_stim'
    ssmod_stim = setfield(ssmod_stim, prerend_fields{f}, getfield(prerend,prerend_fields{f}));
end
defaultparsfields = fieldnames(defaultpars);
for f = 1:length(defaultparsfields)
    if ~isfield(ssmod_stim,defaultparsfields{f})
        ssmod_stim = setfield(ssmod_stim,defaultparsfields{f},getfield(defaultpars,defaultparsfields{f}));
    end
end

% Check that prerendered diam and stim centers match those from run_experiment.   
if ssmod_stim.prerenderedLocAndDiam && (any(ssmod_stim.rf_center ~= stimpref.rf_center))
    error('Specified stim center must match stimpref.rf_center of prerendered gratings.')
elseif ssmod_stim.prerenderedLocAndDiam && (ssmod_stim.max_diam ~= ssmod_stim.diam_minmax(2)) 
    error('Specified largest diameter must match largest diameter of prerendered gratings.')
end

% Check that all stimuli fit on screen and shrink if necessary. 
% edgeVals are values from the specified screen edge listing the aperture
%     radius (in radians) required from stimpref.rf_center to contain each pixel.
perp2stimcenterx = scrnperpx - stimpref.rf_center(1); % x pixels between screen center and stim center
perp2stimcentery = scrnperpy - stimpref.rf_center(2); % y pixels between screen center and stim center
edgeVals.left = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(:,1) + perp2stimcentery*perp2meshy(:,1)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(:,1).^2+perp2meshy(:,1).^2)) );
edgeVals.right = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(:,end) + perp2stimcentery*perp2meshy(:,end)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(:,end).^2+perp2meshy(:,end).^2)) );
edgeVals.top = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(1,:) + perp2stimcentery*perp2meshy(1,:)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(1,:).^2+perp2meshy(1,:).^2)) );
edgeVals.bottom = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(end,:) + perp2stimcentery*perp2meshy(end,:)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(end,:).^2+perp2meshy(end,:).^2)) );
radMaxAllowable = min([edgeVals.left; edgeVals.right; edgeVals.top'; edgeVals.bottom']); % largest radius that can fit on screen

if ssmod_stim.SHOW_CENTER && deg2rad(ssmod_stim.diam_inner/2) > radMaxAllowable % if inner diam does not fit on screen
    error(sprintf(['Inner grating does not fit on screen. Diameters must be \n    ',...
        num2str(2*rad2deg(radMaxAllowable)) ' degrees or smaller to fit on screen.'])) % convert to diameter radians
end
if deg2rad(max(ssmod_stim.diam_minmax)/2) > radMaxAllowable
    commandwindow;
    shrink_ok = input(sprintf(['Largest gratings (%g deg) are cut off by screen. ',... % give option to automatically shrink
        'Enter ''y'' to reduce maximum grating diameter',...
        '\n      from %g degrees to %g degrees.\n'],...
        max(ssmod_stim.diam_minmax), max(ssmod_stim.diam_minmax), 2*rad2deg(radMaxAllowable)),'s');
    if strcmp(shrink_ok,'y')
        ssmod_stim.diam_minmax(2) = 2*rad2deg(radMaxAllowable);
    else
        error('Will not shrink cut-off gratings.')
    end
end

% Create the distributions of outer-grating annulus sizes.
%%% maybe Vaiceliunaite-type distribution of diams is better than exponential; seems
%%% to capture steeper right-side slope which Gao doesn't see
log_diam_minmax = log(ssmod_stim.diam_minmax);
log_diam_range = linspace(log_diam_minmax(1), log_diam_minmax(2), ssmod_stim.n_diams);   
diam_vals = exp(log_diam_range);

if commonvars.send_pulse
    % Stim parameter checks, get pulse signals.
    ssmod_stim_tocheck = ssmod_stim;
    ssmod_stim_tocheck.diam_vals = diam_vals;
    ssmod_stim_tocheck.sf_vals = sf_vals;
    ssmod_stim_tocheck.tf_vals = tf_vals;
    ssmod_stim_tocheck.orient_vals = orient_vals;
    singleVStimOn = 1;  % only one vstimon value will be used
    pulsepars = check_stim_get_pulse(ssmod_stim_tocheck,pulsepars,daq_sess,singleVStimOn,commonvars);
end
    
% Table to contain stimulus parameters. Each row is a trial.
% trialTime variables are for synchronizing with mouse tracking timing. 
% Vivid can read dataset objects but not table objects. 
% Create struct with element for every trial, then shuffle trials to randomize order of parameter combinations. 
sfs_plus_tfs = ssmod_stim.n_tfs + ssmod_stim.n_sfs; % get sum before clearing
n_unique_parsets = ssmod_stim.n_orients*ssmod_stim.n_diams*sfs_plus_tfs; % number of stimulus parameter combinations, including repetitions
pars_for_table = repmat({NaN(n_unique_parsets,1)},1,4);
stimrec_ssmod = table(pars_for_table{:},'VariableNames',{'orient','diam','sf','tf'}); 

% Clear unused variables
% analyze_ss_responses will expect these fields to be removed if the
% respective parameter was not tested.
if ssmod_stim.n_tfs <= 0 % if we aren't testing a range of temporal frequencies 
    ssmod_stim = rmfield(ssmod_stim, {'n_tfs' 'sf_fixed'});
end
if ssmod_stim.n_sfs <= 0 % if we aren't testing a range of spatial frequencies 
    ssmod_stim = rmfield(ssmod_stim, {'n_sfs' 'tf_fixed'});
end

combo = 0; % counter for assigning to successive elements of the stimrec_ssmod structure
for ind_angle = 1:ssmod_stim.n_orients
    for ind_diam = 1:ssmod_stim.n_diams
        if isfield(ssmod_stim, 'n_sfs')
            for ind_sfreq = 1:ssmod_stim.n_sfs  % vary sfreq, keep tfreq fixed
                combo = combo + 1; % assign to next element (trial) of parameters structure  
                stimrec_ssmod.orient(combo) = orient_vals(ind_angle);
                stimrec_ssmod.diam(combo) = diam_vals(ind_diam);    %% for reference during later analysis            
                stimrec_ssmod.sf(combo) = sf_vals(ind_sfreq);    %% for reference during later analysis
                stimrec_ssmod.tf(combo) = ssmod_stim.tf_fixed; % keep tf constant
            end
        end
        if isfield(ssmod_stim, 'n_tfs')
            for ind_tfreq = 1:ssmod_stim.n_tfs % vary tfreq, keep sfreq fixed
                combo = combo + 1;   % assign to next element (trial) of parameters structure  
                stimrec_ssmod.orient(combo) = orient_vals(ind_angle);
                stimrec_ssmod.diam(combo) = diam_vals(ind_diam);    %% for reference during later analysis              
                stimrec_ssmod.sf(combo) = ssmod_stim.sf_fixed;    %% for reference during later analysis   
                stimrec_ssmod.tf(combo) = tf_vals(ind_tfreq); % varying tf  
            end
        end
    end
end
    
grating_table_rows_table = table(NaN(n_unique_parsets,size(outGratTrialTable.grating_table_rows,2)),'VariableNames',{'grating_table_rows'});
outer_apertex_table = table(NaN(n_unique_parsets,1),'VariableNames',{'outer_apertex'}); % outer aperture texture pointers
outer_aperrect_table = table(cell(n_unique_parsets,1),'VariableNames',{'outer_apt_rect'}); % outer aperture rect values
outerGratSrcRect_table = table(NaN(n_unique_parsets,4),'VariableNames',{'outerGratSrcRect'}); % outer grating srcRect values
stimrec_ssmod = [stimrec_ssmod outer_apertex_table outer_aperrect_table outerGratSrcRect_table grating_table_rows_table];
stimrec_ssmod.flipsPerCycle = 1 ./ (commonvars.ifi .* stimrec_ssmod.tf);

% Get outer grating frame indices for each trial
for trial = 1:n_unique_parsets
    match = stimrec_ssmod.orient(trial) == outGratTrialTable.orient &...
            stimrec_ssmod.sf(trial) == outGratTrialTable.sf &...
            stimrec_ssmod.tf(trial) == outGratTrialTable.tf;
    stimrec_ssmod.grating_table_rows(trial,:) = outGratTrialTable.grating_table_rows(match,:);
end

%%% need to edit these back in to make trackmouse functional
% mouseScanStimOn_table = table(NaN(size(diamtable)),'VariableNames',{'mouseScanStimOn'}); % mouse scan at stim start per trial
% mouseScanStimOff_table = table(NaN(size(diamtable)),'VariableNames',{'mouseScanStimOff'}); % mouse scan at stim end per trial

%% Inner Gratings
if ssmod_stim.SHOW_CENTER && ~ssmod_stim.prerenderedLocAndDiam % should probably always prerender inner gratings rather than use this method
    ssmod_stim.diam_inner = input(['Inner gratings not prerendendered in ' prerender_file...
        '\n Type in inner grating diameter (in degrees) then press Enter.']); % set diam_inner here, rather than as a param in run_experiment
    ssmod_stim.amp_inner = input('Type in inner grating amp (max=1) then press Enter.'); % set amp_inner here, rather than as a param in run_experiment
    fprintf('Constructing apertures and inner gratings...')
    apertureTic = tic; 
    
    % Construct aperture for the inner grating. 
    apt_fullscreen = deg2rad(ssmod_stim.diam_inner/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    inner_apt_rect = [min(aptx) min(apty) max(aptx) max(apty)];
    inner_aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of aperture into alpha; RGB=0
    inner_apertex = Screen('MakeTexture', winPtr, inner_aperture);

    % Construct inner gratings of different angles  for each frame. Use the
    % angles, sf_fixed, and tf_fixed already specified by matobj. 
    % (Could prerender these inner gratings along with outer gratings.)
    nInnerFrames = ssmod_stim.nframes; % have not yet implemented reusing frames for this option
    inGratPtrTable  = table(repmat({NaN(1,ssmod_stim.nframes)},n_unique_parsets,1),'VariableNames',{'inner_gratingtex'}); % inner grating texture pointers
    stimrec_ssmod = [stimrec_ssmod inGratPtrTable]; clear inGratPtrTable
    for Angle = 1:length(orient_vals); % Make the inner gratings. Gratings move along the yrotate dimension.
        spatial_period = 1/ssmod_stim.sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*ssmod_stim.tf_fixed;  % convert from cycles per second to radians per second
        thetarad = deg2rad(orient_vals(1,Angle));
        xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
        yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
        eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % angles in radians; pi/2 - acos is faster than asin
        thisAngle_ptrs = inf(1,ssmod_stim.nframes);
        for frame = 1:ssmod_stim.nframes 
            grating = 0.5*WhiteIndex(winPtr) + 0.5*WhiteIndex(winPtr)*ssmod_stim.amp_inner*...
                cos(2*pi*sfrad*(eye2screen) - (frame-1)*ifi*tfrad); % from Marshel et al. 2011
            thisAngle_ptrs(frame) = Screen('MakeTexture', winPtr, grating);
        end
        match = stimrec_ssmod.orient == orient_vals(1,Angle);
        stimrec_ssmod.inner_gratingtex(match) = {thisAngle_ptrs};
    end
elseif ssmod_stim.SHOW_CENTER && ssmod_stim.prerenderedLocAndDiam  % if using prerendered inner gratings... could reuse some inner gratings if tfsIfiMuliples = 1 (not yet implemented)
    fprintf('Constructing outer apertures and loading inner gratings...')
    apertureTic = tic; 
    load(prerender_file,'inGratingTable','inner_apt_rect','inner_aperture','nInnerFrames');
    inGratPtrTable  = table(repmat({NaN(1,nInnerFrames)},n_unique_parsets,1),'VariableNames',{'inner_gratingtex'}); % inner grating texture pointers
    stimrec_ssmod = [stimrec_ssmod inGratPtrTable]; clear inGratPtrTable
    for AngleInd = 1:length(orient_vals); % draw the inner gratings
        for frame = 1:nInnerFrames
            thisAngle_ptrs(frame) = Screen('MakeTexture', winPtr, inGratingTable.grating{AngleInd,frame});
        end
        match = stimrec_ssmod.orient == orient_vals(1,AngleInd);
        stimrec_ssmod.inner_gratingtex(match) = {thisAngle_ptrs};
    end
    inner_apertex = Screen('MakeTexture', winPtr, inner_aperture);
    clear inGratingTable inner_aperture % free up memory after MakeTexture
else
    fprintf('Constructing outer apertures...');
    apertureTic = tic; 
end

%% Construct apertures for the outer grating of sizes listed in diam_outer. 
for diam = 1:length(diam_vals)
    clear aperture

    % R side of inequality lists aperture radius (radians) required from stimpref.rf_center to contain each pixel
    apt_fullscreen = deg2rad(diam_vals(diam)/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    aperture(:,:,4) = uint8(255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx)))); % draw only the nonzero rect of apt into alpha; RGB=0
    outAptRect = [min(aptx) min(apty) max(aptx) max(apty)];
    match = stimrec_ssmod.diam == diam_vals(diam);
    stimrec_ssmod.outer_apt_rect(match) = {outAptRect};
    stimrec_ssmod.outer_apertex(match) = Screen('MakeTexture', winPtr, aperture);
    if ssmod_stim.prerenderedLocAndDiam     % find the srcRect we need to cut out of the full saved outer grating tex for this diam
        stimrec_ssmod.outerGratSrcRect(match,:) = ...% subtract top and left edges
            repmat(outAptRect-[prerend.outerGratFullRect(1:2) prerend.outerGratFullRect(1:2)], length(find(match)),1); 
    end
end
fprintf([' took ' num2str(toc(apertureTic)) 's.\n']);

%% Creat ssmod_stim.repetitions repeats of the each parameter set and shuffle.
stimrec_ssmod = repmat(stimrec_ssmod,ssmod_stim.repetitions,1);
stimrec_ssmod = stimrec_ssmod(randperm(height(stimrec_ssmod)),:);