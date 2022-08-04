%SSMODSTIM_PREPARESTIM: Create and/or load stimulus gratings for ssstim_main.
% Access the prerendered gratings, load variables, and put prerendered parameters into 'stimpars'.
% Note: with large gratings files, it appears faster to load data (aside
% from gratings) rather than accessing it via matfile. 
% last updated 19/04/25 on thermaltake

fprintf(['Loading prerendered parameters and memory-mapping ''' prerender_file '''...'])
matfileTic = tic; % for timing how long matfile() takes
load(prerender_file,'perp2meshx','perp2meshy',...
    'xstraight','ystraight','eye2screen_center','eye2screen_center_pix','scrnperpx','scrnperpy',...
    'sf_vals','tf_vals','orient_vals','prerend','outGratTrialTable');
matobj = matfile(prerender_file);
fprintf([' took ' num2str(toc(matfileTic)) 's.\n']);

% Copy parameters from the prerended gratings file; copy without overwrite
% fields from defaultstimpars.
prerend_fields = fieldnames(prerend);
for f = 1:length(prerend_fields) % analysis looks for a struct named 'stimpars'
    stimpars = setfield(stimpars, prerend_fields{f}, getfield(prerend,prerend_fields{f}));
end

% Check that prerendered diam and stim centers match those from run_experiment.   
if stimpars.prerenderedLocAndDiam && (stimpars.max_diam ~= stimpars.diam_minmax(2)) 
    error('Specified largest diameter must match largest diameter of prerendered gratings.')
end

% Check that all stimuli fit on screen and shrink if necessary. 
% edgeVals are values from the specified screen edge listing the aperture
%     radius (in radians) required from stimpars.stim_center to contain each pixel.
perp2stimcenterx = scrnperpx - stimpars.stim_center(1); % x pixels between screen center and stim center
perp2stimcentery = scrnperpy - stimpars.stim_center(2); % y pixels between screen center and stim center
edgeVals.left = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(:,1) + perp2stimcentery*perp2meshy(:,1)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(:,1).^2+perp2meshy(:,1).^2)) );
edgeVals.right = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(:,end) + perp2stimcentery*perp2meshy(:,end)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(:,end).^2+perp2meshy(:,end).^2)) );
edgeVals.top = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(1,:) + perp2stimcentery*perp2meshy(1,:)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(1,:).^2+perp2meshy(1,:).^2)) );
edgeVals.bottom = acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx(end,:) + perp2stimcentery*perp2meshy(end,:)) ./...
    sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx(end,:).^2+perp2meshy(end,:).^2)) );
radMaxAllowable = min([edgeVals.left; edgeVals.right; edgeVals.top'; edgeVals.bottom']); % largest radius that can fit on screen

if stimpars.SHOW_CENTER && deg2rad(stimpars.diam_inner/2) > radMaxAllowable % if inner diam does not fit on screen
    error(sprintf(['Inner grating does not fit on screen. Diameters must be \n    ',...
        num2str(2*rad2deg(radMaxAllowable)) ' degrees or smaller to fit on screen.'])) % convert to diameter radians
end
if ~stimpars.allow_offscreen_stim && deg2rad(max(stimpars.diam_minmax)/2) > radMaxAllowable
    commandwindow;
    shrink_ok = input(sprintf(['Largest gratings (%g deg) are cut off by screen. ',... % give option to automatically shrink
        'Enter ''y'' to reduce maximum grating diameter',...
        '\n      from %g degrees to %g degrees.\n'],...
        max(stimpars.diam_minmax), max(stimpars.diam_minmax), 2*rad2deg(radMaxAllowable)),'s');
    if strcmp(shrink_ok,'y')
        stimpars.diam_minmax(2) = 2*rad2deg(radMaxAllowable);
    else
        error('Will not shrink cut-off gratings.')
    end
end

% Create the distributions of outer-grating annulus sizes.
%%% maybe Vaiceliunaite-type distribution of diams is better than exponential; seems
%%% to capture steeper right-side slope which Gao doesn't see
log_diam_minmax = log(stimpars.diam_minmax);
log_diam_range = linspace(log_diam_minmax(1), log_diam_minmax(2), stimpars.n_diams);   
diam_vals = exp(log_diam_range);

if commonvars.send_pulse
    % Stim parameter checks, get pulse signals.
    global daq_sess
    stimpars_tocheck = stimpars;
    stimpars_tocheck.diam_vals = diam_vals;
    stimpars_tocheck.sf_vals = sf_vals;
    stimpars_tocheck.tf_vals = tf_vals;
    stimpars_tocheck.orient_vals = orient_vals;
    pulsepars = check_stim_get_pulse(stimpars_tocheck,pulsepars,daq_sess,commonvars);
end
    
% Table to contain stimulus parameters. Each row is a trial.
% trialTime variables are for synchronizing with mouse tracking timing. 
% Vivid can read dataset objects but not table objects. 
% Create struct with element for every trial, then shuffle trials to randomize order of parameter combinations. 
sfs_plus_tfs_plus_orients = stimpars.n_tfs + stimpars.n_sfs + stimpars.n_orients; % get sum before clearing
ntrials = stimpars.repetitions * [stimpars.n_diams*sfs_plus_tfs_plus_orients  + double(stimpars.include_blank_trials)]; % number of stimulus parameter combinations plus blank trials, including repetitions

stimpar_sets = dataset({NaN(ntrials,4),'orient','sfreq','tfreq','diam'}); 
stimpar_sets.isblank = false(ntrials,1);

% Clear unused variables
if stimpars.n_tfs <= 0 % if we aren't testing a range of temporal frequencies 
    stimpars = rmfield(stimpars, {'n_tfs' });
end
if stimpars.n_sfs <= 0 % if we aren't testing a range of spatial frequencies 
    stimpars = rmfield(stimpars, {'n_sfs' });
end
if stimpars.n_orients <= 0 % if we aren't testing a range of orientations 
    stimpars = rmfield(stimpars, {'n_orients'});
end

%%% randomize trial order, making sure we don't start on a blank trial (causes problems with psychtoolbox)
first_stim_not_blank = false; % initialize
while ~first_stim_not_blank
    randorder = randperm(ntrials); % randomized order in which present the stimulus parameter sets
    first_stim_not_blank = ~any(randorder(end-stimpars.repetitions+1:end) == 1); %%% if index 1 is not present in the block from which blank trials will be assigned
end
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
            end
        end
        if isfield(stimpars, 'n_tfs')
            for ind_tfreq = 1:stimpars.n_tfs % vary tfreq
                combo = combo + 1;   % assign to next element (trial) of parameters structure  
                stimpar_sets.orient(randorder(combo)) = stimpars.orient_fixed;
                stimpar_sets.diam(randorder(combo)) = diam_vals(ind_diam);    %% for reference during later analysis              
                stimpar_sets.sfreq(randorder(combo)) = stimpars.sf_fixed;    %% for reference during later analysis   
                stimpar_sets.tfreq(randorder(combo)) = tf_vals(ind_tfreq); % varying tf  
            end
        end
        if isfield(stimpars,'n_orients')
            for ind_orient = 1:stimpars.n_orients
                combo = combo + 1;
                stimpar_sets.orient(randorder(combo)) = orient_vals(ind_orient);
                stimpar_sets.diam(randorder(combo)) = diam_vals(ind_diam);
                stimpar_sets.sfreq(randorder(combo)) = stimpars.sf_fixed;
                stimpar_sets.tfreq(randorder(combo)) = stimpars.tf_fixed;
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
    end
end
grating_table_rows_table = table(NaN(ntrials,size(outGratTrialTable.grating_table_rows,2)),'VariableNames',{'grating_table_rows'});
outer_apertex_table = table(NaN(ntrials,1),'VariableNames',{'outer_apertex'}); % outer aperture texture pointers
outer_aperrect_table = table(cell(ntrials,1),'VariableNames',{'outer_apt_rect'}); % outer aperture rect values
outerGratSrcRect_table = table(NaN(ntrials,4),'VariableNames',{'outerGratSrcRect'}); % outer grating srcRect values
stimpar_sets = [dataset2table(stimpar_sets) outer_apertex_table outer_aperrect_table outerGratSrcRect_table grating_table_rows_table];
stimpar_sets.flipsPerCycle = 1 ./ (commonvars.ifi .* stimpar_sets.tfreq);

% Get outer grating frame indices for each trial
for trial = 1:ntrials
    if ~stimpar_sets.isblank(trial) % if a grating is shown on this trial and therefore we need to fetch a prerendered grating
        match = stimpar_sets.orient(trial) == outGratTrialTable.orient &...
                stimpar_sets.sfreq(trial) == outGratTrialTable.sf &...
                stimpar_sets.tfreq(trial) == outGratTrialTable.tf;
        stimpar_sets.grating_table_rows(trial,:) = outGratTrialTable.grating_table_rows(match,:);
    end
end

%% Inner Gratings
if stimpars.SHOW_CENTER && ~stimpars.prerenderedLocAndDiam % should probably always prerender inner gratings rather than use this method
    stimpars.diam_inner = input(['Inner gratings not prerendendered in ' prerender_file...
        '\n Type in inner grating diameter (in degrees) then press Enter.']); % set diam_inner here, rather than as a param in run_experiment
    stimpars.amp_inner = input('Type in inner grating amp (max=1) then press Enter.'); % set amp_inner here, rather than as a param in run_experiment
    fprintf('Constructing apertures and inner gratings...')
    apertureTic = tic; 
    
    % Construct aperture for the inner grating. 
    apt_fullscreen = deg2rad(stimpars.diam_inner/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    inner_apt_rect = [min(aptx) min(apty) max(aptx) max(apty)];
    inner_aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of aperture into alpha; RGB=0
    inner_apertex = Screen('MakeTexture', winPtr, inner_aperture);

    % Construct inner gratings of different angles  for each frame. Use the
    % angles, sf_fixed, and tf_fixed already specified by matobj. 
    % (Could prerender these inner gratings along with outer gratings.)
    nInnerFrames = stimpars.nframes; % have not yet implemented reusing frames for this option
    inGratPtrTable  = table(repmat({NaN(1,stimpars.nframes)},n_unique_parsets,1),'VariableNames',{'inner_gratingtex'}); % inner grating texture pointers
    stimpar_sets = [stimpar_sets inGratPtrTable]; clear inGratPtrTable
    for orient = 1:length(orient_vals) % Make the inner gratings. Gratings move along the yrotate dimension.
        spatial_period = 1/stimpars.sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*stimpars.tf_fixed;  % convert from cycles per second to radians per second
        thetarad = deg2rad(orient_vals(1,orient));
        xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
        yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
        eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % angles in radians; pi/2 - acos is faster than asin
        this_orient_ptrs = inf(1,stimpars.nframes);
        for frame = 1:stimpars.nframes 
            grating = 0.5*WhiteIndex(winPtr) + 0.5*WhiteIndex(winPtr)*stimpars.amp_inner*...
                cos(2*pi*sfrad*(eye2screen) - (frame-1)*ifi*tfrad); % from Marshel et al. 2011
            this_orient_ptrs(frame) = Screen('MakeTexture', winPtr, grating);
        end
        match = stimpar_sets.orient == orient_vals(1,orient);
        stimpar_sets.inner_gratingtex(match) = {this_orient_ptrs};
    end
elseif stimpars.SHOW_CENTER && stimpars.prerenderedLocAndDiam  % if using prerendered inner gratings... could reuse some inner gratings if tfsIfiMuliples = 1 (not yet implemented)
    fprintf('Constructing outer apertures and loading inner gratings...')
    apertureTic = tic; 
    load(prerender_file,'inGratingTable','inner_apt_rect','inner_aperture','nInnerFrames');
    inGratPtrTable  = table(repmat({NaN(1,nInnerFrames)},n_unique_parsets,1),'VariableNames',{'inner_gratingtex'}); % inner grating texture pointers
    stimpar_sets = [stimpar_sets inGratPtrTable]; clear inGratPtrTable
    for orient_ind = 1:length(orient_vals) % draw the inner gratings
        for frame = 1:nInnerFrames
            this_orient_ptrs(frame) = Screen('MakeTexture', winPtr, inGratingTable.grating{orient_ind,frame});
        end
        match = stimpar_sets.orient == orient_vals(1,orient_ind);
        stimpar_sets.inner_gratingtex(match) = {this_orient_ptrs};
    end
    inner_apertex = Screen('MakeTexture', winPtr, inner_aperture);
    clear inGratingTable inner_aperture % free up memory after MakeTexture
else
    fprintf('Constructing outer apertures...');
    apertureTic = tic; 
end

%% Construct apertures for the outer grating of sizes listed in diam_outer. 
for diam = 1:length(diam_vals)
    clear aperture apt_fullscreen

    % R side of inequality lists aperture radius (radians) required from stimpars.stim_center to contain each pixel
    apt_fullscreen = deg2rad(diam_vals(diam)/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    aperture(:,:,4) = uint8(255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx)))); % draw only the nonzero rect of apt into alpha; RGB=0
    outAptRect = [min(aptx) min(apty) max(aptx) max(apty)];
    match = stimpar_sets.diam == diam_vals(diam);
    stimpar_sets.outer_apt_rect(match) = {outAptRect};
    stimpar_sets.outer_apertex(match) = Screen('MakeTexture', winPtr, aperture);
    if stimpars.prerenderedLocAndDiam     % find the srcRect we need to cut out of the full saved outer grating tex for this diam
        stimpar_sets.outerGratSrcRect(match,:) = ...% subtract top and left edges
            repmat(outAptRect-[prerend.outerGratFullRect(1:2) prerend.outerGratFullRect(1:2)], length(find(match)),1); 
    end
    if length(diam_vals) == 1
        stimpars.stim_area = apt_fullscreen; % save the stimulus mask
    end
end
fprintf([' took ' num2str(toc(apertureTic)) 's.\n']);