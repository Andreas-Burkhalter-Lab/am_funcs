%%%% upper level function for generating warped gratings for specified set of  stim params
% updated 2019-1-16 on thermaltake

%% Make gratings
if exist([getfname(stimset_filename) '.mat'],'file') % delete prerender file of the same name if it already exists
    delete(stimset_filename)
end

% Make sure grating spatial period is between 0 and 180 degrees. 
if any(1./prerend.sf_minmax<0 | 1./prerend.sf_minmax>180)
    error('Grating spatial period must be between 0 and 180 degrees (spatial frequency greater than 0.0056).')
end

% Open PTB window, measure IFI
origSyncLevel = Screen('Preference','SkipSyncTests');
Screen('Preference', 'SkipSyncTests', double(~do_sync_tests));
origDebugLevel = Screen('Preference', 'VisualDebugLevel', ptb_debug_level);
origVerbosity = Screen('Preference', 'Verbosity', ptb_verbosity);
Screen('CloseAll');
[winPtr, winrect] = Screen('OpenWindow',prerend.stimScreen,BlackIndex(prerend.stimScreen));  
screenstats = Screen('Resolution', prerend.stimScreen);
nominal_refresh = 1/Screen('NominalFrameRate',winPtr);
prerend.screen_height = screen_height;
prerend.eye2screen_top_bottom = eye2screen_top_bottom;
if prerend.prerend_ifi_checked
    [prerend.ifi ifi_samp std ] = Screen('GetFlipInterval', winPtr, ifiSamples);
    if abs(prerend.ifi - nominal_refresh) > allowableIfiDifference % if the measured ifi differs too much from the nominal ifi
        error(['Measured interflip interval more than 1e-5 seconds different from nomimal screen refresh interval:',...
            'IFI - Nominal Refresh = %g'],(prerend.ifi - nominal_refresh));
    end
else
    prerend.ifi = nominal_refresh;
    warning('Not computing interflip interval, using ifi = nominal refresh = %g instead.', nominal_refresh);
end
Screen('Close',winPtr);

% Create the distributions of annulus sf, tf, and orients.
% Sf and tf are on log-scales including highest value and orient is on
% a linear scale not including highest value. 
log_sf_minmax = log(prerend.sf_minmax);
log_tf_minmax = log(prerend.tf_minmax);
log_sf_range = linspace(log_sf_minmax(1),log_sf_minmax(2),prerend.n_sfs);
log_tf_range = linspace(log_tf_minmax(1),log_tf_minmax(2),prerend.n_tfs);
sf_vals = exp(log_sf_range);
tf_vals = exp(log_tf_range);
orient_inc = (prerend.orient_range(2) - prerend.orient_range(1)) / prerend.n_orients; % degrees between each orient
orient_vals = prerend.orient_range(1) : orient_inc : prerend.orient_range(2)-orient_inc; % even distribution not including the upper limit

%% Temporal frequency adjustments
%%% Optional 'registration' of tfs to ifi and to each other to allow for
%%% reuse of gratings.
%%% Cycle progression at each flip will be an integer multiple of the respective
%%% value of 'cycleFractionUnit'. 

% Make all tfs that are >= than 1/stimdur integral factors of 
% the ifi to allow us to reuse gratings after a full cycle has elapsed.
% Only round the largest tf if the other tfs are going to be rounded based
% on this tf (if tfsCommonMultiples is on). 
if prerend.tfsIfiMultiples
    prerend.tfValsPreAdjustment = tf_vals; % save original values for reference
    spc_vals = 1./tf_vals; % seconds per cycle
    if prerend.tfsCommonMultiples && ~isempty(tf_vals) && max(tf_vals) >= 1/prerend.stimdur % round only largest tf
        gratsToRound = length(tf_vals); 
    else % round all tfs > 1/stimdur
        [junk gratsToRound] = find(tf_vals >= 1/prerend.stimdur);
        flipsPerCycle = NaN(size(tf_vals));
        cycleFractionUnit = NaN(size(tf_vals));
    end
    for tfInd = gratsToRound
        spcRemainder = rem(spc_vals(tfInd),prerend.ifi); 
        if spcRemainder > 0.5*prerend.ifi
            spc_vals(tfInd) = spc_vals(tfInd) - spcRemainder + prerend.ifi; % round up
        else
            spc_vals(tfInd) = spc_vals(tfInd) - spcRemainder; % round down
        end
    end
    if ~prerend.tfsCommonMultiples % if tfs won't share a common cycle/flip unit, give a unit to each rounded tf
        flipsPerCycle = round(spc_vals ./ prerend.ifi);
        cycleFractionUnit = 1 ./ flipsPerCycle; % the unit for each tf is its own cycles/flip
        unitsPerFlip = ones(size(tf_vals)); % because the unit was customized for this tf and isn't a multiple of any other tf we're using
    end
    tf_vals = 1./spc_vals;
    
    % Make tf_fixed a multiple of the ifi if it's > 1/stimdur
    if prerend.tf_fixed > 1/prerend.stimdur
        prerend.tfFixedPreAdjustment = prerend.tf_fixed; % save original value for reference
        spc_fixed = 1./prerend.tf_fixed; % second per cycle
        spcRemainder = rem(spc_fixed,prerend.ifi);
        if spcRemainder > 0.5*prerend.ifi
            spc_fixed = spc_fixed - spcRemainder + prerend.ifi; % round up
        else
            spc_fixed = spc_fixed - spcRemainder; % round down
        end
        prerend.tf_fixed = 1/spc_fixed;
        flipsPerCycle_fixed = round(spc_fixed ./ prerend.ifi);
        cycleFractionUnit_fixed = 1 ./ flipsPerCycle_fixed; % the unit for each tf is its own cycles/flip
        unitsPerFlip_fixed = 1; % because the unit was customized for this tf and isn't a multiple of any other tf we're using
        unitsPerCycle = flipsPerCycle;
        unitsPerCycle_fixed = flipsPerCycle_fixed; 
    end
else
    cycleFractionUnit_fixed = NaN;% provide filler value
end

% Make temporal frequencies integral factors of the next-lowest tf to
% allow us to reuse gratings across multiple tfs. 
if ~isempty(tf_vals) && prerend.tfsCommonMultiples
    error('check that this option still computes and displays tf correctly')
    if ~isfield(prerend,'tfValsPreAdjustment')
       prerend.tfValsPreAdjustment = tf_vals; % save original values for reference
    end
    spc_vals = 1./tf_vals; % seconds per cycle
    for i = 2:length(spc_vals)
        spcRemainder = rem(spc_vals(end-i+1),spc_vals(end-i+2));
        if  spcRemainder > 0.5 * spc_vals(end-i+2)
            spc_vals(end-i+1) = spc_vals(end-i+1) - spcRemainder + spc_vals(end-i+2);
        else
            spc_vals(end-i+1) = spc_vals(end-i+1) - spcRemainder;
        end
    end
    flipsPerCycle = round(spc_vals ./ prerend.ifi); % make sure it's integer in case prerend.tfsIfiMultiples == 1
    cyclesPerFlip_vals = 1 ./ flipsPerCycle;
    cycleFractionUnit = repmat(cyclesPerFlip_vals(1),1,length(tf_vals)); % same unit for all tfs
    unitsPerFlip = round(cyclesPerFlip_vals ./ cycleFractionUnit);
    unitsPerCycle = round(1 ./ cycleFractionUnit);
    tf_vals = 1./spc_vals;
end

prerend.nframes = round(prerend.stimdur/prerend.ifi); 

if ~(prerend.tfsIfiMultiples || prerend.tfsCommonMultiples) % provide filler values
    flipsPerCycle = repmat(prerend.nframes,length(tf_vals),1); % the whole stimdur is a non-repeated cycle
    flipsPerCycle_fixed = prerend.nframes; % the whole stimdur is a non-repeated cycle
    cycleUnitsElapsed = NaN; cycleFraction = NaN; cycleFractionUnit = NaN(size(tf_vals));
end
    %% 

% Make sure temporal frequencies are positive and don't exceed 0.5/IFI. 
if any(prerend.tf_minmax<0) || any(prerend.tf_minmax>0.5/prerend.ifi)
    error('Stimulus temporal frequency must be between 0 and 0.5/IFI (=%gc/d) or direction may appear reversed.',0.5/ifi)
end

% Set up constant variables for creating gratings.
eye2screen_center = sqrt( eye2screen_top_bottom^2 - (screen_height/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
eye2screen_center_pix = eye2screen_center * (screenstats.height / screen_height); % x0 in Marshel et al; convert to pix
scrnperpx = round(screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
scrnperpy = round(screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
[xcentermesh ycentermesh] = meshgrid(1:screenstats.width,1:screenstats.height);
perp2meshx = scrnperpx-xcentermesh;
perp2meshy = scrnperpy-ycentermesh;
[xstraight ystraight] = meshgrid(-screenstats.width/2 : screenstats.width/2, -screenstats.height/2 : screenstats.height/2);
if prerend.prerenderedLocAndDiam
    perp2stimcenterx = scrnperpx - prerend.stim_center(1); % x pixels between screen center and stim center
    perp2stimcentery = scrnperpy - prerend.stim_center(2); % y pixels between screen center and stim center
end

% Make a table for outer-grating parameters. Create parameter values iteratively.
% Rows in outGratTrialTable describe parameters of the trial. outGratTrialTable.grating_table_rows
% point to the rows of outGratFrameTable.grating that contain the appropriate 
% gratings/frames for this trial. (Rows of outGratTrialTable will later be shuffled but outGratFrameTable 
% will not, so we need to keep track of which parameters point to which row of outGratFrameTable.)
% outGratFrameTable.cycle_fraction is the fraction of a full cycle elapsed
% elapsed at this frame. 
sfs_plus_tfs_plus_orients = prerend.n_sfs + prerend.n_tfs + prerend.n_orients; 
trialNans = NaN(sfs_plus_tfs_plus_orients, 1);
outGratTrialTable = table(trialNans,trialNans,trialNans,trialNans,repmat(trialNans,1,prerend.nframes),...
    'VariableNames',{'orient','sf','tf','flipsPerCycle','grating_table_rows'});
outGratCells = cell(height(outGratTrialTable),prerend.nframes);
frameNans = repmat(trialNans,prerend.nframes,1); % one grating per frame; reduce unused rows later
outGratFrameTable = table(frameNans,frameNans,frameNans,frameNans,frameNans,frameNans,... 
    'VariableNames',{'orient','sf','tfForDebugging','cycleUnitsElapsed','cycleFractionUnit','cycleFractionElapsed'}); % add the gratings later

% Make the outer gratings. 
% Create the prerender file to write gratings into.
if ~prerend.prerenderedLocAndDiam % if drawing to full screen
    prerend = rmfield(prerend,{'stim_center','max_diam','diam_inner'}); %%% remove these fields if they aren't being used
end
save(stimset_filename, '-v7.3','prerend','screenstats','screen_height','eye2screen_top_bottom',...
    'perp2meshx','perp2meshy','xstraight','ystraight','scrnperpx','scrnperpy',...
    'eye2screen_center','eye2screen_center_pix','sf_vals','tf_vals','orient_vals');   
matobj = matfile(stimset_filename, 'Writable', true); % enable saving directly into the file
matobj.outGratCells = cell(size(frameNans)); % need to write into a cell, not a table; transfer to table later
countTrial = 0; % to index into successive outGratTrialTable rows
countGrat = 0; % to index to successively created gratings

%% make gratings for each sf, tf, and orient value
wb = waitbar(0,'Generating warped gratings...');
if prerend.n_sfs > 0
    for sf_ind = 1:length(sf_vals)
        clear sf tf orient spatial_period spatial_period_rad xrotate yrotate eye2screen makeGrating sfrad tfrad...
            grating matchRow cycleUnitsElapsed cycleFraction
        sf = sf_vals(sf_ind);
        tf = prerend.tf_fixed;
        thisOrient = prerend.orient_fixed; %%% 'orient' is a builtin function
        warped_stim_generate_makegrats;
    end
end
if prerend.n_tfs > 0
    for tf_ind = 1:length(tf_vals)
        clear sf tf orient spatial_period spatial_period_rad xrotate yrotate eye2screen makeGrating sfrad tfrad...
            grating matchRow cycleUnitsElapsed cycleFraction
        sf = prerend.sf_fixed;
        tf = tf_vals(tf_ind);
        thisOrient = prerend.orient_fixed;
        warped_stim_generate_makegrats;
    end
end
if prerend.n_orients > 0
    for orient_ind = 1:length(orient_vals)
        clear sf tf orient spatial_period spatial_period_rad xrotate yrotate eye2screen makeGrating sfrad tfrad...
            grating matchRow cycleUnitsElapsed cycleFraction
        sf = prerend.sf_fixed;
        tf = prerend.tf_fixed;
        thisOrient = orient_vals(orient_ind);
        warped_stim_generate_makegrats;
    end
end
if isvalid(wb)
    close(wb)
end

if prerend.prerenderedLocAndDiam % get size of the full saved outer grating; we will cut out parts of this rect for different diam sizes
    prerend.outerGratFullRect = [min(aptx) min(apty) max(aptx) max(apty)];
end
outGratFrameTable(isnan(outGratFrameTable.orient),:) = []; % remove rows which don't describe a grating
matobj.outGratCells = matobj.outGratCells(1:height(outGratFrameTable),:); % remove rows which don't contain a grating
save(stimset_filename,'prerend','outGratFrameTable','outGratTrialTable','-append');  % save the list of parameter combinations and corresponding gratings
    
%% Inner Gratings
% Construct aperture for the inner grating. 
%%% There are few enough inner gratings that we don't need to periodically
%%% write to file.
%%%%% Could reuse some inner gratings if tfsIfiMuliples = 1 (not yet implemented).   
if prerend.prerenderedLocAndDiam
    if prerend.tfsIfiMultiples    
        nInnerFrames = flipsPerCycle_fixed; 
    else
        nInnerFrames = prerend.nframes;
    end
    
    %%% Construct inner aperture.
    apt_fullscreen = deg2rad(prerend.diam_inner/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
        sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
    inner_apt_rect = [min(aptx) min(apty) max(aptx) max(apty)];
    inner_aperture(:,:,4) = uint8(255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx)))); % draw only the nonzero rect of aperture into alpha; RGB=0

    % Construct inner gratings of different orients  for each frame. Use the
    % orients, sf_fixed, and tf_fixed already specified by matobj. 
    inGratingTable  = table(orient_vals',cell(prerend.n_orients,nInnerFrames),...
        'VariableNames',{'orient','grating'}); % inner gratings for each unique frame for each orient
    for thisOrient = 1:prerend.n_orients % Make the inner gratings. Gratings move along the yrotate dimension.
        thetarad = deg2rad(orient_vals(1,thisOrient));
        spatial_period = 1/prerend.sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*prerend.tf_fixed;  % convert from cycles per second to radians per second
        xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
        yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
        eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % orients in radians; pi/2 - acos is faster than asin
        for frame = 1:nInnerFrames
            inGratingTable.grating{thisOrient,frame} = brightfraction * uint8(0.5*WhiteIndex(prerend.stimScreen) + 0.5*WhiteIndex(prerend.stimScreen)*prerend.amp_inner*...
                cos(2*pi*sfrad*(eye2screen) - (frame-1)*prerend.ifi*tfrad)); % from Marshel et al. 2011
        end
    end
    save(stimset_filename,'inGratingTable','inner_apt_rect','inner_aperture','nInnerFrames','-append');
end

%% Cleanup
Screen('Preference','SkipSyncTests',origSyncLevel); % reset to original sync testing setting
Screen('Preference', 'VisualDebugLevel',origDebugLevel); % reset to iriginal PTB debug level
Screen('Preference','Verbosity',origVerbosity); % reset to original PTB warning mode
generationTime = toc; 
save(stimset_filename,'generationTime','-append')
disp(['Generating gratings took ' num2str(generationTime) 's.'])
disp(['Warped gratings drawn into ' [stimset_filename] '.mat'])
commandwindow;
        