function [  ] = generate_ssmod_stim(  )
%GENERATE_SSMOD_STIM Create warped surround-suppression stimuli and save them
%as .mat files. 
%   Call this function to generate stimulus files before presenting warped
%   grating stimuli. Gratings are drawn to the full screen; to specify location and size,
%   create an aperture in the stimulus presentation script to overlay onto the
%   grating. Full-screen 120-frame presentations take ~700ms to
%   load on the stimulus computer, so use isi of 1s+. 
%%% Last updated 1/28/16 on stim comp

tic
%% Paremeters
%%% inner grating sf and tf are fixed; inner and outer amp is fixed;
%%% gratings always iso-oriented
%%% maybe add orient_offset option
screen_height = 13;   %% height in inches of stimulus screen; width = 23.5 inches
eye2screen_top_bottom = 9; %% distance in inches from eye to both top and bottom of stimulus screen

prerend.prerenderedLocAndDiam = 1; % prerender only the portion of the grating filling the required location and diameter
    prerend.rf_center = [854  676]; % only save gratings centered here if prerenderedLocAndDiam = 1
    prerend.max_diam = 64; % only save gratings within the area of this diameter if prerenderedLocAndDiam = 1
    prerend.diam_inner = 13; %% diameter of inner grating in degrees... only has an effect if prerenderedLocAndDiam = 1
    prerend.amp_inner = 1; %% amplitude of inner sine gratings (1=max)... only has an effect if prerenderedLocAndDiam = 1
prerend.sf_fixed = 0.1; % cpd
prerend.tf_fixed = 2; % hz
prerend.sf_minmax = [0.01 1];      %% range of annulus spatial frequencies in cycles per degree (~[0.01 1.6] in Gao)
prerend.tf_minmax = [0.5 2];       %% range of annulus temporal frequencies in hz (~[0.1 13] in Gao fig 1)
prerend.outer_amp = 1;
prerend.orient_range = [0 360];  %%%% degrees; positive is clockwise
prerend.n_sfs = 3;      % number of values of annulus spatial frequency to present
prerend.n_tfs = 2;      % number of values of annulus temporal frequency to present
prerend.n_orients = 1; 
prerend.stimdur = 2; % seconds
prerend.tfsCommonMultiples = 1; % make all tfs factors of the next-smallest tf to reduce # of prerendered gratings
prerend.tfsIfiMultiples = 1; % set all tfs >= than 1/stimdur to integer values of flips/cycle to reduce # of gratings
stimset_filename = 'stimset_prerender_test_sizetuning'; % save gratings and parameter values into a .mat file of this name

% PTB checks
prerend.prerend_ifi_checked = 0; % check that sampled ifi doesn't deviate from nominal ifi by too much
    allowableIfiDifference = 1e-5; % measured ifi can differ from nominal refresh interval by this many seconds
    ifiSamples = 500; % get this many flip samples for ifi checking
do_sync_tests = 0; % sets 'SkipSyncTests' to 0 during execution of this script; setting 'SkipSyncTests' to 2 skips all tests
ptb_debug_level = 1; % 1 to 4 (4 = most thorough)
ptb_verbosity = 0; % 0 = no output, 1 = critical errors, 2 = warnings, 3 = startup info

%% Make gratings
% Make sure grating spatial period is between 0 and 180 degrees. 
if any(1./prerend.sf_minmax<0 | 1./prerend.sf_minmax>180)
    error('Grating spatial period must be between 0 and 180 degrees (spatial frequency greater than 0.0056).')
end

% Open PTB window, measure IFI
origSyncLevel = Screen('Preference','SkipSyncTests');
Screen('Preference', 'SkipSyncTests', double(~do_sync_tests));
origDebugLevel = Screen('Preference', 'VisualDebugLevel', ptb_debug_level);
origVerbosity = Screen('Preference', 'Verbosity', ptb_verbosity);
stimScreen = max(Screen('Screens'));  % create stimuli for highest-index screen
Screen('CloseAll');
[winPtr, winrect] = Screen('OpenWindow',stimScreen,BlackIndex(stimScreen));  
screenstats = Screen('Resolution', stimScreen);
nominal_refresh = 1/Screen('NominalFrameRate',winPtr);
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
    if prerend.tfsCommonMultiples && max(tf_vals) >= 1/prerend.stimdur % round only largest tf
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
        cycleFractionUnit = 1 ./ flipsPerCycle(tfInd); % the unit for each tf is its own cycles/flip
        unitsPerFlip = ones(size(tf_vals)); % because the unit was customized for this tf
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
        unitsPerFlip_fixed = 1; % because the unit was customized for this tf
        unitsPerCycle_fixed = flipsPerCycle_fixed; 
    end
else
    cycleFractionUnit_fixed = NaN;% provide filler value
end

% Make temporal frequencies integral factors of the next-lowest tf to
% allow us to reuse gratings across multiple tfs. 
if prerend.tfsCommonMultiples
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

prerend.nframes = round(prerend.stimdur/prerend.ifi); 

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
    perp2stimcenterx = scrnperpx - prerend.rf_center(1); % x pixels between screen center and stim center
    perp2stimcentery = scrnperpy - prerend.rf_center(2); % y pixels between screen center and stim center
end

% Make a table for outer-grating parameters. Create parameter values iteratively.
% Rows in outGratTrialTable describe parameters of the trial. outGratTrialTable.grating_table_rows
% point to the rows of outGratFrameTable.grating that contain the appropriate 
% gratings/frames for this trial. (Rows of outGratTrialTable will later be shuffled but outGratFrameTable 
% will not, so we need to keep track of which parameters point to which row of outGratFrameTable.)
% outGratFrameTable.cycle_fraction is the fraction of a full cycle elapsed
% elapsed at this frame. 
trialNans = NaN(length(orient_vals) * (length(sf_vals)+length(tf_vals)), 1);
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
for orient = 1:length(orient_vals) % could calculate these values just once for each orient, then call them for inner and outer
    thetarad = deg2rad(orient_vals(orient));
    xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
    yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
    eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % orients in rads; pi/2 - acos is faster than asin
    for sf = 1:length(sf_vals)
        countTrial = countTrial+1;
        spatial_period = 1/sf_vals(sf); % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*prerend.tf_fixed;  % convert from cycles per second to radians per second
        outGratTrialTable.orient(countTrial) = orient_vals(orient);
        outGratTrialTable.sf(countTrial) = sf_vals(sf);
        outGratTrialTable.tf(countTrial) = prerend.tf_fixed;
        outGratTrialTable.flipsPerCycle(countTrial) = flipsPerCycle_fixed;
            for frame = 1:prerend.nframes    
                %% would be faster if prerendering loc/diam to only do the below calculation for the stim region
                % (cut out values of eye2screen ahead of time)
                % -possibly faster general alternative: prerender full screens, then cut out + load stim region after rf mapping
                               %% 
                if prerend.tfsIfiMultiples
                    makeGrating = 0; % turn to true if we need a new grating
                    cycleUnitsElapsed = (frame-1)*unitsPerFlip_fixed;
                    cycleFraction = (frame-1)*unitsPerFlip_fixed*cycleFractionUnit_fixed; % not used for matching because doesn't align perfectly
                    if prerend.tfsIfiMultiples && cycleFraction >= 1; % if we've started over at a new cycle, the grating exists already
                        cycleUnitsElapsed = mod(cycleUnitsElapsed,unitsPerCycle_fixed);
                        cycleFraction = mod(cycleFraction,1); % reset at a new cycle
                    end
                    matchRow = find(and(and(and(orient_vals(orient) == outGratFrameTable.orient,... % check whether this grating already exists
                        sf_vals(sf) == outGratFrameTable.sf),...
                        cycleUnitsElapsed == outGratFrameTable.cycleUnitsElapsed),...
                        cycleFractionUnit_fixed == outGratFrameTable.cycleFractionUnit),1);
                    if matchRow % if this grating already exists
                        outGratTrialTable.grating_table_rows(countTrial,frame) = matchRow; % point to existing grating
                    else   % make a new grating
                        makeGrating = 1; 
                    end
                else % we probably do not already have an identical gcrating
                    makeGrating = 1; 
                end
                if makeGrating
                    countGrat = countGrat+1;
                    outGratFrameTable.orient(countGrat) = orient_vals(orient);
                    outGratFrameTable.sf(countGrat) = sf_vals(sf);
                    outGratFrameTable.tfForDebugging(countGrat) = prerend.tf_fixed; % for debugging only
                    outGratFrameTable.cycleUnitsElapsed(countGrat) = cycleUnitsElapsed;
                    outGratFrameTable.cycleFractionUnit(countGrat) = cycleFractionUnit_fixed;
                    outGratFrameTable.cycleFractionElapsed(countGrat) = cycleFraction; % for debugging only
                    outGratTrialTable.grating_table_rows(countTrial,frame) = countGrat; % point to new grating
                    if prerend.prerenderedLocAndDiam % if we're only rendering gratings with a specified location and diameter
                        grating = uint8(0.5*WhiteIndex(stimScreen) +...
                            0.5*WhiteIndex(stimScreen)*prerend.outer_amp*... 
                            cos(2*pi*sfrad*(eye2screen) - (frame-1)*prerend.ifi*tfrad)); % from Marshel et al. 2011
                        apt_fullscreen = deg2rad(prerend.max_diam/2) > acos( ... % calculate where to draw the aperture/grating
                          (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
                          sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
                        [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture/grating   
                        matobj.outGratCells(countGrat,1) = {uint8(grating(min(apty):max(apty),min(aptx):max(aptx)))}; %keep only the stim region
                    else   % write the fullscreen grating
                        matobj.outGratCells(countGrat,1) = {uint8(0.5*WhiteIndex(stimScreen) +... % write to file
                            0.5*WhiteIndex(stimScreen)*prerend.outer_amp*... 
                            cos(2*pi*sfrad*(eye2screen) - (frame-1)*prerend.ifi*tfrad))}; % from Marshel et al. 2011
                    end
                end
            end
    end
    for tfInd = 1:length(tf_vals)
        countTrial = countTrial+1;
        spatial_period = 1/prerend.sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*tf_vals(tfInd);  % convert from cycles per second to radians per second
        outGratTrialTable.orient(countTrial) = orient_vals(orient);
        outGratTrialTable.sf(countTrial) = prerend.sf_fixed;
        outGratTrialTable.tf(countTrial) = tf_vals(tfInd);
        outGratTrialTable.flipsPerCycle(countTrial) = flipsPerCycle(tfInd);
            for frame = 1:prerend.nframes
                if prerend.tfsIfiMultiples || prerend.tfsCommonMultiples
                    makeGrating = 0; % turn to true if we need a new grating
                    cycleUnitsElapsed = (frame-1)*unitsPerFlip(tfInd);
                    cycleFraction = (frame-1)*unitsPerFlip(tfInd)*cycleFractionUnit(tfInd); % not used for matching because doesn't align perfectly
                    if prerend.tfsIfiMultiples && cycleFraction >= 1; % if we've started over at a new cycle, the grating exists already
                        cycleUnitsElapsed = mod(cycleUnitsElapsed,unitsPerCycle(tfInd)); % for debugging only
                        cycleFraction = mod(cycleFraction,1); % reset at a new cycle
                    end
                    matchRow = find(and(and(and(orient_vals(orient) == outGratFrameTable.orient,... % check whether this grating already exists
                        prerend.sf_fixed == outGratFrameTable.sf),...
                        cycleUnitsElapsed == outGratFrameTable.cycleUnitsElapsed),...
                        cycleFractionUnit(tfInd) == outGratFrameTable.cycleFractionUnit),1);
                    if matchRow % if this grating already exists
                        outGratTrialTable.grating_table_rows(countTrial,frame) = matchRow; % point to existing grating
                    else   % make a new grating
                        makeGrating = 1; 
                    end
                else % we probably do not already have an identical grating
                    makeGrating = 1; 
                end
                if makeGrating
                    countGrat = countGrat+1;
                    outGratFrameTable.orient(countGrat) = orient_vals(orient);
                    outGratFrameTable.sf(countGrat) = prerend.sf_fixed;
                    outGratFrameTable.tfForDebugging(countGrat) = tf_vals(tfInd); % for debugging only
                    outGratFrameTable.cycleUnitsElapsed(countGrat) = cycleUnitsElapsed;
                    outGratFrameTable.cycleFractionUnit(countGrat) = cycleFractionUnit(tfInd);
                    outGratFrameTable.cycleFractionElapsed(countGrat) = cycleFraction; % for debugging only
                    outGratTrialTable.grating_table_rows(countTrial,frame) = countGrat; % point to new grating
                    if prerend.prerenderedLocAndDiam % if we're only rendering gratings with a specified location and diameter
                        grating = uint8(0.5*WhiteIndex(stimScreen) +...
                            0.5*WhiteIndex(stimScreen)*prerend.outer_amp*... 
                            cos(2*pi*sfrad*(eye2screen) - (frame-1)*prerend.ifi*tfrad)); % from Marshel et al. 2011
                        apt_fullscreen = deg2rad(prerend.max_diam/2) > acos( ... % calculate where to draw the aperture/grating
                          (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
                          sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
                        [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture/grating  
                        matobj.outGratCells(countGrat,1) = {uint8(grating(min(apty):max(apty),min(aptx):max(aptx)))}; %keep only the stim region
                    else   % write the fullscreen grating
                        matobj.outGratCells(countGrat,1) = {uint8(0.5*WhiteIndex(stimScreen) +... % write to file
                            0.5*WhiteIndex(stimScreen)*prerend.outer_amp*... 
                            cos(2*pi*sfrad*(eye2screen) - (frame-1)*prerend.ifi*tfrad))}; % from Marshel et al. 2011
                    end
                end
            end
    end
end
if prerend.prerenderedLocAndDiam % get size of the full saved outer grating; we will cut out parts of this for different diam sizes
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
    for orient = 1:prerend.n_orients; % Make the inner gratings. Gratings move along the yrotate dimension.
        thetarad = deg2rad(orient_vals(1,orient));
        spatial_period = 1/prerend.sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*prerend.tf_fixed;  % convert from cycles per second to radians per second
        xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
        yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
        eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % orients in radians; pi/2 - acos is faster than asin
        for frame = 1:nInnerFrames
            inGratingTable.grating{orient,frame} = uint8(0.5*WhiteIndex(stimScreen) + 0.5*WhiteIndex(stimScreen)*prerend.amp_inner*...
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
end
        

