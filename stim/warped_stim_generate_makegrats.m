%%%% warped_stim_generate_makegrats
% called by warped_stim_generate to run grating construction loops
%%%% updated 2019-2-27 on thermaltake
countTrial = countTrial+1;

% set params that apply to all frames of this stimulus
thetarad = deg2rad(thisOrient);
xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % orients in rads; pi/2 - acos is faster than asin
spatial_period = 1/sf; % need to convert to degrees (not 1/deg) before converting to radians
spatial_period_rad = deg2rad(spatial_period); % convert to radians
sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
tfrad = 2*pi*tf;  % convert from cycles per second to radians per second
outGratTrialTable.orient(countTrial) = thisOrient;
outGratTrialTable.sf(countTrial) = sf;
outGratTrialTable.tf(countTrial) = tf;
if tf == prerend.tf_fixed
    outGratTrialTable.flipsPerCycle(countTrial) = flipsPerCycle_fixed;
elseif ~(tf == prerend.tf_fixed) % if tf the parameter being tested on this trial
    outGratTrialTable.flipsPerCycle(countTrial) = flipsPerCycle(tf_ind);
end
        
% generate a grating for each frame of this stimulus=
 for frame = 1:prerend.nframes    
    %% would be faster if prerendering loc/diam to only do the below calculation for the stim region
    % (cut out values of eye2screen ahead of time)
    % -possibly faster general alternative: prerender full screens, then cut out + load stim region after rf mapping
                   %% 
    if prerend.tfsIfiMultiples
        makeGrating = 0; % turn to true if we need a new grating
        if tf == prerend.tf_fixed
            thisTF_unitsPerFlip = unitsPerFlip_fixed;
            thisTF_cycleFractionUnit = cycleFractionUnit_fixed; % cycleFractionUnit = cycles per flip
            thisTF_unitsPercycle = unitsPerCycle_fixed;
        elseif ~(tf == prerend.tf_fixed) % if tf the parameter being tested on this trial
            thisTF_unitsPerFlip = unitsPerFlip(tf_ind);
            thisTF_cycleFractionUnit = cycleFractionUnit(tf_ind); % cycleFractionUnit = cycles per flip
            thisTF_unitsPerCycle = unitsPerCycle(tf_ind);
        end
        cycleUnitsElapsed = (frame-1)*thisTF_unitsPerFlip;
        % cycleFraction = fraction of cycle elapsed
        cycleFraction = (frame-1)*thisTF_unitsPerFlip*thisTF_cycleFractionUnit; % not used for matching because doesn't align perfectly
        if prerend.tfsIfiMultiples && cycleFraction >= 1 % if we've started over at a new cycle, the grating exists already
            cycleUnitsElapsed = mod(cycleUnitsElapsed,unitsPerCycle_fixed);
            cycleFraction = mod(cycleFraction,1); % reset at a new cycle
        end
        matchRow = find(and(and(and(thisOrient == outGratFrameTable.orient,... % check whether this grating already exists
            sf == outGratFrameTable.sf),...
            cycleUnitsElapsed == outGratFrameTable.cycleUnitsElapsed),...
            thisTF_cycleFractionUnit == outGratFrameTable.cycleFractionUnit),1);
        if matchRow % if this grating already exists
            outGratTrialTable.grating_table_rows(countTrial,frame) = matchRow; % point to existing grating
        else   % make a new grating
            makeGrating = 1; 
        end
    else % we do not already have an identical grating
        makeGrating = 1; 
    end
    if makeGrating % grating of this unique specification doesn't exist yet
        countGrat = countGrat+1;
        outGratFrameTable.orient(countGrat) = thisOrient;
        outGratFrameTable.sf(countGrat) = sf;
        outGratFrameTable.tfForDebugging(countGrat) = tf; % for debugging only
        outGratFrameTable.cycleUnitsElapsed(countGrat) = cycleUnitsElapsed;
        outGratFrameTable.cycleFractionUnit(countGrat) = thisTF_cycleFractionUnit;
        outGratFrameTable.cycleFractionElapsed(countGrat) = cycleFraction; % for debugging only
        outGratTrialTable.grating_table_rows(countTrial,frame) = countGrat; % point to new grating

        whiteindex = WhiteIndex(prerend.stimScreen);
        grating = prerend.brightfraction * uint8(0.5*whiteindex +... % brighfraction adjusts brightness of nonzero portions of grating
                0.5*whiteindex*prerend.outer_amp*... 
                cos(2*pi*sfrad*(eye2screen) - (frame-1)*prerend.ifi*tfrad)); % equation from Marshel et al. 2011 - Functional Specialization of Seven Mouse Visual Cortical Areas
        maxPixVal = max(grating(:)); 
        minPixVal = min(grating(:));
        %%% square wave gratings: make top half of pix white, bottom half black
        if prerend.square_wave 
            whitepix = grating > 0.5*maxPixVal;
            grating(whitepix) = maxPixVal;
            grating(~whitepix) = minPixVal;
        end
        
        % save this grating
        if prerend.prerenderedLocAndDiam % if we're only rendering gratings with a specified location and diameter
            apt_fullscreen = deg2rad(prerend.max_diam/2) > acos( ... % calculate where to draw the aperture/grating
              (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
              sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
            [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture/grating   
            matobj.outGratCells(countGrat,1) = {uint8(grating(min(apty):max(apty),min(aptx):max(aptx)))}; %keep only the stim region
        elseif ~prerend.prerenderedLocAndDiam   % write the fullscreen grating
            matobj.outGratCells(countGrat,1) = {grating}; %   write to file
        end
    end
 end

if isvalid(wb)
    waitbar(countTrial/sfs_plus_tfs_plus_orients,wb)
end