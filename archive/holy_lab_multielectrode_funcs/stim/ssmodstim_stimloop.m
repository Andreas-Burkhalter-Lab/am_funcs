%SSMODSTIM_STIMLOOP Present stimuli for each trial for ssmodstim_main.
% last updated 1/28/16 on msi

    %       warn if isi is above the specified isi - find minimum
    %       reliable isi
% input(sprintf([mousewarn, 'Approximate presentation duration = %g minutes from %s.'... 
%     '\nBegin Merec2 recording on diesel, then press Enter on stimulus computer.'],...
%     height(stimrec_ssmod)*(ssmod_stim.stimdur+ssmod_stim.isi)/60),datestr(now));
fprintf('\nNow presenting surround-suppression modulation stimuli.\n');
gratingtex = NaN(1,ssmod_stim.nframes);
vbl = NaN(1,ssmod_stim.nframes + 1);
stimPresentationStartTime = clock; %%% approximate stim start time
isiTic = tic;
for trial = 1:height(stimrec_ssmod)
    outgratrow = stimrec_ssmod.grating_table_rows(trial,:); % get the appropriate rows from which to load/draw outer gratings
    nonRepeatedFlips = min([stimrec_ssmod.flipsPerCycle(trial), ssmod_stim.nframes]); % lesser of a full cycle and a full trial
    outer_gratings = cell(1,nonRepeatedFlips); % only load one repetition of the cycle
    if ~ssmod_stim.grats_in_workspace % if we still need to load the gratings into the workspace   
        rowsToLoad = getLinearSequences(outgratrow(1:nonRepeatedFlips),ssmod_stim.maxFrameSeqLength);
        nextFrame = 1; 
        for seqInd = 1:length(rowsToLoad) % load one full sequence of gratings at a time
            seqEnd = nextFrame+length(rowsToLoad{seqInd})-1; % last frame of this sequence
            outer_gratings(nextFrame:seqEnd) = matobj.outGratCells(rowsToLoad{seqInd},1);
            nextFrame = seqEnd + 1; 
        end
    else
        outer_gratings = outGratCells(outgratrow,1);
    end

    % Draw this trial's outer gratings during the isi if ssmod_stim.makeTexDuringIsi=1   
    if ssmod_stim.makeTexDuringIsi
        for frame = 1:nonRepeatedFlips
            gratingtex(frame) = Screen('MakeTexture',winPtr,outer_gratings{frame}); 
        end
        gratingtex = repmat(gratingtex(1:nonRepeatedFlips),1,ceil(ssmod_stim.nframes/nonRepeatedFlips)); %don't need to cut down
    end    
   
    outer_rect_apt = stimrec_ssmod.outer_apt_rect{trial}; 
    %% First flip
if commonvars.send_pulse
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
        % Send a first stim-off pulse for ISI to make sure we start on Vstimoff
        if trial==1
             daq_sess.queueOutputData(pulsepars.stimoff_signal);
             daq_sess.startForeground(); % pause with startForeground to make sure we set to Vstimoff
        end
    daq_sess.queueOutputData(pulsepars.stimon_signal);
end
    
            % Outer Grating
    if ~ssmod_stim.makeTexDuringIsi % MakeTexture between flips
        if ssmod_stim.grats_in_workspace 
            gratingtex(1) = Screen('MakeTexture',winPtr,outGratCells{stimrec_ssmod.grating_table_rows(trial,1)}); 
        else
            gratingtex(1) = Screen('MakeTexture',winPtr,outer_gratings{1}); 
        end
    end
    Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
    Screen('FillRect', winPtr, [0 0 0 0], outer_rect_apt); % maybe take smaller of winrect and this rect
    Screen('DrawTexture', winPtr, stimrec_ssmod.outer_apertex(1), [], outer_rect_apt); % set alpha within aperture to 1
    Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
    if ssmod_stim.prerenderedLocAndDiam        %% srcRect cut from the non-fullscreen prerendered grating
        Screen('DrawTexture', winPtr, gratingtex(1), stimrec_ssmod.outerGratSrcRect(trial,:), outer_rect_apt) 
    else     % srcRect cut from the fullscreen grating
        Screen('DrawTexture', winPtr, gratingtex(1), outer_rect_apt, outer_rect_apt); 
    end
            % Inner Grating
    if ssmod_stim.SHOW_CENTER
        innerframe = 1;
        Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', winPtr, [0 0 0 0], inner_apt_rect); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', winPtr, inner_apertex, [], inner_apt_rect); % set alpha within aperture to 1
        Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        Screen('DrawTexture', winPtr, stimrec_ssmod.inner_gratingtex{trial}(1), inner_apt_rect, inner_apt_rect) % draw the inner grating  
    end
    if trial > 1 % if this is not the first trial
        while toc(trialElapsed) < ssmod_stim.stimdur + ssmod_stim.isi  %% wait to start this trial until stim duration plus isi have elapsed
        end % alternative: if toc < ssmod_stim.stimdur + ssmod_stim.isi; pause(ssmod_stim.stimdur+ssmod_stim.isi - toc)
    elseif trial==1 && ssmod_stim.trackmouse % if this is the first trial 
        SetMouse(mousevars.mouseScreenCenter(1), mousevars.mouseScreenCenter(2));
        [mousevars.mousePre(1) mousevars.mousePre(2)] = GetMouse; % needs two output arguments to give x and y
        stop(timerfind); % stop any other timers that happen to be running
        start(mouseTimer); % begin measuring mouse movements  
    end

    Screen('WaitBlanking',winPtr); % sync the daq pulsepars to the monitor refresh
    if commonvars.send_pulse
        daq_sess.startBackground(); % if pulsepars sent after flip 1, causes 2 missed flips between flip 1 and 2 (takes ~55ms)  
    end
    vbl(1) = Screen('Flip',winPtr); % Send flip 1 and pulsepars as close in time as possible
    stimrec_ssmod.precedingIsi(trial) = toc(isiTic); % record the actual isi to check if we exceeded the desired isi
    trialElapsed = tic; % start timer for this trial after the first flip
    stimrec_ssmod.stimTimeStart(trial) = now; % time in days for comparing to mouse_movements.scantime
    if ssmod_stim.trackmouse
        stimrec_ssmod.mouseScanStimOn(trial) = mousevars.mouseScanIndex; % mouse scan corresponding to beginning of stim presentation
    end
    
    %% Second and subsequent frames of the trial
    for frame = 2:ssmod_stim.nframes % 
        % Outer Grating
        if ~ssmod_stim.makeTexDuringIsi % MakeTexture between flips
            if ssmod_stim.grats_in_workspace 
                %%% next line probably should be using outGratCells, not outGratFrameTable
                gratingtex(frame) = Screen('MakeTexture',winPtr,outGratFrameTable.grating{stimrec_ssmod.grating_table_rows(trial,frame)}); 
            else
                loadedGratingIndex = rem(frame-1,stimrec_ssmod.flipsPerCycle(trial)) + 1; 
                gratingtex(frame) = Screen('MakeTexture',winPtr,outer_gratings{loadedGratingIndex});
            end
        end
        Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', winPtr, [0 0 0 0], outer_rect_apt); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', winPtr, stimrec_ssmod.outer_apertex(trial), [], outer_rect_apt); % set alpha within aperture to 1
        Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        if ssmod_stim.prerenderedLocAndDiam          %% srcRect cut from the non-fullscreen prerendered grating
            Screen('DrawTexture', winPtr, gratingtex(frame), stimrec_ssmod.outerGratSrcRect(trial,:), outer_rect_apt); 
        else     % srcRect cut from fullscreen grating
            Screen('DrawTexture', winPtr, gratingtex(frame), outer_rect_apt, outer_rect_apt); 
        end
        
        % Inner Grating
        if ssmod_stim.SHOW_CENTER
            innerframe = rem(frame-1,nInnerFrames)+1; % cycle back to the first frame if necessary
            Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
            Screen('FillRect', winPtr, [0 0 0 0], inner_apt_rect); % maybe take smaller of winrect and this rect
            Screen('DrawTexture', winPtr, inner_apertex, [], inner_apt_rect); % set alpha within aperture to 1
            Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
            Screen('DrawTexture', winPtr, stimrec_ssmod.inner_gratingtex{trial}(innerframe), inner_apt_rect, inner_apt_rect) % draw the inner grating  
        end
        
        vbl(frame) = Screen('Flip',winPtr);
        
        if toc(trialElapsed) >= ssmod_stim.stimdur % if we're already past the desired stimulus presentation duration... 
%             disp(['quitframe = ' num2str(frame)]) %% for debugging
            break % go to the isi
        end
    end
    
    % Flip back to black screen when stimulus presentation is done and wait for isi. 
    Screen('FillRect', winPtr, BlackIndex(winPtr));
    if commonvars.send_pulse
        stop(daq_sess); % end stimon pulsepars if it is still running
    end
        if ssmod_stim.trackmouse
        stimrec_ssmod.mouseScanStimOff(trial) = mousevars.mouseScanIndex; % mouse scan corresponding to end of stim presentation
    end
    if commonvars.send_pulse
        daq_sess.queueOutputData(pulsepars.stimoff_signal);
    end
    vbl(frame+1) = Screen('Flip',winPtr); % Send black-screen flip and daq pulsepars as close in time as possible.
    isiTic = tic;
    if commonvars.send_pulse
        daq_sess.startBackground(); % set stimchan to stimoff level to indicate stimulus offset
    end
    stimrec_ssmod.stimTimeEnd(trial) = now; % time in days for comparing to mouse_movements.scantime    
    stimrec_ssmod.frameDuration(trial,1:frame) = diff(vbl(1:frame+1)); % time between flipping to a new grating for this trial; takes 0.03ms
    missedFlipsByGrating = round(stimrec_ssmod.frameDuration(trial,:)/ssmod_stim.ifi) - 1;
    stimrec_ssmod.missedFlips(trial) = length(find(missedFlipsByGrating)); % slightly faster than sum
    if stimrec_ssmod.missedFlips(trial)/(frame+1+stimrec_ssmod.missedFlips(trial)) >= commonvars.missedFlipsWarnProportion
        warning('Trial %g - missed %g flips out of %g interflip intervals.',... % warn if missed flips proportion exceeds threshold
            trial, stimrec_ssmod.missedFlips(trial), frame+1+stimrec_ssmod.missedFlips(trial))
    end
    
    Screen('Close',gratingtex(1:nonRepeatedFlips)); % close all outer grating textures opened during this trial... causes error
end