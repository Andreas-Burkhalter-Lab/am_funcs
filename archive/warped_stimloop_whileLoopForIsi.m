%% this version uses a while loop and tic toc to wait for the isi to elapse; future versions use pause

%SSMODSTIM_STIMLOOP Present stimuli for each trial for ssmodstim_main.
% last updated 18/11/17 on thermaltake

    %       warn if isi is above the specified isi - find minimum
    %       reliable isi
% input(sprintf([mousewarn, 'Approximate presentation duration = %g minutes from %s.'... 
%     '\nBegin Merec2 recording on diesel, then press Enter on stimulus computer.'],...
%     height(stimpar_sets)*(stimpars.stimdur+stimpars.isi)/60),datestr(now));

commandline_output = 1; % show params of each trial as it is presented

disp(sprintf(['Now presenting grating stimuli.\n',...
    'Approximate duration = %g minutes'],ntrials*(stimpars.stimdur+stimpars.isi)/60));
gratingtex = NaN(1,stimpars.nframes);
vbl = NaN(1,stimpars.nframes + 1);
stimPresentationStartTime = clock; %%% approximate stim start time
isiTic = tic;
for trial = 1:ntrials
    if commandline_output
        [table({sprintf('%g/%g',trial,ntrials)},'VariableNames',{'trial'}), stimpar_sets(trial,:)]
    end
    if ~stimpar_sets.isblank(trial) % if a grating is to be shown on this trial
        outgratrow = stimpar_sets.grating_table_rows(trial,:); % get the appropriate rows from which to load/draw outer gratings
        nonRepeatedFlips = min([stimpar_sets.flipsPerCycle(trial), stimpars.nframes]); % lesser of a full cycle and a full trial
        outer_gratings = cell(1,nonRepeatedFlips); % only load one repetition of the cycle
        if ~stimpars.grats_in_workspace % if we still need to load the gratings into the workspace   
            rowsToLoad = getLinearSequences(outgratrow(1:nonRepeatedFlips),stimpars.maxFrameSeqLength);
            nextFrame = 1; 
            for seqInd = 1:length(rowsToLoad) % load one full sequence of gratings at a time
                seqEnd = nextFrame+length(rowsToLoad{seqInd})-1; % last frame of this sequence
                outer_gratings(nextFrame:seqEnd) = matobj.outGratCells(rowsToLoad{seqInd},1);
                nextFrame = seqEnd + 1; 
            end
        else
            outer_gratings = outGratCells(outgratrow,1);
        end

        % Draw this trial's outer gratings during the isi if stimpars.makeTexDuringIsi=1   
        if stimpars.makeTexDuringIsi
            for frame = 1:nonRepeatedFlips
                gratingtex(frame) = Screen('MakeTexture',winPtr,outer_gratings{frame}); 
            end
            gratingtex = repmat(gratingtex(1:nonRepeatedFlips),1,ceil(stimpars.nframes/nonRepeatedFlips)); %don't need to cut down
        end    

        outer_rect_apt = stimpar_sets.outer_apt_rect{trial}; 
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
        if ~stimpars.makeTexDuringIsi % MakeTexture between flips
            if stimpars.grats_in_workspace 
                gratingtex(1) = Screen('MakeTexture',winPtr,outGratCells{stimpar_sets.grating_table_rows(trial,1)}); 
            else
                gratingtex(1) = Screen('MakeTexture',winPtr,outer_gratings{1}); 
            end
        end
        Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
        Screen('FillRect', winPtr, [0 0 0 0], outer_rect_apt); % maybe take smaller of winrect and this rect
        Screen('DrawTexture', winPtr, stimpar_sets.outer_apertex(1), [], outer_rect_apt); % set alpha within aperture to 1
        Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
        if stimpars.prerenderedLocAndDiam        %% srcRect cut from the non-fullscreen prerendered grating
            Screen('DrawTexture', winPtr, gratingtex(1), stimpar_sets.outerGratSrcRect(trial,:), outer_rect_apt) 
        else     % srcRect cut from the fullscreen grating
            Screen('DrawTexture', winPtr, gratingtex(1), outer_rect_apt, outer_rect_apt); 
        end
                % Inner Grating
        if stimpars.SHOW_CENTER
            innerframe = 1;
            Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
            Screen('FillRect', winPtr, [0 0 0 0], inner_apt_rect); % maybe take smaller of winrect and this rect
            Screen('DrawTexture', winPtr, inner_apertex, [], inner_apt_rect); % set alpha within aperture to 1
            Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
            Screen('DrawTexture', winPtr, stimpar_sets.inner_gratingtex{trial}(1), inner_apt_rect, inner_apt_rect) % draw the inner grating  
        end
        if trial > 1 % if this is not the first trial
            while toc(trialElapsed) < stimpars.stimdur + stimpars.isi  %% wait to start this trial until stim duration plus isi have elapsed
            end % alternative: if toc < stimpars.stimdur + stimpars.isi; pause(stimpars.stimdur+stimpars.isi - toc)
        end

        Screen('WaitBlanking',winPtr); % sync the daq pulsepars to the monitor refresh
        if commonvars.send_pulse
            daq_sess.startBackground(); % if pulsepars sent after flip 1, causes 2 missed flips between flip 1 and 2 (takes ~55ms)  
        end
        vbl(1) = Screen('Flip',winPtr); % Send flip 1 and pulsepars as close in time as possible
        stimpar_sets.precedingIsi(trial) = toc(isiTic); % record the actual isi to check if we exceeded the desired isi
        trialElapsed = tic; % start timer for this trial after the first flip
        stimpar_sets.stimTimeStart(trial) = now; % time in days

        %% Second and subsequent frames of the trial
        for frame = 2:stimpars.nframes % 
            % Outer Grating
            if ~stimpars.makeTexDuringIsi % MakeTexture between flips
                if stimpars.grats_in_workspace 
                    %%% next line probably should be using outGratCells, not outGratFrameTable
                    gratingtex(frame) = Screen('MakeTexture',winPtr,outGratFrameTable.grating{stimpar_sets.grating_table_rows(trial,frame)}); 
                else
                    loadedGratingIndex = rem(frame-1,stimpar_sets.flipsPerCycle(trial)) + 1; 
                    gratingtex(frame) = Screen('MakeTexture',winPtr,outer_gratings{loadedGratingIndex});
                end
            end
            Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
            Screen('FillRect', winPtr, [0 0 0 0], outer_rect_apt); % maybe take smaller of winrect and this rect
            Screen('DrawTexture', winPtr, stimpar_sets.outer_apertex(trial), [], outer_rect_apt); % set alpha within aperture to 1
            Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
            if stimpars.prerenderedLocAndDiam          %% srcRect cut from the non-fullscreen prerendered grating
                Screen('DrawTexture', winPtr, gratingtex(frame), stimpar_sets.outerGratSrcRect(trial,:), outer_rect_apt); 
            else     % srcRect cut from fullscreen grating
                Screen('DrawTexture', winPtr, gratingtex(frame), outer_rect_apt, outer_rect_apt); 
            end

            % Inner Grating
            if stimpars.SHOW_CENTER
                innerframe = rem(frame-1,nInnerFrames)+1; % cycle back to the first frame if necessary
                Screen('Blendfunction', winPtr, GL_ONE, GL_ZERO, [0 0 0 1]); % draw only to alpha
                Screen('FillRect', winPtr, [0 0 0 0], inner_apt_rect); % maybe take smaller of winrect and this rect
                Screen('DrawTexture', winPtr, inner_apertex, [], inner_apt_rect); % set alpha within aperture to 1
                Screen('Blendfunction', winPtr, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]); % re-enable drawing to all channels
                Screen('DrawTexture', winPtr, stimpar_sets.inner_gratingtex{trial}(innerframe), inner_apt_rect, inner_apt_rect) % draw the inner grating  
            end

            vbl(frame) = Screen('Flip',winPtr);
        end

        % Flip back to black screen when stimulus presentation is done and wait for isi. 
        Screen('FillRect', winPtr, BlackIndex(winPtr));
        if commonvars.send_pulse
            stop(daq_sess); % end stimon pulsepars if it is still running
        end
        if commonvars.send_pulse
            daq_sess.queueOutputData(pulsepars.stimoff_signal);
        end
        vbl(frame+1) = Screen('Flip',winPtr); % Send black-screen flip and daq pulsepars as close in time as possible.

        if commonvars.send_pulse
            daq_sess.startBackground(); % set stimchan to stimoff level to indicate stimulus offset
        end
        stimpar_sets.stimTimeEnd(trial) = now; % time in days for comparing to mouse_movements.scantime    
        stimpar_sets.frameDuration(trial,1:frame) = diff(vbl(1:frame+1)); % time between flipping to a new grating for this trial; takes 0.03ms
        missedFlipsByGrating = round(stimpar_sets.frameDuration(trial,:)/stimpars.ifi) - 1;
        stimpar_sets.missedFlips(trial) = length(find(missedFlipsByGrating)); % slightly faster than sum
        if stimpar_sets.missedFlips(trial)/(frame+1+stimpar_sets.missedFlips(trial)) >= commonvars.missedFlipsWarnProportion
            warning('Trial %g - missed %g flips out of %g interflip intervals.',... % warn if missed flips proportion exceeds threshold
                trial, stimpar_sets.missedFlips(trial), frame+1+stimpar_sets.missedFlips(trial))
        end
        isiTic = tic;
        pause(stimpars.isi); %%% pause for the interstimulus interval
        Screen('Close',gratingtex(1:nonRepeatedFlips)); % close all outer grating textures opened during this trial... causes error?
    elseif stimpar_sets.isblank(trial) % if screen is to remain blank for this trial
        if commonvars.send_pulse
            stop(daq_sess); %% end any signals that are already being sent to the stim channel
            daq_sess.queueOutputData(pulsepars.stimon_signal);
            daq_sess.startBackground(); % indicate 'start' of blank trial
        end
        pause(stimpars.stimdur); % pause instead of showing stimulus
    end
end
if commonvars.send_pulse
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(pulsepars.stimoff_signal);
    daq_sess.startBackground(); % set stimchan to stimoff level to indicate stimulus offset    
end
pause(stimpars.isi); 