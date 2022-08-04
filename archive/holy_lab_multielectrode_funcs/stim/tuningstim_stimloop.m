%TUNINGSTIM_STIMLOOP Present stimuli for each trial for stimulus preference
%tuning experiments - nonwarped gratings.
%%% Does not currently add values to stimrec (eg. about missed flips, stim
%%% timing)
%%%% Use for orient, sf, tf, and size tuning. 
%%% last updated 1/27/16 on msi
function tuningstim_stimloop(stimrec,stimpars,pulsepars,commonvars)

global daq_sess
nframes = round(stimpars.dur_stim / commonvars.ifi);
vbl = NaN(1,nframes + 1);
stimPresentationStartTime = clock; %%% approximate stim start time
isiTic = tic;

%% Present stimuli
for trial = 1:stimpars.ntrials
    %% First flip
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
        % Send a first stim-off pulse for ISI to make sure we start on Vstimoff
    if trial==1
         daq_sess.queueOutputData(pulsepars.stimoff_signal);
         daq_sess.startForeground(); % pause with startForeground to make sure we set to Vstimoff
    end
    daq_sess.queueOutputData(pulsepars.stimon_signal);
    
         % Determine destinationRect for this stimulus.
    rect_wh = round([stimrec.diam_pix(trial) stimrec.diam_pix(trial)]); % large enough to not cut off the edge of the circle
    destinationRect = CenterRectOnPoint([0 0 rect_wh], stimpars.rf_center(1), stimpars.rf_center(2));
    
    % Create the grating and first screen. 
    % Divide pars.diam_pix by 2 to input as radius.
    % Support width and height must = rect_wh values or else spatial frequency will be incorrect. 
    %%% Use evalc to suppress command-line output.
    evalc(['[gratingid, gratingrect] = CreateProceduralSineGrating(commonvars.winPtr, rect_wh(1), rect_wh(2),',...
        'stimpars.backgroundColorOffset, round(stimrec.diam_pix(trial)/2), stimpars.contrastPreMultiplicator)']);
    Screen('DrawTexture', commonvars.winPtr, gratingid, [], destinationRect, stimrec.orient(trial), [], [],...
        stimpars.modulateColor,[], [], [stimpars.startphase, stimrec.sf_pix(trial), stimpars.amp, 0]);    
    Screen('Flip',commonvars.winPtr);
    
    %% Second through last flips for this trial
    % Present stimulus for this trial.
    for frame = 2:nframes %% flip 2nd and subsequent frames
        Screen('DrawTexture', commonvars.winPtr, gratingid, [], destinationRect, stimrec.orient(trial), [], [],... % maybe get vbl
           [stimpars.modulateColor], [], [],...
           [stimpars.startphase+frame*stimrec.tf_degperflip(trial), stimrec.sf_pix(trial), stimpars.amp, 0]); %deg of a sine wave oscillation, not of visual space
        Screen('Flip', commonvars.winPtr); % maybe use 'DrawingFinished'
    end

    % Flip back to black screen when stimulus presentation is done. 
    Screen('FillRect',commonvars.winPtr,BlackIndex(commonvars.winPtr)); %%
    stop(daq_sess); % end stim chan pulse if it is still running
    daq_sess.queueOutputData(pulsepars.stimoff_signal);
    daq_sess.startBackground(); % set stimchan to Vbase to indicate stimulus offset    
    Screen('Flip',commonvars.winPtr);
    pause(pulsepars.Vstimoff_dur_post); %% pause for intertrial intervial; 'pause' should not disrupt daq or ptb
end