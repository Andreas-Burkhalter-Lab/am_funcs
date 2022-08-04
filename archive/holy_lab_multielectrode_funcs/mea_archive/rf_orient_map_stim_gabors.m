function [stimorder] = rf_orient_map_stim(windowPtr, rect, ifi, daq_sess, pars, pulse, stimcenter)
% Old incomplete version of this function that uses Gabors rather than
% gratings as stimuli. 

%%% Last updated 4/26/2015
%RF_LOC_MAP_STIM % Receptive Field Orientation Mapping Stimuli
% Present patches of different orientations.
% First two args specify the windowPtr and rect of the open screen in which to present stim. 
% 'ifi' is the inter-flip interval computed by psychtoolbox or (for testing only) the nominal refresh interval. 
% 'pars' is the parameters structure for stimuli to be presented. 
% 'daq_sess' is the parameters structure for the daq toolbox session.
% 'pulse' is the parameters structure for stim channel pulses.
% 'stimcenter' is the location on which to center the stimulus grid. 
% output 'stimorder' is the order in which stimuli were presented at each
%   iteration; angle of stim represented by stimorder(x) is 360deg*stimorder(x)/pars.nAngles
% All stimuli will be identical except for orientation.  
% .... Gao et al. had 2-degree-diameter patches or gratings, 2s stim
%      presentation

%% Setup
V1_dur_scans = pulse.V1dur * daq_sess.Rate;    %% convert V1dur to number of scans
V2_dur_scans = (pars.dur_stim - pulseV1dur) * daq_sess.Rate;  %% compute V2_dur and convert it to number of scans
Vbase_dur_post_scans = Vbase_dur_post * daq_sess.Rate; % convert Vbase_dur_post to number of scans
dur_stim_scans = round(pars.dur_stim * daq_sess.Rate); %% convert stim duration from seconds to daq toolbox # scans
dur_stim_frames = round(pars.dur_stim / ifi);  %% convert stim duration from seconds to PTB screen frames
stimchan_Vbase = Vbase * ones(Vbase_dur_post_scans,1); %% signal for setting stim chan voltage back to Vbase

% Assign orientations and stim channel signals for all stimuli. 
stimchan_signal = -inf(V1_dur_scans+V2_dur_scan+Vbase_dur_post_scan, 1); %matrix holding stimchan ID signals; generate error later if not assigned 
stim_angles = -inf(1,nAngles); % initialize; generate error if not properly assigned
for angle_ind = 1:size(stimchan_signal,2) %% create unique stim channel identifier
    stimchan_signal{1,angle_ind} = [ -(floor(angle_ind/10)+1) * ones(V1_dur_scans,1);... %%(tens place +1)=fast neg deflection
        rem(angle_ind-1,10)+1*ones(V2_dur_scans,1)]; %ones place=pos. deflect. until stim end, then back to Vbase
    stim_angles(angle_ind) = pars.angle_range(1) + (angle_ind-1)*diff(pars.angle_range)/(pars.nAngles-1); % both angle_range(1) and (2) get shown
end

stimorder = randperm(pars.nAngles);  % Determine random order in which to send signals. 
destinationRect = CenterRectOnPoint([0 0 pars.rect_wh], stimcenter(1), stimcenter(2)); % get destinationRect for all stimuli

%% Present stimuli
%%%%%% need to empircally check that ifi in this loop = previously
%%%%%% measured; if not, maybe pre-render frames, rather than within this
%%%%%% loop.... maybe tic toc every time, warn if time is different than
%%%%%% expected
for iterations = 1:pars.iterations
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(stimchan_Vbase); % begin iterations at Vbase
    daq_sess.startForeground(); % delay period at Vbase before starting iteration
    for thisstim = 1:length(stimorder)
        angle = stim_angles(stimorder(thisstim)); %% set the angle so we can easily call it during the draw/flip loop

        % Create the gabor and first screen. ............ still need to add the photodiode signal here.
        [gaborid, gaborrect] = CreateProceduralGabor(windowPtr, pars.width, pars.height, pars.nonSymmetric,...
            pars.backgroundColorOffset, pars.disableNorm, pars.contrastPremultiplicator);
        Screen('DrawTexture', windowPtr, gaborid, [], destinationRect, pars.Angle, [], [],...
            [pars.modulateColor],[], kPsychDontDoRotation,...
            [180-pars.phase, pars.freq, pars.sc, pars.contrast, pars.aspectratio, 0, 0, 0]);    

        % Send signal to stim channel indicating the grid location of the
        % stimulus and stim onset time. Send as close to first flip time as possible. 
        stop(daq_sess); %% end any signals that are already being sent to the stim channel
        daq_sess.queueOutputData(stimchan_signal{stimorder(1,thisstim),stimorder(2,thisstim)});%% queue stim channel signal
        daq_sess.startBackground();    %% send stim channel signal    

        % Send the stimulus with the next orientation. &&& ?check whether DrawingFinished helps? 
        Screen('Flip', windowPtr); % flip first frame
        for frame = 2:dur_stim_frames %% draw and flip 2nd and subsequent frames
            Screen('DrawTexture', windowPtr, gaborid, [], destinationRect, angle, [], [],... % angle set above
               [pars.modulateColor], [], kPsychDontDoRotation,...
               [180-pars.phase+frame*pars.speed,pars.freq, pars.sc, pars.contrast, pars.aspectratio, 0, 0, 0]);
            Screen('Flip', windowPtr); 
        end
        
        % Flip back to black screen when stimulus presentation is done. 
        Screen('FillRect',windowPtr,BlackIndex(windowPtr)); %% photodiode 'stim-off' signal shoud go here
        stop(daq_sess); % end stim chan pulse if it is still running
        daq_sess.queueOutputData(stimchan_Vbase);
        daq_sess.startBackground(); % set stimchan to Vbase to indicate stimulus offset    
        Screen('Flip',windowPtr);
        pause(pulse_pars.Vbase_dur_post); %% pause for intertrial intervial; 'pause' should not disrupt daq or ptb
    end
end
    
end

