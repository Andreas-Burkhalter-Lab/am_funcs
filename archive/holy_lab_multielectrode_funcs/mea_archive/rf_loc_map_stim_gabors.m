function [stimorder] = rf_loc_map_stim(windowPtr, rect, ifi, daq_sess, pars, pulse, screenstats, gridcenter)
% Old incomplete version of this function that uses Gabors rather than
% gratings as stimuli.
%%% Last updated 5/1/2015

%RF_LOC_MAP_STIM % Receptive Field Location Mapping Stimuli

% Present patches spaced out in a grid. 
% First two args specify the windowPtr and rect of the open screen in which to present stim. 
% 'ifi' is the inter-flip interval computed by psychtoolbox or (for testing only) the nominal refresh interval. 
% 'pars' is the parameters structure for stimuli to be presented. 
% 'daq_sess' is the parameters structure for the daq toolbox session.
% 'pulse' is the parameters structure for stim channel pulses.
% 'gridcenter' is the location on which to center the stimulus grid. 
% output 'stimorder' is the order in which stimuli were presented at each
%    iteration; stimorder(1,:) is grid columns, stimorder(2,:) is grid columns
% All stimuli will be identical except for location on the screen. 
% .... Gao et al. had 5-degree-diameter patches or gratings, 2s stim
%      presentation

%% Setup
% Note: when a queued 'pulse' to the stim channel ends, the stim channel
% voltage remains at last specified voltage of the pulse. The command
% 'stop([[[[sessionname]]]])' terminates the pulse and leaves voltage at the last
% valued specified before the pulse is terminated. 
V1_dur_scans = pulse.V1_dur * daq_sess.Rate;    %% convert V1dur to number of scans
V2_dur_scans = (pars.dur_stim - pulse.V1_dur) * daq_sess.Rate;  %% compute V2_dur and convert it to number of scans
Vbase_dur_post_scans = pulse.Vbase_dur_post * daq_sess.Rate; % convert Vbase_dur_post to number of scans
dur_stim_scans = round(pars.dur_stim * daq_sess.Rate); %% convert stim duration from seconds to daq toolbox # scans
dur_stim_frames = round(pars.dur_stim / ifi);  %% convert stim duration from seconds to PTB screen frames
stimchan_Vbase = pulse.Vbase * ones(Vbase_dur_post_scans,1); %% signal for setting stim chan voltage back to Vbase

% Assign screen locations and stim channel signals for all stimuli. 
stimchan_signal = cell(pars.rows, pars.columns); % cell array for holding daq signals identifying row and column of stimulus on grid
stim_centers = cell(pars.rows, pars.columns);  %% center location of stimuli

for row = 1:size(stimchan_signal,1)
    for column = 1:size(stimchan_signal,2)
        % unique stim channel identifier... row = fast neg deflection, column = pos deflection until stim end
        stimchan_signal{row,column} =  [-row*ones(V1_dur_scans,1); column*ones(V2_dur_scans,1)]; 
        stim_centers{row,column} =...   %% assign grid location on the screen
            [gridcenter(1) + pars.horz_spacing * (column - pars.columns/2 - 0.5),... % x center
            gridcenter(2) + pars.vert_spacing * (row - pars.rows/2 - 0.5)]; % y center
    end
end

% Determine random order in which to send stimuli. 
% First row of stimorder is the row index of the stimulus, second row is
% the column index of the stimulus. Stimuli will be displayed in the order
% of the columns in stimorder. 
[stimorder(1,:) stimorder(2,:)] = ind2sub([pars.rows pars.columns],...   %% remap to indices
    randperm(pars.rows * pars.columns)); %% reshuffle 1:number of grid spaces

%% Present stimuli
%%%%%% need to empirically check that ifi in this loop = previously
%%%%%% measured; if not, maybe pre-render frames, rather than within this
%%%%%% loop.... maybe tic toc every time, warn if time is different than
%%%%%% expected
%%% could compare vbl of first and last screens to some timer collected by
%%% the daq.... use offset to adjust analyze program... maybe not necessary
%%% if photodiode timing works... use photodiode signal to indicate stimoff
%%% time as well as stimon
for iteration = 1:pars.iterations
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(stimchan_Vbase); % begin iterations at Vbase
    daq_sess.startForeground(); % delay period at Vbase before starting iteration
    for thisstim = 1:size(stimorder,2)
        % Determine destinationRect for this stimulus.
        destinationRect = CenterRectOnPoint([0 0 pars.rect_wh],...
            stim_centers{stimorder(1,thisstim),stimorder(2,thisstim)}(1),...  %% x location of stim center
            stim_centers{stimorder(1,thisstim),stimorder(2,thisstim)}(2));    %% y location of stim center

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

        % Send the stimulus in the next grid location. 
        Screen('Flip', windowPtr); % flip first frame
        for frame = 2:dur_stim_frames %% flip 2nd and subsequent frames
            Screen('DrawTexture', windowPtr, gaborid, [], destinationRect, pars.Angle, [], [],... % maybe get vbl
               [pars.modulateColor], [], kPsychDontDoRotation,...
               [180-pars.phase+frame*pars.speed,pars.freq, pars.sc, pars.contrast, pars.aspectratio, 0, 0, 0]);
            Screen('Flip', windowPtr); % maybe use 'DrawingFinished'
        end
        
        % Flip back to black screen when stimulus presentation is done. 
        Screen('FillRect',windowPtr,BlackIndex(windowPtr)); %% photodiode 'stim-off' signal shoud go here
        stop(daq_sess); % end stim chan pulse if it is still running
        daq_sess.queueOutputData(stimchan_Vbase);
        daq_sess.startBackground(); % set stimchan to Vbase to indicate stimulus offset    
        Screen('Flip',windowPtr);
        pause(pulse.Vbase_dur_post); %% pause for intertrial intervial; 'pause' should not disrupt daq or ptb
    end
end
    
   
    
end

