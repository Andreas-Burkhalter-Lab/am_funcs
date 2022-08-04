function [stimorder stim_centers grid_center] = rfstim_stimloop(rf_stim, pulsepars, commonvars, grid_center)

%%%% change of angles not working (only presents at anglezero)

%RFSTIM_STIMLOOP % Receptive Field Location Mapping Stimuli
%%% Last updated 1/15/15 on stim computer
% Present circular gratings spaced out in a grid. 
%%%%%%         Input: 
% % 'pars' is the parameters structure for stimuli to be presented. 
% % 'pulsepars' is the parameters structure for stim channel pulseparss.
% % 'grid_center' is the manually (or automatically) determined grid center on the screen. 
%%%%%%         Output:
% % 'stimorder' is the order in which stimuli were presented at each iteration; 
% %     stimorder(1,:) is grid rows, stimorder(2,:) is grid columns. 
% % 'stim_centers' consists of the coordinates of the stimulus centers in the grid.
% % 'grid_center' is the center point of the stimulus grid; it will be
% %    changed in this function if the grid does not fit on the screen.

% A PTB window in which stimuli are to be presented should be open when
%       this function is called. 
% All stimuli will be identical except for location on the screen. 
% .... Gao et al. had 2-degree-diameter patches or gratings, 2s stim
%      presentation


%% Setup
% Note: when a queued 'pulsepars' to the stim channel ends, the stim channel
% voltage remains at last specified voltage of the pulsepars. The command
% 'stop([[[[sessionname]]]])' terminates the pulsepars and leaves voltage at the last
% valued specified before the pulsepars was terminated. 
winPtr = commonvars.winPtr;
global daq_sess
winrect = commonvars.winrect;
V1_dur_scans = pulsepars.V1_dur * daq_sess.Rate;    %% convert V1dur to number of scans
V2_dur_scans = (rf_stim.dur_stim - pulsepars.V1_dur) * daq_sess.Rate;  %% compute V2_dur and convert it to number of scans
Vstimoff_dur_post_scans = pulsepars.Vstimoff_dur_post * daq_sess.Rate; % convert Vstimoff_dur_post to number of scans
dur_stim_scans = round(rf_stim.dur_stim * daq_sess.Rate); %% convert stim duration from seconds to daq toolbox # scans
dur_stim_frames = round(rf_stim.dur_stim / commonvars.ifi);  %% convert stim duration from seconds to PTB screen frames
stimchan_Vstimoff = pulsepars.Vstimoff * ones(Vstimoff_dur_post_scans,1); %% signal for setting stim chan voltage back to Vstimoff
% Angles = linspace(0, 360 - 360/rf_stim.nAngles, rf_stim.nAngles) + rf_stim.zeroAngle; % evenly space angles offset by rf_stim.zeroAngle

% Check that no stimuli in the grid go off of the screen. 'Overshoot'
% values are true if stimuli will be drawn off-screen.
grid_edge.left = grid_center(1) - ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.h_spacing_pix*(rf_stim.columns - 1));
grid_edge.right = grid_center(1) + ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.h_spacing_pix*(rf_stim.columns - 1));
grid_edge.top = grid_center(2) - ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.v_spacing_pix*(rf_stim.rows - 1));
grid_edge.bottom = grid_center(2) + ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.v_spacing_pix*(rf_stim.rows - 1));
grid_overshoot.left = grid_edge.left < winrect(1);
grid_overshoot.right = grid_edge.right > winrect(3);
grid_overshoot.top = grid_edge.top < winrect(2);
grid_overshoot.bottom = grid_edge.bottom > winrect(4);

if grid_overshoot.left
    pixelsToMove.right = ceil(winrect(1) - grid_edge.left);
    input(sprintf(['Grid must be moved %g pixels to the right to fit on the screen.'... 
        ' Press Enter to move the stimulus grid right.'],pixelsToMove.right));
    grid_center(1) = grid_center(1) + pixelsToMove.right;
elseif grid_overshoot.right % should not have both left and right overshoot (grid width checked in rf_mapping)
    pixelsToMove.left = ceil(-(winrect(3) - grid_edge.right));
    input(sprintf(['Grid must be moved %g pixels to the left to fit on the screen.'... 
        ' Press Enter to move the stimulus grid left.'],pixelsToMove.left));
    grid_center(1) = grid_center(1) - pixelsToMove.left;
end
if grid_overshoot.top
    pixelsToMove.down = ceil(winrect(2) - grid_edge.top); 
    input(sprintf(['Grid must be moved %g pixels downward to fit on the screen.'... 
        ' Press Enter to move the stimulus grid downward.'],pixelsToMove.down));
    grid_center(2) = grid_center(2) + pixelsToMove.down;
elseif grid_overshoot.bottom % should not have both top and bottom overshoot (grid height checked in rf_mapping)
    pixelsToMove.up = ceil(-(winrect(4) - grid_edge.bottom)); 
    input(sprintf(['Grid must be moved %g pixels upward to fit on the screen.'... 
        ' Press Enter to move the stimulus grid upward.'],pixelsToMove.up));
    grid_center(2) = grid_center(2) - pixelsToMove.up;
end

% Readjust values of grid_edge.
grid_edge.left = grid_center - ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.h_spacing_pix*(rf_stim.columns - 1));
grid_edge.right = grid_center + ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.h_spacing_pix*(rf_stim.columns - 1));
grid_edge.top = grid_center - ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.v_spacing_pix*(rf_stim.rows - 1));
grid_edge.bottom = grid_center + ceil(0.5*rf_stim.diam_pix + 0.5*rf_stim.v_spacing_pix*(rf_stim.rows - 1));

% Assign screen locations and stim channel signals for all stimuli. 
stimchan_signal = cell(rf_stim.rows, rf_stim.columns); % cell array for holding daq signals identifying row and column of stimulus on grid
stim_centers = cell(rf_stim.rows, rf_stim.columns);  %% center location of stimuli

for row = 1:size(stimchan_signal,1)
    for column = 1:size(stimchan_signal,2)
        % unique stim channel identifier... row = fast neg deflection, column = pos deflection until stim end
        stimchan_signal{row,column} =  [-row*ones(V1_dur_scans,1); column*ones(V2_dur_scans,1)]; 
        stim_centers{row,column} =...   %% assign grid location on the screen
            [grid_center(1) + rf_stim.h_spacing_pix * (column - rf_stim.columns/2 - 0.5),... % x center
            grid_center(2) + rf_stim.v_spacing_pix * (row - rf_stim.rows/2 - 0.5)]; % y center
    end 
end

% Determine random order in which to send stimuli. 
% First row of stimorder is the row index of the stimulus, second row is
% the column index of the stimulus. Stimuli will be displayed in the order
% of the columns in stimorder. 
[stimorder(1,:) stimorder(2,:)] = ind2sub([rf_stim.rows rf_stim.columns],...   %% remap to indices
    randperm(rf_stim.rows * rf_stim.columns)); %% reshuffle 1:number of grid spaces



stim_presentation_start_time = clock; % to save into logfile
rf_stim.logfile = [rf_stim.filetag '_' date];
if commonvars.overwrite_check && (exist(rf_stim.logfile,'file') || exist(strcat(rf_stim.logfile,'.mat'),'file'))
    overwrite = input(sprintf('File ''%s'' already exists. Enter ''y'' to overwrite.\n', rf_stim.logfile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
else
    fprintf('\n Stimulus log information saved into %s.\n',rf_stim.logfile)
end
save(rf_stim.logfile);

input('Begin Merec2 recording on diesel, then press Enter on stimulus computer.');
fprintf(['\nNow presenting RF location-mapping stimuli...\n',...
    'Approximate duration = %g minutes'],...
    (rf_stim.dur_stim+rf_stim.isi)*rf_stim.rows*rf_stim.columns*rf_stim.repetitions/60);

%% Present stimuli
%%%%%% need to empirically check that ifi in this loop = previously
%%%%%% measured; if not, maybe pre-render frames, rather than within this
%%%%%% loop.... maybe tic toc every time, warn if time is different than
%%%%%% expected
%%% could compare vbl of first and last screens to some timer collected by
%%% the daq.... use offset to adjust analyze program...
for iteration = 1:rf_stim.repetitions
    stop(daq_sess); %% end any signals that are already being sent to the stim channel
    daq_sess.queueOutputData(stimchan_Vstimoff); % begin iterations at Vstimoff
    daq_sess.startForeground(); % delay period at Vstimoff before starting iteration
    for thisstim = 1:size(stimorder,2)
        % Determine destinationRect for this stimulus.
        rect_wh = round([rf_stim.diam_pix rf_stim.diam_pix]); % large enough to not cut off the edge of the circle
        destinationRect = CenterRectOnPoint([0 0 rect_wh],...
            stim_centers{stimorder(1,thisstim),stimorder(2,thisstim)}(1),...  %% x location of stim center
            stim_centers{stimorder(1,thisstim),stimorder(2,thisstim)}(2));    %% y location of stim center

        % Create the grating and first screen. 
        % Divide rf_stim.diam_pix by 2 to input as radius.
        % Support width and height must = rect_wh values or else spatial frequency will be incorrect. 
        %%% Use evalc to suppress command-line output.
        evalc(['[gratingid, gratinggrect] = CreateProceduralSineGrating(winPtr, rect_wh(1), rect_wh(2),',...
            'rf_stim.backgroundColorOffset, round(rf_stim.diam_pix/2), rf_stim.contrastPreMultiplicator)']);
        Screen('DrawTexture', winPtr, gratingid, [], destinationRect, rf_stim.orient, [], [],... % do not change angle
            rf_stim.modulateColor,[], kPsychDontDoRotation, [rf_stim.startphase, rf_stim.sf_pix, rf_stim.amp, 0]);    

        % Send signal to stim channel indicating the grid location of thes
        % stimulus and stim onset time. Send as close to first flip time as possible.
        % Send daq pulsepars right after, rather than before the first flip, or else
        % the pulsepars will be sent up to ifi seconds too early. 
        stop(daq_sess); %% end any signals that are already being sent to the stim channel
        daq_sess.queueOutputData(stimchan_signal{stimorder(1,thisstim),stimorder(2,thisstim)});%% queue stim chan signal
        Screen('Flip', winPtr); % flip first frame
        daq_sess.startBackground();    %% send stim channel signal    

        % Send the stimulus in the next grid location. 
        for frame = 2:dur_stim_frames %% flip 2nd and subsequent frames
            Screen('DrawTexture', winPtr, gratingid, [], destinationRect, rf_stim.orient, [], [],... % maybe get vbl
               [rf_stim.modulateColor], [], kPsychDontDoRotation,...
               [rf_stim.startphase+frame*rf_stim.tf_degperflip, rf_stim.sf_pix, rf_stim.amp, 0]); %deg of a sine wave oscillation, not of visual space
            Screen('Flip', winPtr); % maybe use 'DrawingFinished'
        end
        
        % Flip back to black screen when stimulus presentation is done. 
        Screen('FillRect',winPtr,BlackIndex(winPtr)); %% photodiode 'stim-off' signal should go here
        stop(daq_sess); % end stim chan pulsepars if it is still running
        daq_sess.queueOutputData(stimchan_Vstimoff);
        daq_sess.startBackground(); % set stimchan to Vstimoff to indicate stimulus offset    
        Screen('Flip',winPtr);
        pause(pulsepars.Vstimoff_dur_post); %% pause for intertrial intervial; 'pause' should not disrupt daq or ptb
    end
end
    
   
    
end

