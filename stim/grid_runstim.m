function grid_runstim(stimpars, pulsepars, setupvars)

%%%% change of angles not working (only presents at anglezero)

%RFSTIM_STIMLOOP % Receptive Field Location Mapping Stimuli
%%% Last updated 2018-11-27 on thermaltake
% Present circular gratings spaced out in a grid. 

% All stimuli will be identical except for location on the screen. 
% .... Gao et al. had 2-degree-diameter patches or gratings, 2s stim
%      presentation

%%%% updated 2019-1-18 on thermaltake

%% Setup
commandline_output = 1; % show params of each trial as it is presented
commandwindow;
stimpars.warping = 'nonwarped'; % indicate nonwarped stimuli
stimpars.experiment_type = 'rf_mapping';
stimpars.computer = getenv('computername');
if setupvars.save_stim_file
    stimpars.savefile = strcat('stimdata_',strrep(datestr(now),':',';'),'_rf'); % name of stim log file
else
    fprintf('Not saving stim file.');
    stimpars.savefile = '';
end

commonvars = pulse_and_stim_setup(setupvars);
winPtr = commonvars.winPtr;
winrect = commonvars.winrect;

stimdur_frames = round(stimpars.stimdur / commonvars.ifi);  %% convert stim duration from seconds to PTB screen frames

try stimpars = rmfield(stimpars,'grid_center'); end
stimpars.grid_center = flipud(fliplr(stimpars.grid_center_yx)); % psychtoolbox uses xy, so fliplr or flipud
if stimpars.manual_grid_center
    stimpars.grid_center = follow_cursor;
elseif isempty(stimpars.grid_center) % use screen center as grid center if not otherwise specified
    stimpars.grid_center(1) = round(0.5*[winrect(3) - winrect(1)]);
    stimpars.grid_center(2) = round(0.5*[winrect(4) - winrect(2)]);
end
if commonvars.send_pulse
%     clear global daq_sess
    global daq_sess
    pulsepars.stimon_dur_scans = stimpars.stimdur * daq_sess.Rate;    %% get stim-on duration in scans for the stim channel pulsepars
    pulsepars.stimoff_dur_scans = max([stimpars.isi * daq_sess.Rate, 1]); %% get stim-off duration in scans for the stim channel pulsepars
    pulsepars.stimon_signal = pulsepars.Vstimon * ones(pulsepars.stimon_dur_scans,1); % signal to send at stim onset
    pulsepars.stimoff_signal = pulsepars.Vstimoff * ones(pulsepars.stimoff_dur_scans,1); % signal to send at stim offset
    stimdur_scans = round(stimpars.stimdur * daq_sess.Rate); %% convert stim duration from seconds to daq toolbox # scans

end

% Convert units (degs to pix, hz to degs/flip, cycles/deg to pix/cycle). We are assuming that 
% aspect ratio of screen resolution = aspect ratio of physical screen size, and therefore 
% assuming that the same vert and horz diameter can be used for a circular stimulus. 
%%% Note: deg2pixels is only valid for lengths with one end at the screen center.  
stimpars.diam_pix = deg2pixels(stimpars.diam, commonvars); 
stimpars.h_spacing_pix = deg2pixels(stimpars.h_spacing, commonvars);
stimpars.v_spacing_pix = deg2pixels(stimpars.v_spacing, commonvars);
stimpars.sf_pix = 1 / deg2pixels(1/stimpars.sf, commonvars); % spatial period in pixels
stimpars.tf_degperflip = 360 * commonvars.ifi * stimpars.tf; %%%% deg along a sine wave oscillation, not deg in visual space
stimpars.edge_refresh = (stimpars.tf/stimpars.sf_pix); % velocity in pixels per second... maybe warn if too low

% Check that grid is not too large to fit on the screen.
grid_width = stimpars.h_spacing_pix*(stimpars.columns - 1) + stimpars.diam_pix; % pixels
grid_height = stimpars.v_spacing_pix*(stimpars.rows - 1) + stimpars.diam_pix; % pixels
reduce_h_spacing = []; reduce_v_spacing = []; % initialize
if grid_width > winrect(3) - winrect(1) % if grid is too wide
    max_h_spacing = floor([winrect(3)-winrect(1)-stimpars.diam_pix]/stimpars.columns); % pixels
    reduce_h_spacing = input(sprintf(['Stimulus grid is too wide to fit on screen.\n'...
        'Enter ''y'' to reduce horizontal center-to-center stimulus spacing from %g pixels to %g pixels.\n'],...
        stimpars.h_spacing_pix, max_h_spacing),'s');
    if strcmp(reduce_h_spacing,'y')
        stimpars.h_spacing_pix = max_h_spacing;
    else
        error('Not reducing horizontal center-to-center stimulus spacing. Quitting rf_mapping...')
    end
end
if grid_height > winrect(4) - winrect(2) % if grid is too tall
    max_v_spacing = floor([winrect(4)-winrect(2)-stimpars.diam_pix]/stimpars.rows); % pixels
    reduce_v_spacing = input(sprintf(['Stimulus grid is too tall to fit on screen.\n'...
        '   Enter ''y'' to reduce vertical center-to-center stimulus spacing from %g pixels to %g pixels.\n'],...
        stimpars.v_spacing_pix, max_v_spacing),'s');
    if strcmp(reduce_v_spacing,'y')
        stimpars.v_spacing_pix = max_v_spacing;
    else
        error('Not reducing vertical center-to-center stimulus spacing. Quitting rf_mapping...')
    end
        
end

if strcmp(reduce_h_spacing,'y') & strcmp(reduce_v_spacing,'y') % if both hrz and vert spacing had to be reduced
    disp('Stimulus grid will be centered at screen center; skipping manual grid centering.')
else
    if strcmp(reduce_h_spacing,'y')
        fprintf('\n Grid will be horizontally centered at screen horizontal center.');
        stimpars.grid_center(1) = round(winrect(3) - winrect(1));
    elseif strcmp(reduce_v_spacing,'y')
        fprintf('\n Grid will be vertically centered at screen vertical center.');
        stimpars.grid_center(2) = round(winrect(4) - winrect(2));
    end
end

% Check that no stimuli in the grid go off of the screen. 'Overshoot'
% values are true if stimuli will be drawn off-screen.
grid_edge.left = stimpars.grid_center(1) - ceil(0.5*stimpars.diam_pix + 0.5*stimpars.h_spacing_pix*(stimpars.columns - 1));
grid_edge.right = stimpars.grid_center(1) + ceil(0.5*stimpars.diam_pix + 0.5*stimpars.h_spacing_pix*(stimpars.columns - 1));
grid_edge.top = stimpars.grid_center(2) - ceil(0.5*stimpars.diam_pix + 0.5*stimpars.v_spacing_pix*(stimpars.rows - 1));
grid_edge.bottom = stimpars.grid_center(2) + ceil(0.5*stimpars.diam_pix + 0.5*stimpars.v_spacing_pix*(stimpars.rows - 1));
grid_overshoot.left = grid_edge.left < winrect(1);
grid_overshoot.right = grid_edge.right > winrect(3);
grid_overshoot.top = grid_edge.top < winrect(2);
grid_overshoot.bottom = grid_edge.bottom > winrect(4);

if grid_overshoot.left
    pixelsToMove.right = ceil(winrect(1) - grid_edge.left);
    input(sprintf(['\n Grid must be moved %g pixels to the right to fit on the screen.'... 
        ' Press Enter to move the stimulus grid right.'],pixelsToMove.right));
    stimpars.grid_center(1) = stimpars.grid_center(1) + pixelsToMove.right;
elseif grid_overshoot.right % should not have both left and right overshoot (grid width checked in rf_mapping)
    pixelsToMove.left = ceil(-(winrect(3) - grid_edge.right));
    input(sprintf(['\n Grid must be moved %g pixels to the left to fit on the screen.'... 
        ' Press Enter to move the stimulus grid left.'],pixelsToMove.left));
    stimpars.grid_center(1) = stimpars.grid_center(1) - pixelsToMove.left;
end
if grid_overshoot.top
    pixelsToMove.down = ceil(winrect(2) - grid_edge.top); 
    input(sprintf(['\n Grid must be moved %g pixels downward to fit on the screen.'... 
        ' Press Enter to move the stimulus grid downward.'],pixelsToMove.down));
    stimpars.grid_center(2) = stimpars.grid_center(2) + pixelsToMove.down;
elseif grid_overshoot.bottom % should not have both top and bottom overshoot (grid height checked in rf_mapping)
    pixelsToMove.up = ceil(-(winrect(4) - grid_edge.bottom)); 
    input(sprintf(['\n Grid must be moved %g pixels upward to fit on the screen.'... 
        ' Press Enter to move the stimulus grid upward.'],pixelsToMove.up));
    stimpars.grid_center(2) = stimpars.grid_center(2) - pixelsToMove.up;
end

% Readjust values of grid_edge.
grid_edge.left = stimpars.grid_center - ceil(0.5*stimpars.diam_pix + 0.5*stimpars.h_spacing_pix*(stimpars.columns - 1));
grid_edge.right = stimpars.grid_center + ceil(0.5*stimpars.diam_pix + 0.5*stimpars.h_spacing_pix*(stimpars.columns - 1));
grid_edge.top = stimpars.grid_center - ceil(0.5*stimpars.diam_pix + 0.5*stimpars.v_spacing_pix*(stimpars.rows - 1));
grid_edge.bottom = stimpars.grid_center + ceil(0.5*stimpars.diam_pix + 0.5*stimpars.v_spacing_pix*(stimpars.rows - 1));

% Assign screen locations and stim channel signals for all stimuli. 
rows_times_columns = stimpars.rows * stimpars.columns;
ntrials = rows_times_columns * stimpars.repetitions;
nanrows = NaN(rows_times_columns,1);
stimpar_sets = table(nanrows,     nanrows,       nanrows,         nanrows,...
    'VariableNames', { 'row',     'column',  'stim_center_x',  'stim_center_y'});  
                    
stim_centers_yx_pix = cell(stimpars.rows, stimpars.columns);  %% visual representation of center locations of stimuli

stimcount = 0;
for row = 1:stimpars.rows
    for column = 1:stimpars.columns
        stimcount = stimcount+1; 
        stimpar_sets.row(stimcount) = row;
        stimpar_sets.column(stimcount) = column;
        stimpar_sets.stim_center_x(stimcount) = stimpars.grid_center(1) + stimpars.h_spacing_pix * (column - stimpars.columns/2 - 0.5); %% assign grid location on the screen
        stimpar_sets.stim_center_y(stimcount) = stimpars.grid_center(2) + stimpars.v_spacing_pix * (row - stimpars.rows/2 - 0.5); % y center    
        stim_centers_yx_pix{row,column} = [stimpar_sets.stim_center_y(stimcount), stimpar_sets.stim_center_x(stimcount)];
    end 
end
stimpar_sets = repmat(stimpar_sets,stimpars.repetitions,1); % copy stimpar_sets for each repetition
if stimpars.random_order % shuffle order of stim if applicable
    stimpar_sets = stimpar_sets(randperm(ntrials),:);
end

stimpars.stim_centers_yx_pix = stim_centers_yx_pix;
stimpars.stim_presentation_start_time = clock; % to save into logfile
stimpars.commonvars = commonvars;
if commonvars.save_stim_file
    if stimpars.commonvars.send_pulse
        stimpars.pulse = rmfield(pulsepars,{'stimon_signal','stimoff_signal'});
    end
    save(stimpars.savefile,'stimpars','stimpar_sets');
end
fprintf(['\nNow presenting RF location-mapping stimuli...\n',...
    'Approximate duration = %g minutes \n \n'],...
    (stimpars.stimdur+stimpars.isi)*stimpars.rows*stimpars.columns*stimpars.repetitions/60);

%% Present stimuli
for itrial = 1:ntrials
    % if first trial, make sure stim chan starts at stimoff value
    if itrial == 1
        if commonvars.send_pulse
            stop(daq_sess); %% end any signals that are already being sent to the stim channel
            daq_sess.queueOutputData(pulsepars.stimoff_signal); % begin iterations at Vstimoff
            daq_sess.startForeground(); % delay period at Vstimoff before starting iteration
        end
    end
    % Determine destinationRect for this stimulus.
    rect_wh = round([stimpars.diam_pix stimpars.diam_pix]); % large enough to not cut off the edge of the circle
    destinationRect = CenterRectOnPoint([0 0 rect_wh],...
        stimpar_sets.stim_center_x(itrial),...  %% x location of stim center
        stimpar_sets.stim_center_y(itrial));    %% y location of stim center

    % Create the grating and first screen. 
    % Divide stimpars.diam_pix by 2 to input as radius.
    % Support width and height must = rect_wh values or else spatial frequency will be incorrect. 
    %%% Use evalc to suppress command-line output.
    evalc(['[gratingid, gratinggrect] = CreateProceduralSineGrating(winPtr, rect_wh(1), rect_wh(2),',...
        'stimpars.backgroundColorOffset, round(stimpars.diam_pix/2), stimpars.contrastPreMultiplicator)']);
    Screen('DrawTexture', winPtr, gratingid, [], destinationRect, stimpars.orient, [], [],... % do not change angle
        stimpars.modulateColor,[], kPsychDontDoRotation, [stimpars.startphase, stimpars.sf_pix, stimpars.amp, 0]);    

    if stimpars.keypress_between_trials % wait for keypress
        go_on = '-';
        while ~strcmp(go_on,'')
            go_on = input(['[' num2str(stimpar_sets.row(itrial)) ' ' num2str(stimpar_sets.column(itrial)),...
                ']... = row-column grid index of upcoming grating. Press Enter to show grating.'],'s');
        end
    elseif ~stimpars.keypress_between_trials & commandline_output % commandline output
        [table({sprintf('%g/%g',itrial,ntrials)},'VariableNames',{'trial'}), stimpar_sets(itrial,:)]
    end

    if setupvars.send_pulse
        % Send signal to stim channel as close to first flip time as possible.
        % Send daq pulsepars right after, rather than before the first flip, or else
        % the pulsepars will be sent up to ifi seconds too early. 
        stop(daq_sess); %% end any signals that are already being sent to the stim channel
        daq_sess.queueOutputData(pulsepars.stimon_signal);%% queue stim chan signal
        Screen('Flip', winPtr); % flip first frame
        daq_sess.startBackground();    %% send stim channel signal    
    end

    % Send the stimulus in the next grid location. 
    for frame = 2:stimdur_frames %% flip 2nd and subsequent frames
        Screen('DrawTexture', winPtr, gratingid, [], destinationRect, stimpars.orient, [], [],... % maybe get vbl
           [stimpars.modulateColor], [], kPsychDontDoRotation,...
           [stimpars.startphase+frame*stimpars.tf_degperflip, stimpars.sf_pix, stimpars.amp, 0]); %deg of a sine wave oscillation, not of visual space
        Screen('Flip', winPtr); % maybe use 'DrawingFinished'
    end

    % Flip back to black screen when stimulus presentation is done. 
    Screen('FillRect',winPtr,BlackIndex(winPtr)); %% photodiode 'stim-off' signal should go here
    if setupvars.send_pulse
        stop(daq_sess); % end stim chan pulsepars if it is still running
        daq_sess.queueOutputData(pulsepars.stimoff_signal);
        daq_sess.startBackground(); % set stimchan to Vstimoff to indicate stimulus offset    
    end
    Screen('Flip',winPtr);
    if ~stimpars.keypress_between_trials
        pause(stimpars.isi); %% pause for intertrial intervial; 'pause' should not disrupt daq or ptb
    end
end
clear global daq_sess
disp('Grid stim presentation complete.')
