%SSMODSTIM_MAIN Experimental stimuli for testing the modulation of
% surround suppression and temporal and spatial frequency. 
% Use generate_ssmod_stim.m beforehand to generate warped surround stimuli. 
%%% last updated 1/28/16 on msi

%%%%% maybe map 'average' preferred sf and tf of all units rather than
%%%%% starting with fixed value; Vaiceliunaite stared with fixed value   
%%%%%%% possibility: analyze mean tf and sf tuning and hwhm; compare to
%%%%%%% wq's and my collected data to determine likelihood that shank is in
%%%%%%% patch/nonpatch; if not in desired location, penetrate in new area
%%%%%%% (far enough away from first penetratation)

%% Setup
winPtr = commonvars.winPtr;
ifi = commonvars.ifi;
mouseTimer = commonvars.mouseTimer;
commonvars = rmfield(commonvars,'mouseTimer');
global daq_sess

% matfileTic = NaN;
ssmodstim_preparestim; % create and/or load stimulus gratings 

stimrec_ssmod = repmat(stimrec_ssmod, ssmod_stim.repetitions, 1); % copy the stim combinations ssmod_stim.repetitions times
stimrec_ssmod = stimrec_ssmod(randperm(height(stimrec_ssmod)),:); %% shuffle the trial order
stimrec_ssmod_temp = table2dataset([stimrec_ssmod(:,'diam') stimrec_ssmod(:,'orient') stimrec_ssmod(:,'sf') stimrec_ssmod(:,'tf')]); % save now in case trial loop gets interrupted
ssmod_stim.warping = 'warped'; % indicate that warped stimuli were used
ssmod_stim.logfile = [ssmod_stim.filetag '_log_' date];
if commonvars.overwrite_check && (exist(ssmod_stim.logfile,'file') || exist(strcat(ssmod_stim.logfile,'.mat'),'file'))
    overwrite = input(sprintf('File ''%s'' already exists. Enter ''y'' to overwrite.\n', ssmod_stim.logfile),'s');
    if ~strcmp(overwrite,'y')
        error('Will not overwrite file.')
    end
end
save(ssmod_stim.logfile,'stimrec_ssmod_temp','prerender_file','ssmod_stim','pulsepars','commonvars',...
 'diam_vals','sf_vals','tf_vals','orient_vals');

if ssmod_stim.grats_in_workspace
    fprintf('Loading prerendered warped gratings...')
    gratingLoadTic = tic;
    load(prerender_file,'outGratCells');
    fprintf([' took ' num2str(toc(gratingLoadTic)) 's.\n']);
end

%% Mouse tracking
% need to figure out which screen measurement mouse is on and Hide this
% mouse - need to turn off mouse accel and make sure that distance
% measurements are reliable at different velocities
%%%% maybe create calibration program for mouse at a given angle - find
%%%% direction of proper motion and pixel/cm ratio

%% need to incorporate directionality (take only projection of motion along a predefined axis)
%%% keep raw mouse movement data, save calibration results and calibration
%%% file name; if doesn't take too long, perform the projection in the
%%% analysis script
% replace while loop with pause based on tic toc remaining time in trial
% (timer should keep going during pause)

if ssmod_stim.trackmouse
    global mouse_movements mousevars 
    mouseAccelerationCheck; % make sure that mouse movement:cursor movement ratio is linear and fixed
    totaldur = height(stimrec_ssmod) * (ssmod_stim.stimdur + ssmod_stim.isi); % total approximate stim presentation duration in seconds
    mouse_movements = NaN(round(commonvars.mouseBufferFactor * totaldur * commonvars.mouseScanRate),2);
    mousevars.timeSincePreviousMouseScan = NaN(size(mouse_movements,1),1); % recording of time between samples for calculating velocity
    mouseTimer.TimerFcn = @(h,~)trackmouse_callback(h);
    mouseTimer.TasksToExecute = round(size(mouse_movements,1)); % stop at this point in case error before stop(mouseTimer)
    mousevars.mouseScreenCenter = commonvars.mouseScreenCenter;
    mousevars.mouseScanIndex = 0; % counter for mouse scans over the course of this experiment
    trial = 0; 
    mousevars.mouseTic = uint64(0); % time between mouse samples... may vary significantly from mouseTimer.Period
    mousewarn = 'UNPLUG ALL MICE EXCEPT FOR MEASEUREMENT MOUSE.\n'; % maybe use a program to check how many mice are plugged in?
else
    mousewarn = [];
end

ssmodstim_stimloop; %% present stimuli

if ssmod_stim.trackmouse %% probably needs editing to be functional
    stop(mouseTimer);
    if any(mouse_movements(:,1) >= mousevars.mouseScreenCenter(1)) || any(mouse_movements(:,2) >= mousevars.mouseScreenCenter(2))
        warning('Some mouse displacements were cut off by the screen.')
    end
    mouse_movements = mouse_movements(1:stimrec_ssmod.mouseScanStimOff(end),:); % get rid of unused rows in mouse_movements
    mousevars.timeSincePreviousMouseScan = mousevars.timeSincePreviousMouseScan(1:stimrec_ssmod.mouseScanStimOff(end),:); % get rid of unused rows
    mousemovediag = sqrt(mouse_movements(:,1).^2 + mouse_movements(:,2).^2); % 2D mouse displacements
    mousePixPerSecond = mousemovediag ./ mousevars.timeSincePreviousMouseScan; 
    mouse_movements = table(mousePixPerSecond,mousevars.timeSincePreviousMouseScan,mouse_movements(:,1),mouse_movements(:,2),mousemovediag,...
        'VariableNames',{'pixPerSec','tSinceLastScan','x_dis','y_dis','dis2d'}); 
    save(ssmod_stim.logfile,'mouse_movements','-append')
end
stimrec_ssmod.precedingIsi(1) = NaN; % first trial has no preceding ISI

stimrec_ssmod.grating_table_rows = []; % clear unnecessary variables from stimrec_ssmod
stimrec_ssmod.outer_apertex = [];
stimrec_ssmod.outer_apt_rect = [];
if ssmod_stim.SHOW_CENTER
    stimrec_ssmod.inner_gratingtex = [];
end
stimrec_ssmod = table2dataset(stimrec_ssmod); % vivid matlab (2012) can read datasets, not tables
save(ssmod_stim.logfile,'stimrec_ssmod','stimPresentationStartTime','-append')
stimrec_ssmod_mm = matfile(ssmod_stim.logfile, 'Writable', true); 
stimrec_ssmod_mm.stimrec_ssmod_temp = []; % delete the temporary stim record
disp(sprintf(['Surround-suppression stimulus presentation complete; plug back in non-measurement mice.',...
    '\nStimulus and mouse movement saved in ' ssmod_stim.logfile '.mat.'...
    '\nEnd recording on acquisition computer.']));