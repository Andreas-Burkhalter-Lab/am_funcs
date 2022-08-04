% % % function mouse_movements = track_mouse_movements
function track_mouse_movements(calibration_file)
%TRACK_MOUSE_MOVEMENTS Track changes in mouse position for a fixed length
%of time and return timing and amplitude of x and y deflections.
%%% last edited 10/1/15 on stim comp

%% Parameters 
% dur = duration to scan in seconds; if infiniteDuration = 1, this parameter only
%   affects mouse_movements preallocation
%%% multiply expected number of scans per trial by 'sampleBufferFactor value when 
%%%   pre-allocating mouse_movements values
infiniteDuration = 1; % if on, do not end mousetimer until user gives input
dur =500; 
mouseScanRate = 300; % rate of recording mouse position in hz
xResetDist = 500; % recenter the mouse if it gets farther than this many pixels from horizontal center
yResetDist = 300; % recenter the mouse if it gets farther than this many pixels from vertical center
sampleBufferFactor = 1.2; 
savedata = 1; % if on, save data into specified file
    savefile = ['mousemoves_' date]; % if savedata = 1, save into this file

%% Setup
if ~exist('calibration_file','var') || isempty(calibration_file)
    [calb_fname calb_fpath] = uigetfile('*.mat','Select mouse movements calibration file.');
    calibration_file = fullfile(calb_fpath,calb_fname);
end
load(calibration_file,'beltMovementDirection','pixelsPerCm'); % load the main direction of track movement found duration calibration
global mouse_movements mousevars
screenstats = Screen('Resolution',1);

mouse_movements = NaN(round(sampleBufferFactor * dur * mouseScanRate),2);

mousevars.scantime = NaN(size(mouse_movements,1),1); % store the time of each scan in this vector
mousevars.mouseScreenCenter = round([screenstats.width/2 screenstats.height/2]);
mousevars.mouseScanIndex = 0; % counter for mouse scans over course of experiment
[mousevars.mousePre(1) mousevars.mousePre(2)] = GetMouse;
mousevars.xResetDist = xResetDist;
mousevars.yResetDist = yResetDist;

%% add to ss_getevents the mouse_movements results during each trial, save in the .trialdata file

%%% params below were used when mouse tracking was controlled in the same
%%% Matlab window as stimulus presentation; now using 2 Matlab windows,
%%% syncing via 'now'
% mousevars.mouseTic = uint64(0); % time between mouse samples... may vary significantly from mouseTimer.Period
% mousevars.timeSincePreviousMouseScan = NaN(size(mouse_movements,1),1); % recording of time between samples for calculating velocity


mouseTimer = timer;
mouseTimer.ExecutionMode = 'fixedSpacing';
mouseTimer.Period = 1/mouseScanRate;
mouseTimer.TimerFcn = @(h,~)trackmouse_callback(h);
if infiniteDuration
    mouseTimer.TasksToExecute = inf; % keep getting mouse position until stop(mouseTimer)
else
    mouseTimer.TasksToExecute = round(dur/mouseTimer.Period); % end at the specified duration
end

%% Run
mouseAccelerationCheck; % make sure that mouse movement:cursor movement ratio is linear and fixed
SetMouse(mousevars.mouseScreenCenter(1), mousevars.mouseScreenCenter(2));
stop(timerfind);
start(mouseTimer); % begin measuring mouse movements
if infiniteDuration % stop only at user input or after expending the preallocated mouse_movements indices
    stopMouseTimer = [];
    while ~strcmp(stopMouseTimer,'y')
        stopMouseTimer = input('Enter ''y'' to stop recording mouse_movements. ','s');
    end
else % end after specified duration
    pause(dur);
end
stop(mouseTimer); 

%% Postprocessing
%dispInTrackDirection_table = projection of movements movements onto major direction of track movement found during calibration
%%% Do this projection to eliminate side-to-side motion of the belt and
%%% keep only straight forward and backward motion. Positive values indicate forward
%%% motion of the animal (top of track moves backward), negative values 
%%% indicate backward animal motion. Velocity calculations for
%%% 'pixPerSec' will therefore also be the projected value in the major 
%%% belt movement direction.    
% Projection may be time-consuming because of doing dot product on large
% number of scans; could only do projection on nonzero-displacement scans. 
mousePostProcessingTic = tic;
fprintf('Finding velocity in major direction of belt movements... ')
if any(mouse_movements(:,1) >= mousevars.mouseScreenCenter(1)) || any(mouse_movements(:,2) >= mousevars.mouseScreenCenter(2))
    warning('Some mouse displacements were cut off by a screen edge.')
end
mouse_movements = mouse_movements(1:mousevars.mouseScanIndex,:); % get rid of unused rows in mouse_movements
mousevars.scantime = mousevars.scantime(1:mousevars.mouseScanIndex,:); % get rid of unused rows
mouse_movements_table = table(mouse_movements); % so we can use rowfun... may be unnecessarily time-consuming
dispInTrackDirection = table2array(rowfun(@(x)dot(x,beltMovementDirection),mouse_movements_table)); % projected distance
secsSinceLastScan = [NaN; 24*60*60*-diff(mousevars.scantime)]; % NaN for first scan; convert days to seconds
mousePixPerSecond = dispInTrackDirection ./ secsSinceLastScan; % velocity in the major direction of track movement
mouse_movements = table(transpose(1:mousevars.mouseScanIndex),mousevars.scantime,mousePixPerSecond,secsSinceLastScan,...
    mouse_movements(:,1),mouse_movements(:,2),... %% could remove x_dis_pix and y_dis_pix to make file smaller
    'VariableNames',{'scanIndex','scantime','pixPerSec','secsSinceLastScan','x_dis_pix','y_dis_pix'});
nonzero_rows = ~isnan(mouse_movements{:,'x_dis_pix'}) | ~isnan(mouse_movements{:,'x_dis_pix'}); 
mouse_movements = mouse_movements(nonzero_rows,:); % remove rows with no movement
mouse_movements = table2dataset(mouse_movements); % convert to dataset for vivid
lastMouseScan = mousevars.mouseScanIndex;

fprintf('took %g seconds.\n',toc(mousePostProcessingTic))
if savedata
    save(savefile,'mouse_movements','mouseScanRate','beltMovementDirection','pixelsPerCm','lastMouseScan',...
        'calibration_file','savefile');
    fprintf('Mouse movements saved into ''%s''.\n',savefile)
end

clear global mouse_movements mousevars 

    
end