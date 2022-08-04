function [beltMovementDirection pixelsPerCm] = calibrate_measurement_mouse(savefile_append)
%CALIBRATE_MEASUREMENT_MOUSE Determine the mean direction of motion around
%the track (unitMeanDisplacement), then determine the pixel/cm ratio for 
%the measurement mouse. 
% Optionally provide an input filename suffix to save the full
% calibration results.
%%%  last updated 10/2/15

%% Parameters
trackLength = 28; % length once around the track in cm - measured with measuring tape and string 9/14/15
timesAroundTrack = 10; % number of times to rotate the track during calibration
mouseScanRate = 300; % hz
maxDur = 600; % only allocate enough indices in mouse_movements to record for this long in secs; stop recording at this point
xResetDist = 500; % recenter the mouse if it gets farther than this many pixels from horizontal center
yResetDist = 300; % recenter the mouse if it gets farther than this many pixels from vertical center

%% Setup
clear global mouse_movements mousevars
global mouse_movements mousevars
totalIndices = round(maxDur * mouseScanRate);
mouse_movements = NaN(totalIndices,2); % allocate indices to record for maxDur seconds
commandScreenStats = Screen('Resolution',1); % stats from the command screen  - where the measurement mouse cursor will be

mouseTimer = timer; 
mouseTimer.ExecutionMode = 'fixedSpacing';
mouseTimer.Period = round(1000/mouseScanRate)/1000;
mouseTimer.TimerFcn = @(h,~)trackmouse_callback(h);
mouseTimer.TasksToExecute = totalIndices; % stop at this point in case error before stop(mouseTimer)

mousevars.xResetDist = xResetDist;
mousevars.yResetDist = yResetDist;
mousevars.scantime = NaN(size(mouse_movements,1),1); % store the time of each scan in this vector
mousevars.mouseScreenCenter = [commandScreenStats.width/2 commandScreenStats.height/2]; % center of command screen [x y]
mousevars.mouseScanIndex = 0; % counter for mouse scans over the course of this experiment
mousevars.mouseTic = uint64(0); % time between mouse samples... may vary significantly from mouseTimer.Period

%% Run calibration.
input(['Unplug all non-measurement mice, make sure that mouse acceleration is off,'...
    '\n     then press Enter to begin calibration.']) 

mouseAccelerationCheck; % make sure that mouse movement:cursor movement ratio is linear and fixed
SetMouse(mousevars.mouseScreenCenter(1), mousevars.mouseScreenCenter(2));
[mousevars.mousePre(1) mousevars.mousePre(2)] = GetMouse; % needs two output arguments to give x and y
stop(timerfind); % stop any other timers that happen to be running
start(mouseTimer); % begin measuring mouse movements  

input(sprintf(['\nRecording will last for %g seconds.',...
    '\nRotate the track a net %g times, then press Enter.'], maxDur, timesAroundTrack));

%% Finish calibration.
stop(mouseTimer);
if any(mouse_movements(:,1) >= mousevars.mouseScreenCenter(1)) ||...
        any(mouse_movements(:,2) >= mousevars.mouseScreenCenter(2))
    warning('Some mouse displacements were cut off by the screen.')
end
mouse_movements = mouse_movements(1:mousevars.mouseScanIndex,:); % get rid of unused rows in mouse_movements
mousevars.scantime  = mousevars.scantime (1:mousevars.mouseScanIndex,:); % get rid of unused rows

%%% Get unit vector of mean displacement direction so we can later use the
%%% dot product to calculate the projection of a given movement onto this
%%% direction. 
% To get unit vectors:
%%% 1... xOverY = xSum/ySum = xUnit/yUnit...... unit vector keeps same slope
%%% 2... xUnit = xOverY*yUnit
%%% 3... 1^2 = sqrt(xUnit^2 + yUnit^2)..... unit vector has length 1
%%% 4... 1 = xUnit^2 + yUnit^2
%%% 5... 1 = yUnit^2 +(yUnit^2)*xOverY^2....... replacement
%%% 6... 1 = yUnit^2 * (1+xOverY^2)
%%% *7... yUnit = 1 / sqrt(1+xOverY^2)
%%% 8... 1 = xUnit^2 + 1/(1+xOverY)..... replace 7. into 4.
%%% 9... xUnit^2 = 1 - 1/(1+xOverY)
%%% *10.. xUnit = sqrt(1 - 1/(1+xOverY))
xSum = nansum(mouse_movements(:,1)); % net displacement in the x direction in pixels
ySum = nansum(mouse_movements(:,2)); % net displacement in the y direction in pixels
xSign = sign(xSum);
ySign = sign(ySum);
xOverY  = xSum/ySum; % slope of the vector of total displacements
xUnit = sqrt(1 - 1/(1+abs(xOverY))); % abs to avoid imaginary answers
yUnit = 1 / sqrt(1+abs(xOverY)); % abs to avoid imaginary answers
beltMovementDirection = [xSign*xUnit ySign*yUnit]; % unit mean displacement in pixels; apply the sign

% Get pixelsPerCm.
%%% netDisplacement could also be calculated as sqrt(xSum^2 + ySum^2)
%%% because unitMeanDisplacement is parallel to [xSum ySum]. 
netDisplacement = dot([xSum ySum],beltMovementDirection); % pixels; project sum vector onto (parallel) unit vector
traversedDistance = timesAroundTrack * trackLength; % total distance in cm traversed by the mouse in the direction of the track
pixelsPerCm = netDisplacement / traversedDistance; 

% Save results if savefile_append is specified. 
if exist('savefile_append','var') && ~isempty(savefile_append)
    datecheck = clock;
    save(['mouse_calibration_' savefile_append])
end

clear global mouse_movements mousevars

%% use external callback function instead

%% Timer callback function that measures mouse displacement between time intervals
% Uses and modifies variables within the scope of the larger function. 
% % % % function trackmouse_callback(h)
% % % %     if mouseScanIndex > totalIndices
% % % %         error(['Exceeded number of pre-allocated indices in mouse_movements',...
% % % %             '(%g minutes, %g scans)'],maxDur,totalIndices);
% % % %     end
% % % %     
% % % %     mouseScanIndex = mouseScanIndex + 1;
% % % %     [xNow yNow] = GetMouse;
% % % %     timeSincePreviousMouseScan(mouseScanIndex) = toc(mouseTic);
% % % %     mouseTic = tic;
% % % %     
% % % %     % If the mouse moved since the last scan, record this scan number (col 1)
% % % %     % within the trial (col 1) and the x (col 2) and y (col 3) displacements. 
% % % %     if any([xNow yNow] ~= mousePre); % if there was a mouse displacement
% % % %         mouse_movements(mouseScanIndex,1:2) = [xNow yNow]-mousePre;
% % % %     end
% % % %     
% % % %     % Recenter mouse position if necessary and set current mouse position
% % % %     % as mousePre for the next scan. 
% % % %     if abs(xNow - mouseScreenCenter(1)) > xResetDist || abs(yNow - mouseScreenCenter(2)) > yResetDist
% % % %         SetMouse(mouseScreenCenter(1),mouseScreenCenter(2));
% % % %         mousePre = mouseScreenCenter;
% % % %     else
% % % %         mousePre = [xNow yNow];
% % % %     end
% % % %     
% % % % end


end 