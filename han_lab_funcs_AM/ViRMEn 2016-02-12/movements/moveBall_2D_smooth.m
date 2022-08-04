function velocity = moveBall_2D_smooth(vr)
%for use with initializeDAQ_smooth_EH

velocity = [0 0 0 0];
if ~isfield(vr,'scaling')
    %1st is yaw, 2nd is forward
    %vr.scaling = [30 30];%orig
     vr.scaling = [0.10 30];%gains changed. strange. 2.75 rotations forward(1.8m). 360deg is 7.5 yaw
    %vr.scaling = [.8 90];%2.0 rotations forward(1.3m). 360deg is 7.5 yaw
    %vr.scaling = [.8 75];%2.3 rotations forward(1.5m). 360deg is 7.5 yaw

    
end

% Read data from NIDAQ
data = peekdata(vr.ballMovement,50); % (.050s)*(1000samples/s) = 50 samples, so 50ms smoothing window

% Remove NaN's from the data (these occur after NIDAQ has stopped)
f = isnan(mean(data,2));
data(f,:) = [];
data = mean(data,1)';
data(isnan(data)) = 0;

velocity(1) = data(vr.ballForwardChannel)*vr.scaling(1)*-sin(vr.position(4)); %% right CHECK % commented to make 1D
% velocity(4) = data(vr.ballForwardChannel)*vr.scaling(1)*cos(vr.position(4)); % forward CHECK
% velocity(2) = data(vr.ballRotationChannel)*vr.scaling(2); % angular CHECK BUT DRIFT
velocity(2) = data(vr.ballForwardChannel)*vr.scaling(2)*cos(vr.position(4)); % forward CHECK
velocity(4) = data(vr.ballRotationChannel)*(vr.scaling(1)); % angular CHECK BUT DRIFT
% velocity(2) = data(vr.ballRotationChannel)*vr.scaling(1)*cos(vr.position(4));
% velocity(1) = vr.moveData(vr.ballForwardChannel)*vr.scaling(1)*-sin(vr.position(4));
% velocity(4) = -data(vr.ballForwardChannel)*vr.scaling(2);%ptr(2)*vr.scaling(2)/500;

% output x, y, and angle via the DAQ to the clampex computer
outputSingleScan(vr.moveRecordingSession, [(mod(vr.position(4)+pi,2*pi)-pi)*vr.angleScaling, vr.position(1)*vr.xScaling+vr.xOffset, vr.position(2)*vr.yScaling+vr.yOffset]);
