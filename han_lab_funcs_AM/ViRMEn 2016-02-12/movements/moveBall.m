function velocity = moveBall(vr)

velocity = [0 0 0 0];
if ~isfield(vr,'scaling')
    %1st is yaw, 2nd is forward
    %vr.scaling = [30 30];%orig
     vr.scaling = [1.1 90];%gains changed. strange. 2.75 rotations forward(1.8m). 360deg is 7.5 yaw
     %170206 changed mouse position and gains. yaw gain had changed during
     %experiment.  mysterious. original gains 0.8 70
    %vr.scaling = [.8 90];%2.0 rotations forward(1.3m). 360deg is 7.5 yaw
    %vr.scaling = [.8 75];%2.3 rotations forward(1.5m). 360deg is 7.5 yaw

    
end
%ptr = get(0,'pointerlocation')-scr(3:4)/2;
%scaling is 7.5 rotations is 360 degrees and 2.75 forward is 1.8 m
vr.moveData = inputSingleScan(vr.moveSession);
%x
%velocity(1) = vr.moveData(vr.ballForwardChannel)*vr.scaling(1)*-sin(vr.position(4));
velocity(1) = vr.moveData(vr.ballForwardChannel)*vr.scaling(1)*-sin(vr.position(4));
%y
%velocity(2) = vr.moveData(vr.ballForwardChannel)*vr.scaling(1)*cos(vr.position(4));
velocity(2) = vr.moveData(vr.ballForwardChannel)*vr.scaling(2)*cos(vr.position(4));
%velocity(4) = vr.moveData(vr.ballRotationChannel)*.1;%ptr(2)*vr.scaling(2)/500;
%velocity(4) = vr.moveData(vr.ballRotationChannel);
velocity(4) = vr.moveData(vr.ballRotationChannel)*(vr.scaling(1));
%velocity(4) = vr.moveData(vr.ballRotationChannel)*.1;%problems turning
outputSingleScan(vr.moveSession, [(mod(vr.position(4)+pi,2*pi) - pi)*vr.angleScaling,vr.position(1)*vr.xScaling+vr.xOffset,vr.position(2)*vr.yScaling+vr.yOffset]);
%velocity(1:2) = [cos(vr.position(4)) -sin(vr.position(4)); sin(vr.position(4)) cos(vr.position(4))]*velocity(1:2)';