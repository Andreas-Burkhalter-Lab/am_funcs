function code = LinearTrack_withNidaqv1
% linearTrackTwoDirections   Code for the ViRMEn experiment linearTrackTwoDirections.
%   code = linearTrackTwoDirections   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)
initializeDAQ(vr);
%Variables

vr.timeSolenoid = 14; %in milliseconds
vr.endZone = 300;
vr.beginZone = 0;

%vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
%vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
%vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
%vr.currentGoal = 1;

vr.numRewards = 0;
vr.startTime = now;


function datamissed(varargin)


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)
if (vr.position(2)>vr.endZone)
    %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
    vr.isReward = 1;
    reward(vr,vr.timeSolenoid);
    %%Teleport added
    vr.position(2) = vr.beginZone;
    vr.dp(:) = 0;
else
    vr.isReward = 0;
end

if vr.isReward
    vr.numRewards = vr.numRewards + 1;  
    
end



% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
vr.window.Dispose;
