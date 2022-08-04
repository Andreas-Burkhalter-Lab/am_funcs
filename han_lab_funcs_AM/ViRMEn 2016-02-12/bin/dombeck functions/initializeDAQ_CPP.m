function vr = initializeDAQ_CPP(vr)
%Added JPB BAR, Nidaq data session1
vr.moveSession = daq.createSession('ni');
vr.waterSession = daq.createSession('ni');

%INPUTS
vr.ballForward = addAnalogInputChannel(vr.moveSession, 'Dev1', 'ai0', 'Voltage');
vr.ballRotation = addAnalogInputChannel(vr.moveSession, 'Dev1', 'ai1', 'Voltage');

%OUPUTS
vr.waterReward = addAnalogOutputChannel(vr.waterSession, 'Dev1','ao0','Voltage');
vr.ballYaw = addAnalogOutputChannel(vr.moveSession, 'Dev1','ao1','Voltage');
vr.ballX = addAnalogOutputChannel(vr.moveSession, 'Dev1','ao2','Voltage');
vr.ballY = addAnalogOutputChannel(vr.moveSession, 'Dev1','ao3','Voltage');

%Variables to be used for output
vr.ballForwardChannel = 1;
vr.ballRotationChannel = 2;


% vr.xScaling = 18/str2double(vr.exper.variables.trackWidth);
vr.xScaling = 18/((str2double(vr.exper.variables.trackLength))/3);
vr.xOffset = 0; %Use all of dynamice range -10 V to 10V, offset to start at -9
vr.yScaling = 18/str2double(vr.exper.variables.trackLength);
vr.yOffset = -9;%Use all of dynamice range -10 V to 10V, offset to start at -9
vr.angleScaling = 9/pi;

%vr.highVoltage = 5;
%vr.highVoltage = 9.6;%grass stim sd9 that is triggered by 9.5V
vr.highVoltage = 9.98;%grass stim sd9.VR2 requires higher voltage to trigger
vr.lowVoltage = 0;