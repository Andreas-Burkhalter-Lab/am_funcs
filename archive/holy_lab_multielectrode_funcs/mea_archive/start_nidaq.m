%% Start Nidaq
% Setup to use nidaq usb card, create simple AO signal, send signal. The
% card must be plugged in.

s = daq.createSession('ni');               %% add the device
addAnalogOutputChannel(s,'Dev1','ao0','Voltage');    %% add analogue output channel 1
d = (0:0.001:8)';                           %% create simple AO signal
s.queueOutputData(d);           %% queue output signal
% s.startBackground();       %% run output signal in the background

% %  s.IsContinuous = 1;