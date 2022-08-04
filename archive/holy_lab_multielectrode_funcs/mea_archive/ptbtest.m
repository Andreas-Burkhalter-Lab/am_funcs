%% ptbtest
% Test stimuli from Psychtoolbox.
% Send AI signal to dumbbox before flipping stimulus screen while recording
% screen on dumbbox with photo counter; measure precision of latency
% between the two signals. 

run Ryan_init         %% setup second screen
s = daq.createSession('ni');               %% add nidaq card
addAnalogOutputChannel(s,'Dev1','ao0','Voltage');    %% add analogue output channel 1

loops = 1000;       %% number of test iterations
delay = 1;      %% wait this much time in seconds between trials to clearly separate them
vbl = zeros(loops,1);
stim_onset = zeros(loops,1);
d = [-10; 10];                           %% create simple AO signal
s.queueOutputData(d);           %% queue output signal

for i = 1:loops
    Screen(windowPtr,'FillRect',[255 255 255]);   %% prepare white screen      
    s.startBackground();      %% signal analogue signal
    vbl(i) = Screen(windowPtr,'Flip');      %% flip to white screen
    Screen(windowPtr,'FillRect',[0 0 0]);    %% prepare black screen
    Screen(windowPtr,'Flip');               %% flip to black screen
    pause(delay/2);
    s.queueOutputData(d);           %% queue output signal during delay
    pause(delay/2)
end

Screen('CloseAll')