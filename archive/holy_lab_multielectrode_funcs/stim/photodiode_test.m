function [] = photodiode_test()
% Compare response speed of photodiode and daq signal

%% DO NOT MODIFY THIS FILE.
% This script was used for the recording 'photodiode_vs_daq.merec' on
% vivid.
%%% daq after flip, single scan

dur_white = 0.2; % seconds
dur_black = 0.2; % seconds
iterations = 4000; 
rect = [0 0 200 200]; % upper left corner 

if isempty(Screen('Windows'))
    win = Screen('OpenWindow', max(Screen('Screens')));
else 
    win = max(Screen('Windows'));
end
    
if ~exist('daq_sess','var')
    daq_sess = daq.createSession('ni');
    addAnalogOutputChannel(daq_sess,'Dev1','ao0','Voltage');  
end
    
for iter = 1:iterations
    Screen('FillRect', win, BlackIndex(win), rect);
    Screen('Flip', win);
    outputSingleScan(daq_sess,0);
    pause(dur_black);
    
    Screen('FillRect', win, WhiteIndex(win), rect);
    Screen('Flip', win);
    outputSingleScan(daq_sess,1);    
    pause(dur_white);
end