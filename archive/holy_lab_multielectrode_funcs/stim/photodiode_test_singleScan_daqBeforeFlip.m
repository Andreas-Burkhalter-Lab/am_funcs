function [] = photodiode_test_singleScan_daqBeforeFlip()
% Compare response speed of photodiode and daq signal
%%% updated 9/21/15

%%% maybe could instead use 'dontclear' to synchronize daq to the retrace
%%% rather than (deprecated) WaitBlanking

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
    Screen('WaitBlanking',win);
    outputSingleScan(daq_sess,0);
    Screen('Flip', win);
    pause(dur_black);
    
    Screen('FillRect', win, WhiteIndex(win), rect);
    Screen('WaitBlanking',win);
    outputSingleScan(daq_sess,1);  
    Screen('Flip', win);
    pause(dur_white);
end