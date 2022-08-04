%%%%%%%%%%% send test pulse
%%%%% upadted 18-10-18
%% params

device_name = 'Dev2';
output_chan = 'ao0';
val = 5;
dur_seconds = 10;
q = [val * ones(1000*dur_seconds,1); 0];
% q = [-1*ones(500,1); 0*ones(500,1); 1*ones(500,1); 2*ones(500,1);3*ones(500,1)];  
% q = [-val; val* ones(500,1);*ones(49,1) ]; 
% q = repmat([val; -val], 500, 1);


%% Analog Outputs
global daq_sess
clear global daq_sess
daq_sess = daq.createSession('ni');
addAnalogOutputChannel(daq_sess,device_name,output_chan,'Voltage'); 
daq_sess.IsContinuous = 0;

    
stop(daq_sess);


% outputSingleScan(s,0);            %%% set channel to zero V
% q = (0:0.0001:val)';
% q = zeros(5000,1);
% q = [(-val:0.001:val)']; %val*ones(1000,1)];
% q = [-val:val/5000:val]';
daq_sess.queueOutputData(q); 
% daq_sess.startBackground;
daq_sess.startForeground;



%% Digital Outputs
% addDigitalChannel(s,'Dev1','port1/line1','OutputOnly');
% outputSingleScan(s,0)
% outputSingleScan(s,1)