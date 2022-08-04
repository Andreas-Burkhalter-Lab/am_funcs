if ~exist('s','var')
    daq_sess = daq.createSession('ni');
    addAnalogOutputChannel(daq_sess,'Dev1','ao0','Voltage'); 
    daq_sess.IsContinuous = 1;
end
    

stop(daq_sess);
%% Analog Outputs
val = 05;
% q = val * ones(1000,1);
q = [-1*ones(500,1); 0*ones(500,1); 1*ones(500,1); 2*ones(500,1);3*ones(500,1)];  
% q = [-val; val* ones(500,1);*ones(49,1) ]; 
% q = repmat([val; -val], 500, 1);

% outputSingleScan(s,0);            %%% set channel to zero V
% q = (0:0.0001:val)';
% q = zeros(5000,1);
% q = [(-val:0.001:val)']; %val*ones(1000,1)];
% q = [-val:val/5000:val]';
daq_sess.queueOutputData(q); 
daq_sess.startBackground;
% s.startForeground;



%% Digital Outputs
% addDigitalChannel(s,'Dev1','port1/line1','OutputOnly');
% outputSingleScan(s,0)
% outputSingleScan(s,1)