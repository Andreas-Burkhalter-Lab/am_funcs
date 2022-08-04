function [varargout]= Do_average_inward(chan)
% input argument is recording channel (1 or 2); defaults to 1
%%% last edited AM 2/12/16
error('charge integration window not working - see AM to use this file')

charge_integration_window = 0.05; % seconds

files=list_files('./','*.xsg');
firstfile = deblank(files(1,:));
firstfile = load(firstfile,'-mat');
% npulses = length(firstfile.header.ephys.ephys.pulseParameters);
% charge_integral = inf(npulses,size(files,1));
charge_integral = inf(1,size(files,1));

multipulsewarn = 0;
for ii=1:size(files)
    IPSC=deblank(files(ii,:));
    d = load(IPSC, '-mat');
%     if length(d.header.ephys.ephys.pulseParameters)~=npulses
%         error('All traces must contain the same number of pulses.')
%     end
    if length(d.header.ephys.ephys.pulseParameters) > 1;
        multipulsewarn = 1;
    end
    
    if ~exist('chan','var') || chan==1
        tr = d.data.ephys.trace_1;
    elseif chan==2
        tr = d.data.ephys.trace_2;
    else
        error('Did not recognize channel name.')
    end
    
    % Get integral
    pars = d.header.ephys.ephys.pulseParameters{1,1};
    pulseStart = 1 + pars.sampleRate*pars.squarePulseTrainDelay;
    pulseEnd = pulseStart + pars.squarePulseTrainWidth*pars.sampleRate;
    pulseLength = pulseEnd-pulseStart;
    baseline_mean = mean(tr(1:pulseStart-1));
    charge_integral(1,ii) = sum(tr(pulseStart:pulseEnd)) - pulseLength*baseline_mean;
    
    %%% below not yet working for multiple pulses; don't know how to find
    %%% onset time of pulses after the first
%     for pulse = 1:npulses
%         pars = d.header.ephys.ephys.pulseParameters{1,pulse};
%         pulseStart = 1 + pars.sampleRate*(pars.squarePulseTrainDelay +... %start of first pulse
%             (pars.squarePulseTrainNumber-1)*pars.squarePulseTrainISI); % add offset for pulses after the first
%         pulseEnd = pulseStart + pars.squarePulseTrainWidth*pars.sampleRate; % maybe should use ISI instead of width?
%     end
        
    trace_data(:,ii)=detrend(tr);
    mdata=mean(trace_data,2);
   % plot(detrend(mdata(500:2000,:)));

end

trace_peaks = min(trace_data); % use min for downward deflections
output.mean_traces = mdata;
output.all_traces = trace_data;
output.mean_amplitude = min(mdata); % use min for downward deflections
peak_std = std(trace_peaks);
output.peak_std = peak_std;
output.largest_amplitude = min(trace_peaks); % use min for downward deflections
output.charge_integral = charge_integral;
output.mean_charge_integral = mean(charge_integral);

disp(['Charge integral window = ' num2str(charge_integration_window) 'sec.'])
if multipulsewarn
    warning('Found multiple pulses per trace; only calculating charge integral for the first pulse per trace.')
end

plot(mdata);

if nargout>=1  
    varargout{1} = output;
end

save avg_result output;