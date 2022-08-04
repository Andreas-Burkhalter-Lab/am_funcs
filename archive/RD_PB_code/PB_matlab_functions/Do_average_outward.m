function [varargout]= Do_average_outward(chan)
% input argument is recording channel (1 or 2); defaults to 1
%%% last edited AM 2/12/16

%%%% test on: C:\Ephus\R210\Ephus\Data\PawanData\Second
%%%% Project\PV_cre\2016\Feb 2016\020816\PB0001\flashes\2.5 HZ....... test
%%%% on 3 or 4 files in this folder

charge_integration_window_override = 1; % set to 1 to use custom window
    charge_integration_window = 0.4; % seconds... only has effect if charge_integration_window_override==1

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
    
    %% following line is only to correct for ephus getting the wrong train
    %% start time
% % % % % %     pars.squarePulseTrainDelay = 0.1;
    
    pulseStart = 1 + pars.sampleRate*pars.squarePulseTrainDelay;
    if charge_integration_window_override
        pulseEnd = pulseStart + charge_integration_window*pars.sampleRate;
    else
        pulseEnd = pulseStart + pars.squarePulseTrainWidth*pars.sampleRate;
    end
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
        
    trace_data(:,ii)=detrend(tr); % remove linear trend, get trace relative to its baseline
    mdata=mean(trace_data,2); % get mean trace
   % plot(detrend(mdata(500:2000,:)));

end

trace_peaks = max(trace_data); % use max for upward deflections
output.mean_traces = mdata;
output.all_traces = trace_data;
output.mean_amplitude = max(mdata); % use max for upward deflections
peak_std = std(trace_peaks);
output.peak_std = peak_std;
output.largest_amplitude = max(trace_peaks); % use max for upward deflections
output.charge_integral = charge_integral;
output.mean_charge_integral = mean(charge_integral);


if charge_integration_window_override
    windowstring = '(user-specified value).';
else
    windowstring = '(stimulation length from ephus data).';
end
disp(['Charge integral window = ' num2str(charge_integration_window) 'sec' windowstring '.'])
if multipulsewarn
    warning('Found multiple pulses per trace; only calculating charge integral for the first pulse per trace.')
end

plot(mdata);

if nargout>=1  
    varargout{1} = output;
end

save avg_result output;