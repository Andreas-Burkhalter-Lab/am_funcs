function varargout = RD_spike_latency( listname )
%Analyze, in paired recordings, the relationship between presynaptic input 
% vs. postsynaptic first-spike latency. Optionally provide an output
% argument, which will be the full results dataset.
%%% Presynaptic vs. postsynaptic channels can be defined in
%%% RD_lists_spikeLatency.m. 
%%%%%% Last updated 10/12/15 on recording comp

%% User Parameters
% First-spike latency will be found at the first peak following the first
% threshold crossing in the postsynaptic trace. 
spikeThreshold = 0; % threshold for finding the first spike (in same units as postsynaptic traces)
checkspikes = 1; % if on, display plot of all traces to make sure that first-spike latency was detected correctly
    checkspikes_rows  = 4; % if checkspikes = 1, this is the maximum number of rows of traces to display (all traces will be displayed)
print_latencies = 1; % if on, output latencies in seconds to the command line
save_results = 0; % if on, save results in a file with name listed on next line
    save_filename = 'RD_spike_latency_analyzed.mat'; % if save_results = 1, this is the filename for saving the analysis results

%% Load traces 
close all
if ~exist('listname','var') || isempty(listname)
    nfiles = input('How many files to analyze for spike latency vs. input strength? ');
    fprintf('Files being averaged:\n');
    for fileind = 1:nfiles
        filelist{fileind, 1} = uigetfile('*.xsg');
        chancheck = questdlg(sprintf('File '' %s \n... Which channel contains the postsynaptic trace?',...
            filelist{fileind,1}),'Pick Trace', 'Channel 1', 'Channel 2', 'Channel 1');
        filelist{fileind, 2} = str2num(chancheck(end));
        fprintf('%s [%s]\n',filelist{fileind, 1}, chancheck);
       
    end
else
    filelist = RD_lists_spikeLatency(listname);
    nfiles = size(filelist,1);
end

[junk filenames] = cellfun(@fileparts,(filelist(:,1)),'UniformOutPut',0);

dsNans = NaN(nfiles,1);
dsCell = cell(nfiles,1);
latency_results = dataset(dsCell,dsNans,dsCell,dsCell,dsCell,'VarNames',...
    {'pulse_name','spike_latency_sec','filename','presyn_trace','postsyn_trace'});
%     {'pulse_name','presyn_input_amp','spike_latency_sec','filename','presyn_trace','postsyn_trace'});

for traceInd = 1:length(filelist)
    load(filelist{traceInd,1},'-mat');
    chanPostsyn = filelist{traceInd,2};
    if chanPostsyn == 1
        chanPresyn = 2;
    elseif chanPostsyn ==2
        chanPresyn = 1;
    else
        error('Postsynaptic channel must be either 1 or 2.')
    end
    tracePostsyn = eval(['data.ephys.trace_' num2str(chanPostsyn)]);
    tracePresyn = eval(['data.ephys.trace_' num2str(chanPresyn)]);
    samprate = header.ephys.ephys.sampleRate;
    
    % Find first-spike latency.
    if ~isempty(find(tracePostsyn>=spikeThreshold)); % if postsyn trace has a threshold crossing
        threshCrossIndex = find((tracePostsyn>=spikeThreshold),1);
        firstPeakIndex = threshCrossIndex - 1 +...                               %% add threshCross time to...
                         find( diff(tracePostsyn(threshCrossIndex:end))<=0, 1); %% ... latency between threshCross and first non-upward slope
        stimOnsetTime = header.ephys.ephys.pulseParameters{1,chanPresyn}.squarePulseTrainDelay; % stim latency from recording start in seconds
        stimOnsetIndex = round(stimOnsetTime*samprate); % convert to scan index to compare with spike timing
        firstSpikeLatency = (firstPeakIndex - stimOnsetIndex) / samprate;
    else % if postsyn trace does not have a threshold crossing
        firstSpikeLatency = NaN;
    end
%     presyn_input_amp = header.ephys.ephys.pulseParameters{1,chanPresyn}.amplitude; % not sure what this value is - does not correspond to picoAmps
             
    %% Put results into a dataset.
    %%% NOTE: picoAmps from pulse_name does NOT match amplitude (units unknown) 
    %%% listed in header.ephys.ephys.pulseParameters{1,chanPresyn}.amplitude     
    latency_results.pulse_name{traceInd} = header.ephys.ephys.pulseParameters{1,chanPresyn}.name;
%     latency_results.presyn_input_amp(traceInd) = presyn_input_amp;
    latency_results.spike_latency_sec(traceInd) = firstSpikeLatency;
    latency_results.filename{traceInd} = filelist{traceInd,1};
    latency_results.presyn_trace{traceInd} = tracePresyn;
    latency_results.postsyn_trace{traceInd} = tracePostsyn;
    
    %% Display first-spike timing for visual checking.
    if checkspikes
        ncolumns = ceil(nfiles/checkspikes_rows);
        subplot(checkspikes_rows,ncolumns,traceInd)
        plotwidth = min(1500,ncolumns*400);
        set(gcf,'Position',[100 100 plotwidth 900]);
        hold on
        timeAxis = linspace(0,length(tracePostsyn)/samprate,length(tracePostsyn));
        plot(timeAxis,tracePostsyn);
        if isnan(firstSpikeLatency) % if no threshold crossing for this trace
            legend('No spikes found')
        else % if there was a spike in this trace
            firstPeakVoltage = tracePostsyn(firstPeakIndex);
            scatter(firstPeakIndex/samprate,firstPeakVoltage,'r');
%             legend('First spike')
        end
            
        ylabel('Voltage')
        xlabel('Time (seconds)')
%         xlabel('Time (scans)')
        title(filenames{traceInd});
        suptitle('Postsynaptic traces with first spikes marked')
        
    end
end

if nargout > 0
    varargout{1} = latency_results;
end

if print_latencies
    fprintf('\n Postsynaptic spike latencies (seconds):\n')
    disp(latency_results.spike_latency_sec)
end
    
% Save results.
if save_results
    save(save_filename,'latency_results','filelist','nfiles','spikeThreshold','listname');
end

    