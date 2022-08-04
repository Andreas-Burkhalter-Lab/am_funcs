function [ raw_traces charge_raw  ] = PB_avg_xsg( listname )
%PB_AVG_XSG Perform 60hz noise filtering on a specified list of .xsg trace files,
%then average all traces. 
%  Also display integral between current and baseline current for stim
%  duration.
%  Traces must all have equal duration and sample rate. 
%%%%%% Last updated 6/20/15 ephus comp

%% Parameters
plot_mean.show = 0;                          % plot the mean filtered traces
    plot_mean.pickXAxis = 0;                   % pick the limits of the x and y axes by clicking on the plot 
plot_mean_over_individual.show = 1;          % overlay mean trace in bold over individual traces as thin lines
    plot_mean_over_individual.pickColor = 1; % choose colors for the mean and individual traces
    plot_mean_over_individual.pickXAxis = 1; % pick the limits of the x and y axes by clicking on the plot
plot_separately = 0;                         % plot each filtered trace separately
plot_current_integral = 0;                   % plot current integrals as a bar graph
print_current_integral = 0;                  % display current integral values on the command line

width_override = 0.050;   %% if not empty, use this value as the time window in s for current integration

%% Analysis
filterfreq = 60; % frequency to filter out in hz
if ~exist('listname','var')
    nfiles = input('How many files to filter and average? ');
    fprintf('Files being averaged:\n');
    for fileind = 1:nfiles
        filelist{fileind, 1} = uigetfile('*.xsg');
        chancheck = questdlg(sprintf('File '' %s \nWhich channel to take trace from?', filelist{fileind,1}),...
            'Pick Trace', 'Channel 1', 'Channel 2', 'Channel 1');
        filelist{fileind, 2} = str2num(chancheck(end));
        fprintf('%s [%s]\n',filelist{fileind, 1}, chancheck);
    end
else
    filelist = PB_lists_traces2average(listname);
end

[junk filenames] = cellfun(@fileparts,(filelist(:,1)),'UniformOutPut',0);

for thisfile = 1:size(filelist,1)
    tracestruct(thisfile) = load(filelist{thisfile,1},'-mat');
    raw_traces(:,thisfile) = eval(['tracestruct(thisfile).data.ephys.trace_' num2str(filelist{thisfile,2})]);
    samprate(thisfile) = tracestruct(thisfile).header.ephys.ephys.sampleRate;
end
if any(samprate(1)~=samprate)
    error('All sample rates must be equal.')
end

% Filter and average the traces.
% Resampling shouldn't be necessary because the sample rate (usually 10khz)
% is already high. 
oddRoundPeriod = 2*round((samprate(1)/filterfreq+1)/2) - 1; % use this period for Savitzky-Golay filtering
filttrace = sgolayfilt(raw_traces,1,oddRoundPeriod); % filter all traces individually
meanfilttrace = mean(filttrace,2);      % take the mean of the filtered traces

% Get integral of current above baseline during stim duration.
for thisfile = 1:size(filelist,1)
        tfirst = tracestruct(thisfile).header.ephys.ephys.pulseParameters{filelist{thisfile,2}}.squarePulseTrainDelay;
    if isempty(width_override)
        tlast = tfirst + tracestruct(thisfile).header.ephys.ephys.pulseParameters{filelist{thisfile,2}}.squarePulseTrainWidth;
    else
        tlast = tfirst + width_override;
    end
    stimwindow(thisfile,:) = round(samprate(thisfile)*[tfirst tlast] + [1 0]); % window start and stop in scans
    baseline_raw(thisfile) = mean(filttrace(1:stimwindow(thisfile,1)-1,thisfile));
    baseline(thisfile) = mean(filttrace(1:stimwindow(thisfile,1)-1,thisfile));
    charge_raw(thisfile) = sum(raw_traces(stimwindow(thisfile,1):stimwindow(thisfile,2),thisfile) - baseline_raw(thisfile));
    intgr(thisfile) = sum(filttrace(stimwindow(thisfile,1):stimwindow(thisfile,2),thisfile) - baseline(thisfile));
end

%% Plotting
% Plot the averaged filtered trace.
timeaxis = (1:size(meanfilttrace,1))/samprate(1);   %% time axis with units in seconds
if plot_mean.show
    figure
    plot(timeaxis, meanfilttrace)
    xlabel('Time (s)')
    ylabel('Voltage')
    title('Averaged Filtered Trace')
    ylabel('Current')
    if plot_mean.pickXAxis
        pick_x_axis
    end
end

if plot_mean_over_individual.show
    if plot_mean_over_individual.pickColor
        colMean = uisetcolor('Color for mean trace');
        colIndividual = uisetcolor('Color for individual traces');
    else
        colMean = [1 1 1];
        colIndividual = [0.7 0.7 0.7];
    end
    figure; hold on
    plot(timeaxis',filttrace,'Color',colIndividual,'LineWidth',1)
    plot(timeaxis',meanfilttrace,'Color',colMean,'LineWidth',2)
    xlabel('Time (ms)')
    ylabel('Current')
    legend('Individual Filtered Trials','Mean Filtered Trials')
    if plot_mean_over_individual.pickXAxis
        pick_x_axis
    end
end

if plot_separately
    figure
    rows = ceil(size(filelist,1)/4);
    for thisfile = 1:(size(filelist,1))
        subplot(rows,4,thisfile)
        plot(filttrace(:,thisfile))
        title(filenames{thisfile},'Interpreter','none')
        xlabel('Time')
        ylabel('Current')
    end
end

if plot_current_integral
    figure
    bar(intgr)
    ylabel('Current Integral')
    set(gca,'XTickLabel',filenames)
end

if print_current_integral
    disp(sprintf('\nCurrent Integral        Filename'))
    for i = 1:size(filelist,1)
        disp([num2str(intgr(i)) '                ' filenames{i}])
    end
end

end

%%
function [] = pick_x_axis()
    disp('Pick X-axis lower bound.')
    [xlow ~] = ginput(1);
    disp('Pick X-axis upper bound.')
    [xhigh ~] = ginput(1);
    xlim([xlow xhigh]);
end


% % % % % %% Use a Butterworth filter; not as effective as Savitzky-Golay
% % % % % filt60hz = designfilt('bandstopiir','FilterOrder',2, ...
% % % % %                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
% % % % %                'DesignMethod','butter','SampleRate',samprate(1));
% % % % %                   
% % % % % % Filter the traces.
% % % % % for thisfile = 1:length(filelist)
% % % % %     butterLoop(:,thisfile) = filtfilt(filt60hz, traces(:,thisfile)); 
% % % % % end
% % % % % avgtrace = mean(butterLoop,2);