function [  ] = RD_avg_xsg( listname )
%RD_AVG_XSG Perform 60hz noise filtering on a specified list of .xsg trace files,
%then average all traces. 
%  Also display integral between current and baseline current for stim
%  duration.
%  Traces must all have equal duration and sample rate. 
%%%%%% Last updated 9/17/15 on recording comp

close all

%% Parameters
plot_mean = 1;              % plot the mean filtered traces
plot_separately = 1;        % plot each filtered trace separately
plot_picoAmp_seconds = 1;   % plot current integrals as a bar graph in units of picoAmp*seconds
plot_mean_picoAmps = 1;     % plot bar graph of mean current (in picoAmps) throughout the stimulation period
print_picoAmp_seconds = 1;  % display current integral values on the command line - units of picoAmp*seconds
print_mean_picoAmps = 1;    % display mean current values on the command line - units of picoAmp*seconds

width_override = 0.075;   %% if not empty, use this value as the time window in s for current integration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE: the two options below output values in arbitary units, NOT
%%% picoAmps or picoAmp*seconds. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_current_integral = 0;  % plot current integrals as a bar graph - NOTE: these values are not in meaningful units
print_current_integral = 0; % display current integral values on the command line - NOTE: these values are not in meaningful units



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
    filelist = RD_lists_traces2average(listname);
end

[junk filenames] = cellfun(@fileparts,(filelist(:,1)),'UniformOutPut',0);

for thisfile = 1:size(filelist,1)
    tracestruct(thisfile) = load(filelist{thisfile,1},'-mat');
    traces(:,thisfile) = eval(['tracestruct(thisfile).data.ephys.trace_' num2str(filelist{thisfile,2})]);
    samprate(thisfile) = tracestruct(thisfile).header.ephys.ephys.sampleRate;
end
if any(samprate(1)~=samprate)
    error('All sample rates must be equal.')
end

% Filter and average the traces.
% Resampling shouldn't be necessary because the sample rate (usually 10khz)
% is already high. 
oddRoundPeriod = 2*round((samprate(1)/filterfreq+1)/2) - 1; % use this period for Savitzky-Golay filtering
filttrace = sgolayfilt(traces,1,oddRoundPeriod); % filter all traces individually
meanfilttrace = mean(filttrace,2);      % take the mean of the filtered traces

% Get integral of current above baseline during stim duration.
for thisfile = 1:size(filelist,1)
        tfirst = tracestruct(thisfile).header.ephys.ephys.pulseParameters{1,filelist{thisfile,2}}.squarePulseTrainDelay;
    if isempty(width_override)
        tlast = tfirst + tracestruct(thisfile).header.ephys.ephys.pulseParameters{filelist{thisfile,2},1}.squarePulseTrainWidth;
    else
        tlast = tfirst + width_override;
    end
    stimwindow(thisfile,:) = round(samprate(thisfile)*[tfirst tlast] + [1 0]); % window start and stop in scans
    windowsize_seconds(thisfile) = [ stimwindow(thisfile,2) - stimwindow(thisfile,1) ] / samprate(thisfile);
    baseline(thisfile) = mean(filttrace(1:stimwindow(thisfile,1)-1,thisfile));
    
    % intgr_raw is in meaningless units; use meanCurrentInStimWindow  for
    % mean picoAmps current during stimulus, use picoAmpSecsInStimWindow
    % for current integral in picoAmp*seconds. 
    intgr_raw(thisfile) = sum(filttrace(stimwindow(thisfile,1):stimwindow(thisfile,2),thisfile) - baseline(thisfile));
end

picoAmpSecsInStimWindow = intgr_raw ./ samprate; % current integral during stim
meanCurrentInStimWindow = intgr_raw ./ (samprate.*windowsize_seconds); % average current over time during stim

%% Plotting
% Plot the averaged filtered trace.
if plot_mean
    figure
    timeaxis = (1:size(meanfilttrace,1))/samprate(1);   %% time axis with units in seconds
    plot(timeaxis, meanfilttrace)
    xlabel('Time (s)')
    ylabel('Voltage')
    title('Averaged Filtered Trace')
    xlabel('Time')
    ylabel('Current')
    set(gcf,'Position',[650 300 620 500])
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
    set(gcf,'Position',[1300 300 620 500])     
end

if plot_picoAmp_seconds
    figure
    bar(picoAmpSecsInStimWindow)
    ylabel('Current Integral (picoAmp*seconds)')
    set(gca,'XTickLabel',filenames)
    set(gcf,'Position',[0 300 620 500]);
end

if plot_mean_picoAmps
    figure
    bar(meanCurrentInStimWindow)
    ylabel('Mean Current During Stim Window (picoAmps)')
    set(gca,'XTickLabel',filenames)
    set(gcf,'Position',[0 300 620 500]);
end



%%% Plotting in meaningless units
if plot_current_integral
    figure
    bar(intgr_raw)
    ylabel('Current Integral (ARBITRARY UNITS)')
    set(gca,'XTickLabel',filenames)
    set(gcf,'Position',[0 300 620 500]);
end

%% Print values to command line.
if print_picoAmp_seconds
%     disp(sprintf('\nCurrent Integral (picoAmp*seconds)        Filename'))
    disp(sprintf('\nCurrent Integral (picoAmp*seconds)        '))
    for i = 1:size(filelist,1)
        disp([num2str(picoAmpSecsInStimWindow(i))])
    end
end

if print_mean_picoAmps
%     disp(sprintf('\nMean Current (picoAmps)       Filename'))
    disp(sprintf('\nMean Current (picoAmps)       '))
    for i = 1:size(filelist,1)
        disp([num2str(meanCurrentInStimWindow(i))])
    end
end


%%%% meaningless units
if print_current_integral
    disp(sprintf('\nCurrent Integral (arbitrary units)       Filename'))
    for i = 1:size(filelist,1)
     %   disp([num2str(intgr_raw(i)) '                ' filenames{i}])  %%
        disp([num2str(intgr_raw(i))])
    end
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