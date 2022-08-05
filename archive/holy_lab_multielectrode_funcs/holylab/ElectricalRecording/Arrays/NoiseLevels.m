function NoiseLevels(options)
% NOISELEVELS: Compute rms noise & power spectra for multichannel records
% NoiseLevels by itselft runs the default program
% noiselevels(options) allows default behavior to be changed, by setting
% fields in the structure options:
%   toffset: skips forward toffset seconds into the file (default 0)
%   duration: number of seconds of data to load (default 2)
%   meansubtract: if true, the mean signal across all channels is subtracted
%       before computed the rms (default false)
%   showchannel: user-controlled channel for plotting waveform & power
%       spectrum.
%
% An alternative syntax: noiselevels(toffset)
[file,path] = uigetfile('*.merec','Select raw data file');
tmax = 2;                % Time in seconds to load
toffset = 0;
meansubtract = 0;
if (nargin > 0)
    if isnumeric(options)
        toffset = options;
    elseif isstruct(options)
        if isfield(options,'toffset')
            toffset = options.toffset;
        end
        if isfield(options,'duration')
            tmax = options.duration;
        end
        if isfield(options,'meansubtract')
            meansubtract = options.meansubtract;
        end
    else
        error('Unrecognized calling syntax')
    end
end
[d,h] = loadmc([path,file],[toffset,toffset+tmax],struct('tovolts',1));
if meansubtract
    d = d - repmat(mean(d),size(d,1),1);
end
n = std(d');
figure
subplot(2,1,1)
plot(h.channels,n,'*')
axis tight
xlabel('Channel #');
ylabel('RMS signal (V)');
title(h.usrhdr);
subplot(2,1,2);
nsort = sort(n);
frac = (1:length(nsort))/length(nsort);
plot(nsort,frac)
xlabel('RMS signal (V)');
ylabel('Fraction');
% Power spectrum of the channel with median rms noise
[sn,indx] = sort(n);
repchan = indx(round(length(n)/2));
if (nargin > 0)
    if isfield(options,'showchannel')
        repchan = find(h.channels == options.showchannel);
        if (length(repchan) ~= 1)
            error('Desired channel was not recorded');
        end
    end
end
figure
subplot(2,1,1)
t = (0:size(d,2)-1)/h.scanrate;
plot(t,d(repchan,:));
axis tight
ylabel('Voltage (V)');
xlabel('Time (s)');
title(sprintf('Signal on channel %d',h.channels(repchan)));
subplot(2,1,2);

psd(d(repchan,:),16384,h.scanrate);
set(gca,'XScale','log');
xlabel('Frequency (Hz)');
% Find out if user want to do the correlation calculation
% (takes some time)
ButtonName = questdlg('Do you want to calculate cross-channel correlations?');
switch ButtonName
        case 'Yes'
        % Plot correlations between channels
        figure
        [ac,cc] = allcorranalog(d);
        imagesc(cc,[0 1])
        axis square
        colormap(1-gray);
        set(gca,'YDir','normal');
        title('Cross-correlations between channels');
        xlabel('Channel INDEX #');
        ylabel('Channel INDEX #');
end %switch
