function [] = mapopen( file, channels, plottype )
%MAPOPEN Plot an .xsg file acquired during mapping
%   Last updated by AM 5/4/15.
% channels = array of channel numbers to plot; defaults to all recorded channels
% plottype = 'heatmap' (default) or 'trace'
%%% need to figure out where correct time markers are; not clear what data.ephys.dataEventTimestamps_1 is 
if ~exist('file','var')
    'Which mapper file to plot?'
    file = uigetfile('-.xsg');
end
load(file, '-mat');
gridorder = header.mapper.mapper.mapPatternArray;
snip_secs = header.mapper.mapper.version;
samprate = header.mapper.mapper.sampleRate;
snip_scans = round(snip_secs*samprate);

if ~exist('channels','var')
    channels = 1:length(find(header.ephys.ephys.acqOnArray));
end

for chan = 1:channels
    gridout{chan} = -inf(size(gridorder));
    wave = eval(['data.ephys.trace_' num2str(channels(chan))]); % long trace
    tstamps = eval(['data.ephys.dataEventTimestamps_' num2str(channels(chan))]);
    for gridspace = 1:numel(gridorder)
        traceind = gridorder(gridspace); % which snippet to take from the long trace
        gridout{chan}(gridspace) = wave(round(tstamps(traceind)) : tstamps(traceind)+snip_scans); % get the snippet
    end
    
    % Plotting
    if ~strcmp(plottype, 'trace') % default to heat map
        baseline = median(wave); % could use more sophisticated way to find baseline
        responses = cell2mat(cellfun(@sum,gridout{chan})) - baseline*snip_scans;
      
        figure
        imagesc(responses);
        
    else
        for row = 1:size(gridout,1)
            for col = 1:size(gridout,2)
                subplot(size(gridout{chan},1), size(gridout{chan},2), col + (row-1)*size(gridout{chan},2))
                plot(gridout{chan}{row,column})
            end
        end
    end
        
    title(eval(['data.ephys.amplifierName_' num2str(chan)]))
    

% clearvars data header

end

