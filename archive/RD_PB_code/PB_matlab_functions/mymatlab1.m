function mymatlab1(chan)
% input argument is recording channel (1 or 2); defaults to 1
%https://openwiki.janelia.org/wiki/plugins/viewsource/viewpagesrc.action?pageId=2130030
% d = load('AA0003AAAA0159.xsg', '-mat');
%%% last edited AM 7/29/15
filename = uigetfile('*.xsg','Select .xsg file.');
d = load(filename,'-mat');
    if ~exist('chan','var') || chan==1
        tr = d.data.ephys.trace_1;
    elseif chan==2
        tr = d.data.ephys.trace_2;
    else
        error('Did not recognize channel name.')
    end
%sr = d.header.ephys.ephys.sampleRate;
plot(tr);