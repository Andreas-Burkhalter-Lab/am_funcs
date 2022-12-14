function [traceOut timeAxis] = rfMapSeeTrial(chans,iteration,row,column,event_bins,merecfile,...
    show_plot,bufferSecs,max_vsteps)
%%% Get and plot the .merec traces that occured during the trial described by
%%% iteration, row, and column during rf_mapping. Must have the event_bins
%%% dataset generated by rf_getevents (in rf_main) open
% chans = vector of channels to get traces from
% event_bins = variable of this name created by rf_getevents
% merecfile = file containing full raw trace of recording
% show_plot --- set to zero to not plot traces
% bufferSecs = time to also show before and after onset and offset
% max_vsteps = the same parameter from rf_getevents
%%% timeAxis = time axis for plotting traces

% Use defaults if unspecified.
if ~exist('bufferSecs') || isempty(bufferSecs)
    bufferSecs = 0;
end
if ~exist('max_vsteps') || isempty(max_vsteps)
    max_vsteps = 10;
end
if ~exist('chans') || isempty(chans)
    chan = input('Which channels to get traces from?');
end
if ~exist('merecfile') || isempty(merecfile)
    merecfile = uigetfile('Select .merec file.');
end
if ~exist('show_plot') || isempty(show_plot)
    show_plot = 1;
end

mObj = merecmm(merecfile);
rowcolumn = max_vsteps*(row) + column;
matchrow = event_bins.rowcolumn == rowcolumn &...
           event_bins.iteration == iteration;
bufferScans = bufferSecs*mObj.scanrate;
stimStartScan = event_bins.startScan(matchrow);
plotStartScan = stimStartScan-bufferScans;
stimEndScan = stimStartScan + event_bins.durScans(matchrow);
plotEndScan = stimEndScan+bufferScans;
traceOut = mObj(chans,plotStartScan:plotEndScan);
timeAxis = [plotStartScan:plotEndScan]./mObj.scanrate; % timeAxis units in seconds

if show_plot
    close all
    figure
    plot(timeAxis',traceOut')
    xlabel('Secs from beginning of recording')
    ylabel('Voltage')
    hold all
    stimStartSec = stimStartScan/mObj.scanrate;
    stimEndSec = stimEndScan/mObj.scanrate;
    scatter([stimStartSec stimEndSec],[0 0],'r'); % mark stim on an off
end
