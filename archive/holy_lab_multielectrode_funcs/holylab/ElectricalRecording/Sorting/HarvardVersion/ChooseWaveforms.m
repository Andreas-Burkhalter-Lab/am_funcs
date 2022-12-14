function fig = ChooseWaveforms(spikes,sniprange,thresh)
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

if (nargin < 3)
        peakval = sort(abs(spikes(-sniprange(1)+1,:)));
        thresh = peakval(ceil(length(peakval)/2));
end
% emode: controls the erase mode for the waveforms
% The two real choices are: 'xor' and 'normal'
% normal looks better, but is slower on operations that
% delete only a few waveforms
emode = 'xor';
nspikes = size(spikes,2);
h0 = figure('Units','points', ...
        'Color',[0.8 0.8 0.8], ...
        'Position',[294 276 604 291], ...
        'Tag','ChooseWaveformsWindow',...
        'KeyPressFcn','ChooseWfmsCallback Delete');
hwf = axes('Parent',h0, ...
        'Units','pixels', ...
        'CameraUpVector',[0 1 0], ...
        'Color',[1 1 1], ...
        'Position',[43 39 328 233], ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'ZColor',[0 0 0]);
hthresh = axes('Parent',h0, ...
        'CameraUpVector',[0 1 0], ...
        'Color',[1 1 1], ...
        'Position',[0.6837748344370861 0.4604810996563574 0.2764900662251656 0.4707903780068728], ...
        'XColor',[0 0 0], ...
        'YColor',[0 0 0], ...
        'ZColor',[0 0 0]);
h1 = uicontrol('Parent',h0, ...
        'Units','points', ...
        'Position',[413 38 69 30], ...
        'String','Cancel', ...
        'Tag','CancelButton', ...
        'Callback','ChooseWfmsCallback Cancel');
h1 = uicontrol('Parent',h0, ...
        'Units','points', ...
        'Position',[511 38 69 30], ...
        'String','Done', ...
        'Tag','DoneButton', ...
        'Callback','ChooseWfmsCallback Done');
htext = uicontrol('Parent',h0, ...
        'Units','points', ...
        'Position',[412 81 166 16], ...
        'Style','text', ...
        'Tag','NumWaveforms');
if nargout > 0, fig = h0; end
% Histogram the peak values for threshold setting
axes(hthresh)
offset = -sniprange(1)+1;
[n,x] = hist(abs(spikes(offset,:)),round(nspikes^(1/3)));
bar(x,n+1);
set(gca,'YScale','log');
yrange = get(gca,'YLim');
% Plot the threshold-setting line
hline = line([thresh,thresh],yrange,...
        'LineStyle','-',...
        'EraseMode','xor',...
        'ButtonDownFcn','ChooseWfmsCallback ThreshStart');
xlabel('Peak value');
ylabel('#/bin');
set(hthresh,'Tag','ThreshAxes');
% Now plot all the spike waveforms
axes(hwf)
TAxis = sniprange(1):sniprange(2);
LineH = plot(TAxis,spikes,'SelectionHighlight','off','EraseMode',emode,'Tag','wfm',...
        'ButtonDownFcn','ChooseWfmsCallback SelectLine','Visible','off');
set(gca,'XLim',sniprange);
% Only show the ones above threshold
onindx = find(abs(spikes(-sniprange(1)+1,:)) >= thresh);
set(LineH(onindx),'Visible','on');
% Plot width-selection lines
yrange = get(gca,'YLim');
line([sniprange(1)+1 sniprange(1)+1],yrange,'Color','k','LineStyle',':','EraseMode','xor',...
        'ButtonDownFcn','ChooseWfmsCallback WidthStart','Tag','WidthLine');
line([sniprange(2)-1 sniprange(2)-1],yrange,'Color','k','LineStyle',':','EraseMode','xor',...
        'ButtonDownFcn','ChooseWfmsCallback WidthStart','Tag','WidthLine');
% Give the waveforms an index#, in case we want to refer to them that way
indxnv = 1:length(LineH);
indxn = num2cell(indxnv);
%for i = 1:length(LineH)
%        indxn{i} = i;
%end
set(LineH,{'UserData'},indxn');
% Update the text
nwfm = length(onindx);
set(htext,'String',sprintf('%d waveforms',nwfm));
% Set callback properties
set(hwf,'Tag','WaveformsAxes', ...
        'ButtonDownFcn','ChooseWfmsCallback SelectRegion');
%set(LineH,'SelectionHighlight','off','EraseMode',emode,'Tag','wfm',...
%        'ButtonDownFcn','ChooseWfmsCallback SelectLine');
% Remember some useful #s
%setappdata(hwf,'PeakPos',-sniprange(1)+1);
%setappdata(hwf,'oldthresh',thresh);
setappdata(h0,'PeakPos',-sniprange(1)+1);
setappdata(h0,'oldthresh',thresh);
