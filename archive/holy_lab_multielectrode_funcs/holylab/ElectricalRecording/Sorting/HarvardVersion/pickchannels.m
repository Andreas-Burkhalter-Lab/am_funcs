function goodi = pickchannels(arg)
% PICKCHANNELS: choose a subset of channels, based on spike shape
% goodindx = pickchannels(snips)
% where
%   snips is a cell array of length nchans, each containing the matrix
%     of snippets on the given channel
% and
%   goodindx is the vector of selected channel indices
%
% Other functions (such as VALIDATECHANNELS) give more sophisticated
% information about spike shape.
% 
% See also: VALIDATECHANNELS, VALCHANRUN.
if (ischar(arg))
  % Callback processing
  hcheck = findobj(gcf,'Style','checkbox');
  val = get(hcheck,'Value');
  chani = get(hcheck,'UserData');
  gooditmp = get(gcf,'UserData');
  for i = 1:length(val)
    if (val{i} == 1)
      gooditmp(end+1) = chani{i};
    end
  end
  set(gcf,'UserData',gooditmp);
  switch arg
    case 'Done',
      uiresume
    case 'Quit',
      delete(hcheck); % sends a signal that we're quitting
      uiresume
    otherwise
      error('Unrecognized calling syntax');
  end
  return
else
  snips = arg;
  options.nplot = 50;
  options.nrows = 4;
  options.ncol = 4;
  nchans = length(snips);
  maxtoplot = options.nrows*options.ncol;
  chanistart = 1;
  goodi = [];
  while (chanistart <= nchans)
    chaniend = min([nchans,chanistart+maxtoplot-1]);
    channels = chanistart:chaniend;
    figure('Position',[63    48   800 600],'Visible','off',...
      'BackingStore','off','Units','pixels','UserData',goodi,...
      'CloseRequestFcn','pickchannels Quit');
    % Set up the 'Done' button
    hbutton = uicontrol('Style','pushbutton','Callback','pickchannels Done',...
      'String','Done','Position',[750 5 80 25]);
    if (chaniend < nchans)
      set(hbutton,'String','Next');
    end
    hcheck = [];
    for i = 1:length(channels)
      subplot(options.nrows,options.ncol,i);
      % Adjust size of axis to make room for checkbox
      set(gca,'Units','pixels');
      pos = get(gca,'Position');
      set(gca,'Position',pos - [0 0 20 0]);
      set(gca,'Units','normalized')
      pos = [pos(1)+pos(3)-10 pos(2)+pos(4)/2];
      [sniplen,nsnips] = size(snips{channels(i)});
      if nsnips > 0
        plot(snips{channels(i)}(:,1:min([options.nplot nsnips])));
        axis tight
        hcheck(i) = uicontrol('Style','checkbox','Max',1,'Min',0,'Value',0,...
          'Position',[pos 20 20],'UserData',channels(i));
      else
        cla
        hcheck(i) = uicontrol('Style','checkbox','Max',1,'Min',0,'Value',0,...
          'Enable','off','Position',[pos 20 20],'UserData',channels(i));
      end
    end
    set(gcf,'Visible','on')
%    return   % exit before calling uiwait, for debugging
    uiwait

    % When it gets here, the user has clicked 'Next/Done' or quit the figure

    % Retrieve the data
    goodi = get(gcf,'UserData');
    if any(~ishandle(hcheck))
      delete(gcf)
      return; % User closed the figure
    end
    delete(gcf)
    chanistart = chanistart + maxtoplot;
  end
end
