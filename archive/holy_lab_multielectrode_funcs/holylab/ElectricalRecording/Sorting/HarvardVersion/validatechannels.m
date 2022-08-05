function goodchannels = validatechannels(varargin)
% VALIDATECHANNELS: gui for selecting active channels

% Calling syntaxes:
%   validatechannels(ephys) sets up the figure and calls uiwait
%   validatechannels(string) handles the callbacks

if (nargin == 1 & ischar(varargin{1}))
  cmd = varargin{1};
  switch cmd
    case 'Done',
      uiresume
      return
    otherwise
      error('Unrecognized calling syntax');
  end
elseif (nargin == 1 & isstruct(varargin{1}))
  options.mingood = 10;
  options.anglethresh = 0.6;
  options.nplot = 50;
  options.maxchantoplot = 8;
  % This is the first call, do the automatic analysis
  ephysin = varargin{1};
  sniprange = ephysin(1).sniprange;
  allchannels = ephysin(1).channels;
  % Flip all spikes to positive polarity
  polsnips = cell(1,length(allchannels));
  for chani = 1:length(allchannels)
    tmpsnips = cell(1,length(ephysin));
    for j = 1:length(ephysin)
      tmpsnips{j} = ephysin(j).snippets{chani};
    end
    snips = [tmpsnips{:}];
    [sniplen,nsnips] = size(snips);
    if ~isempty(snips)
      pol = sign(snips(-sniprange(1)+1,:));
      polsnips{chani} = snips.*repmat(pol,sniplen,1);
    end
  end
  % Compute the mean snippet across channels (an
  % estimate of what spikes are supposed to look like)
  % Use only "exemplary channels" in computing the mean
  uiwait(warndlg('Pick some exemplary channels, to learn what spikes look like'));
  echans = pickchannels(polsnips);
  if isempty(echans)
    error('Need to pick some exemplary channels');
  end
  mnsnip = mean([polsnips{echans}],2);
  mnsnip = mnsnip/sqrt(sum(mnsnip.^2));
  figure
  plot(mnsnip)
  axis tight
  title('Mean snippet')
  % Now start drawing figures
  uiwait(warndlg('Now select active channels'));
  chanistart = 1;
  goodchannels = [];
  while (chanistart <= length(allchannels))
    ntodo = min([length(allchannels)-chanistart+1,options.maxchantoplot]);
    channels = allchannels(chanistart:chanistart+ntodo-1);
    nchan = length(channels);
    nrows = ceil(nchan/2);
    figure('Position',[63    48   800 600],'Visible','off',...
      'BackingStore','off','Units','pixels');
    % Set up the 'Done' button
    uicontrol('Style','pushbutton','Callback','validatechannels Done',...
      'String','Done','Position',[750 5 80 25]);
    hcheck = [];
    for i = 1:nchan
      chani = find(allchannels == channels(i));
      snips = polsnips{chani};
      [sniplen,nsnips] = size(snips);
      % Position for checkbox
      if (mod(i,2))
        subplot(nrows,6,3*i-2);
        set(gca,'Units','pixels');
        pos = get(gca,'Position');
        pos = [pos(1)-60 pos(2)+pos(4)/2];
      else
        subplot(nrows,6,3*i);
        set(gca,'Units','pixels');
        pos = get(gca,'Position');
        pos = [pos(1)+pos(3)+30 pos(2)+pos(4)/2];
      end
      if ~isempty(snips)
        snipproj = sum(snips.*repmat(mnsnip,1,nsnips));
        costh = snipproj ./ sqrt(sum(snips.^2));
        goodsnips = find(costh > options.anglethresh);
        badsnips = setdiff(1:nsnips,goodsnips);
        isgood = (length(goodsnips) > options.mingood & ...
          length(goodsnips) > length(badsnips));
        subplot(nrows,6,3*i-2);
        hist(costh);
        subplot(nrows,6,3*i-1)
        if ~isempty(goodsnips)
          plot(snips(:,goodsnips(1:min([options.nplot length(goodsnips)]))))
          axis tight
        else
          cla
        end
        subplot(nrows,6,3*i)
        if ~isempty(badsnips)
          plot(snips(:,badsnips(1:min([options.nplot length(badsnips)]))))
          axis tight
        else
          cla
        end
        % Put in the check box
        hcheck(i) = uicontrol('Style','checkbox','Max',1,'Min',0,'Value',isgood,...
          'Position',[pos 20 20]);
      else
        subplot(nrows,6,3*i-2)
        cla
        subplot(nrows,6,3*i-1)
        cla
        subplot(nrows,6,3*i)
        cla
        hcheck(i) = uicontrol('Style','checkbox','Max',1,'Min',0,'Value',0,...
          'Enable','off','Position',[pos 20 20]);
      end
    end
    set(gcf,'Visible','on')
    %return   % exit before calling uiwait, for debugging
    uiwait
    % When it gets here, the user has clicked 'Done'
    % Retrieve the data
    goodflag = get(hcheck,'Value');
    goodf = zeros(size(goodflag));
    for i = 1:length(goodflag)
      goodf(i) = goodflag{i};
    end
    good = find(goodf);
    goodchannels(end+1:end+length(good)) = channels(good);
    close(gcf)
    chanistart = chanistart+options.maxchantoplot;
  end
else
  error('Unknown calling syntax');
end
