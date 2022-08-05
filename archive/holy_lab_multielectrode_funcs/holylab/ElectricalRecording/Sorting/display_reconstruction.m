function display_reconstruction(varargin)
% DISPLAY_RECONSTRUCTION: show waveform reconstructions from snippets
%
% This function creates the reconstruction in a somewhat sophisticated way:
% at low resolution (large time scales), each snippet is represented simply
% by a vertical line from its minimum to its maximum value. At high
% resolution, each snippet is drawn. The advantages of this more complex
% approach are both in drawing speed and in memory consumption.
%
% This function is also capable of coloring each snippet according to its
% cluster assignment.
%
%
% Syntaxes:
%
%  display_reconstruction(shc,snipminmax)
%    The simplest call. shc is a channel-specific sortheader, of the type
%    you obtain with SNIPFILE2SORTHEADER followed by SORTHEADER_IMPORTCHAN.
%    snipminmax is a cell array of snippet min/max pairs, of the type you
%    obtain with SORTHEADER_SNIPMINMAX.
%    Important: this syntax requires that you supplying min/max values for
%    every snippet on the chosen channel. See below for alternatives.
%
%   display_reconstruction(shc,snipminmax,clusternum)
%     This syntax allows you to colorize snippets according to their
%     clusternum.  clusternum{i} is a vector of integer labels (starting at
%     0), one for each snippet.
%
%   display_reconstruction(shc,snipminmax,clusternum,snipIndex)
%     Use this syntax if you have min/max data on only a subset of spikes.
%     snipIndex is a cell array, with snipIndex{i} containing the list of
%     snippets described in file shc(i) by snipminmax{i} and clusternum{i}.
%     Note that clusternum can be given as an empty matrix or cell array,
%     in which case each snippet is assigned to cluster 0.
%
%   display_reconstruction(hax,...)
%     Use this syntax to target the output to a particular axis with handle
%     hax.
%
%   display_reconstruction(hax,'XLim',xlim,'YLim',ylim)
%     This implements the "zoom" functionality.  The x-limits are examined
%     to determine whether the resolution justifies displaying each snippet
%     waveform. If not, min/max lines are drawn.
%     Note that when snippet waveforms are drawn, this function will
%     display _every_ snippet waveform that occurs within the selected
%     timerange, even if omitted in snipIndex. Snippets omitted in
%     snipIndex are drawn in gray with a dashed line.
%
% This function plays well with SLIDERWINDOW if you set
%    options.axisupdatefcn = @display_reconstruction
% in your call to SLIDERWINDOW.
%
% See also: SORTHEADER_SNIPMINMAX, SORTHEADER_IMPORTCHAN, SLIDERWINDOW.

% Copyright 2006 by Timothy E. Holy

  new_window = false;
  clusternum = {};
  snipIndex = {};
  if isstruct(varargin{1})
    new_window = true;
    hax = gca;
    shc = varargin{1};
    snipminmax = varargin{2};
    if (nargin > 2)
      clusternum = varargin{3};
    end
    if (nargin > 3)
      snipIndex = varargin{4};
    end
  elseif (ishandle(varargin{1}) && isstruct(varargin{2}))
    new_window = true;
    hax = varargin{1};
    shc = varargin{2};
    snipminmax = varargin{3};
    if (nargin > 3)
      clusternum = varargin{4};
    end
    if (nargin > 4)
      snipIndex = varargin{5};
    end
  else
    if (~ishandle(varargin{1}) || ~ischar(varargin{2}))
      error('Input syntax not recognized');
    end
    hax = varargin{1};
  end
  % Defer any additional argument parsing for later

  if new_window
    %
    %% We're starting a new graph
    %
    cla(hax,'reset')  % clear any existing objects & properties
    setappdata(hax,'shc',shc);
    n_files = length(shc);
    % Supply any default values and put into cell array format
    if ~iscell(snipminmax)
      snipminmax = {snipminmax};
    end
    if isempty(snipIndex)
      for fileIndex = 1:n_files
        n_snips = length(shc(fileIndex).sniptimes);
        if (size(snipminmax{fileIndex},2) ~= n_snips)
          error('When snipIndex is empty, you must supply data on every snippet');
        end
        snipIndex{fileIndex} = 1:n_snips;
      end
    end
    if ~iscell(snipIndex)
      snipIndex = {snipIndex};
    end
    setappdata(hax,'snipIndex',snipIndex);
    if isempty(clusternum)
      for fileIndex = 1:n_files
        clusternum{fileIndex} = zeros(1,length(snipIndex{fileIndex}));
      end
    end
    if ~iscell(clusternum)
      clusternum = {clusternum};
    end
    % Calculate the maximum cluster number
    max_clusternum = -Inf;
    for fileIndex = 1:n_files
      if ~isempty(clusternum{fileIndex})
        max_clusternum = max(max_clusternum,max(clusternum{fileIndex}));
      end
    end
    setappdata(hax,'clusternum',clusternum);
    setappdata(hax,'max_clusternum',max_clusternum);
    % Calculate many things dealing with timing
    % Get the file start times
    tstart = sortheader_absolute_starttime(shc);
    % Get the snippet times
    t = sortheader_absolute_sniptime(shc,snipIndex);
    % Calculate the relative times
    for fileIndex = 1:n_files
      t{fileIndex} = t{fileIndex} - min(tstart);
    end
    tstart = tstart - min(tstart);
    setappdata(hax,'file_starttime',tstart);
    % Calculate file duration
    tlen = [shc.nscans] ./ [shc.scanrate];
    setappdata(hax,'file_endtime',tstart+tlen);
    % Set the original x limits to the duration of the experiment
    trange = [0 max(tstart+tlen)];
    trange = trange + [-1 1]*diff(trange)/100;  % make 1% wider
    setappdata(hax,'XLim0',trange)
    % Set y limits
    ymin = Inf;
    ymax = -Inf;
    for fileIndex = 1:n_files
      ymin = min(ymin,min(snipminmax{fileIndex}(:)));
      ymax = max(ymax,max(snipminmax{fileIndex}(:)));
    end
    yrange = [ymin ymax];
    yrange = yrange + [-1 1]*diff(yrange)/20; % make 5% wider
    setappdata(hax,'YLim0',yrange);
    % Draw the min/max lines (we do this even if the initial display will
    % be "zoomed in" to show waveforms, because those drawn lines are
    % going to be the only repository of 'sniptime' and 'snipminmax')
    hMinMaxLines = nan(n_files,max_clusternum+1);
    for fileIndex = 1:n_files
      clabel = agglabel(clusternum{fileIndex}+1); %+1 for noise=0
      n_clusters = length(clabel);
      for clusterIndex = 1:n_clusters
        if ~isempty(clabel{clusterIndex})
          col = unique_color(clusterIndex,max_clusternum+1);
          tsnip_clust = t{fileIndex}(clabel{clusterIndex});
          % Instead of creating a separate line object for each min/max
          % pair, create a single line object and use NaNs to create
          % breaks in the line. This is much faster, since one only
          % needs to create a single handle
          x = tsnip_clust([1 1 1],:);
          y = [snipminmax{fileIndex}(:,clabel{clusterIndex}); ...
            nan(1,length(tsnip_clust))];
          hMinMaxLines(fileIndex,clusterIndex) = ...
            line(x(:),y(:),'Color',col,'Parent',hax);
        end
      end
    end
    nanFlag = isnan(hMinMaxLines);
    hMinMaxLines = hMinMaxLines(~nanFlag(:));
    setappdata(hax,'minmaxlines',hMinMaxLines);
    % Now call to set limits, visibility properties, etc.
    display_reconstruction(hax,'XLim',trange,'YLim',yrange)
  else
    %
    %% We're doing drawing and/or resetting axis limits
    %
    shc = getappdata(hax,'shc');
    n_files = length(shc);
    xlim0 = getappdata(hax,'XLim0');
    ylim0 = getappdata(hax,'YLim0');
    new_xlim = xlim0;
    new_ylim = ylim0;
    %% More argument parsing
    for argIndex = 2:2:nargin
      switch lower(varargin{argIndex})
        case 'xlim'
          new_xlim = varargin{argIndex+1};
        case 'ylim'
          new_ylim = varargin{argIndex+1};
        otherwise
          error(['Property name ' varargin{argIndex} ' not recognized']);
      end
    end
    %% Determine whether we're in min/max mode or in waveform mode
    spiketiming = (shc(1).sniprange(1):shc(1).sniprange(2))/shc(1).scanrate;
    minmaxmode = true;
    if (diff(new_xlim) < 1000*diff(spiketiming([1 end])))
      minmaxmode = false;
    end
    % Clear any old snippet lines
    hSnipLines = getappdata(hax,'sniplines');
    if ~isempty(hSnipLines)
      delete(hSnipLines);
      setappdata(hax,'sniplines',[]);
    end
    if minmaxmode
      % Make min/max lines visible
      hMinMaxLines = getappdata(hax,'minmaxlines');
      set(hMinMaxLines,'Visible','on');
    else
      % Make min/max lines invisible
      hMinMaxLines = getappdata(hax,'minmaxlines');
      set(hMinMaxLines,'Visible','off');
      % Draw snippet lines
      % Fetch waveforms in the visible region
      % Note we fetch all waveforms, not just the ones that are indexed
      file_starttime = getappdata(hax,'file_starttime');
      file_endtime = getappdata(hax,'file_endtime');
      snipIndex = getappdata(hax,'snipIndex');
      clusternum = getappdata(hax,'clusternum');
      max_clusternum = getappdata(hax,'max_clusternum');
      snipIndexRegion = cell(1,n_files);
      sniptime = cell(1,n_files);
      for fileIndex = 1:n_files
        % Find the appropriate part of the file
        if ~isempty(IntersectIntervals(...
            [file_starttime(fileIndex) file_endtime(fileIndex)],...
            new_xlim))
          trange_relative_in_scans = (new_xlim - file_starttime(fileIndex)) * ...
            shc(fileIndex).scanrate;
          snipIndexRegion{fileIndex} = find(...
            shc(fileIndex).sniptimes >= trange_relative_in_scans(1) & ...
            shc(fileIndex).sniptimes < trange_relative_in_scans(2));
          % Translate cluster numbers to new region
          clusternumtmp = nan(1,length(snipIndexRegion{fileIndex}));
          [commontimes,timeIndex,timeRegionIndex] = ...
            intersect(snipIndex{fileIndex},snipIndexRegion{fileIndex});
          clusternumtmp(timeRegionIndex) = clusternum{fileIndex}(timeIndex);
          clusternum{fileIndex} = clusternumtmp;
          % Prepare the snip times
          sniptimetmp = ...
            double(shc(fileIndex).sniptimes(snipIndexRegion{fileIndex}));
          sniptime{fileIndex} = sniptimetmp'/shc(fileIndex).scanrate + ...
            file_starttime(fileIndex) - min(file_starttime);
        end
      end
      sniptmp = sortheader_readsnips(shc,snipIndexRegion);
      % Loop over files and draw colorized spikes
      hSnipLines = [];
      for fileIndex = 1:n_files
        % First draw the "extra" (unclassified) spikes in a dashed gray
        % line
        nanFlag = isnan(clusternum{fileIndex});
        if any(nanFlag)
          n_to_draw = sum(nanFlag);
          twave = repmat(spiketiming',1,n_to_draw) + ...
            repmat(sniptime{fileIndex}(nanFlag),length(spiketiming),1);
          hSnipLines(end+1) = line(twave,sniptmp{fileIndex}(:,nanFlag),...
            'Color',[0.5 0.5 0.5],'LineStyle','--','Parent',hax);
        end
        % Now draw the other clusters. Offset by 1 so clusternum = 0
        % gets drawn as the first color
        tmpclusternum = clusternum{fileIndex};
        tmpclusternum(nanFlag) = -1;
        clabel = agglabel(tmpclusternum+1);
        for clusterIndex = 1:length(clabel)
          if (~isempty(clabel{clusterIndex}) && ~isempty(sniptime{fileIndex}))
            twave = repmat(spiketiming',1,length(clabel{clusterIndex})) + ...
              repmat(sniptime{fileIndex}(clabel{clusterIndex}),length(spiketiming),1);
            hSnipLines(end+1:end+length(clabel{clusterIndex})) = ...
              line(twave,sniptmp{fileIndex}(:,clabel{clusterIndex}),...
              'Color',unique_color(clusterIndex,max_clusternum+1),...
              'Parent',hax);
          end
        end % loop over clusters
      end % loop over files
      setappdata(hax,'sniplines',hSnipLines);
    end % if minmaxmode
    %% Set axis limits
    set(hax,'XLim',new_xlim,'YLim',new_ylim);
    % Might want to call "resetplotview(hax,'SaveCurrentView')"
    % here if we want this to interact nicely with the zoom functionality
    % of matlab.
  end % if new_window
