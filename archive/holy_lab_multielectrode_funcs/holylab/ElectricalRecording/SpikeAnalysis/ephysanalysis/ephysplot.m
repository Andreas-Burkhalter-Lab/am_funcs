function [out] = ephysplot(ephys,plotparams,hax)

% EPHYSPLOT: plot ephys data in a variety of ways
% Syntax:
%   out = ephysplot(ephys,plotparams,hax)
% where
%   ephys contains the input ephys data;
%   plotparams is a structure which specifies the plot type (see
%     EPHYSPLOTPARAMS); 
%   hax (optional) is the axis handle in which to create the plot
%     (defaults to current axis).
%   out is a structure containing random outputs necessary for eg
%     the annotation of publication plots
%
% See also: EPHYS, EPHYSPLOTPARAMS, EPHYSGUI.

% Copyright 2001 Timothy E. Holy
%
% Changelog
%  2004-05-03: Convert from using plotparams.number of channelnumber or
%               cellnumber
%  2006-03-26: Added in allowances for single celltimes vectors that aren't
%               within redundant cell braces
%  2006-12-04: (RCH) Adding in plotparams that allow more control of axis
%               labels (primarily useful for eg small bar plots of stimulus 
%               times along the bottom of a plot).  
%  2007-02-11: (RCH) overlay legend put under the control of plotparam
%               showtags
%  2007-03-05: (RCH) added option of outputting a structure; not elegent
%               useful in a way that decreases chances of publication
%               plotting errors - and shouldn't be destructive of anything
%  2007-06-08: (RCH) added option of a stimulus time offest.  Somewhat
%               dangerous (!!), but sometimes needed...

% Defaults
plotparams = default(plotparams,'tagFontSize',NaN);
plotparams = default(plotparams,'YTickOff',0);
plotparams = default(plotparams,'stimOnAlwaysSecond',0);
plotparams = default(plotparams,'autoYaxis',0);
plotparams = default(plotparams,'offsettime',0);
% currently, only used if plotting stimulus; could be added to other
% fields if so desired
% If the user didn't set a minimum time for computing rmax, then set that
% time to 0.
if isfield(plotparams,'raster_rmax_times')
  if length(plotparams.raster_rmax_times) < 2
    plotparams.raster_rmax_times(2) = 0;
  end
end

% Initialize output structure
out = struct;

% Axis preliminaries
if (nargin < 3)
  hax = gca;
end
if (~ishandle(hax) | ~strcmp(get(hax,'Type'),'axes'))
  error('Input axis is not valid');
end
% Delete any legends that are hanging around in this axis
hleg = findobj(0,'Tag','legend');
for i = 1:length(hleg)
  ud = get(hleg(i),'UserData');
  if (ud.PlotHandle == hax)
    delete(hleg(i));
  end
end
% This is a cla for axis hax
delete(findobj(get(hax,'Children'),'flat','HandleVisibility','on'))

% Tags preliminaries. Must junk those tags which
% do not pull up any data, and warn the user.
if ~isfield(plotparams,'tags')
  return;
end
if ischar(plotparams.tags)
  tags = {plotparams.tags};
else
  tags = plotparams.tags;
end
ntags = length(tags);
tagindex = cell(1,ntags);
emptytags = zeros(1,ntags);
for i = 1:ntags
  tagindex{i} = strmatch(tags{i},{ephys.tag},'exact')';
  if isfield(plotparams,'alltags')
    tagnum(i) = strmatch(tags{i},plotparams.alltags,'exact');
  end
  emptytags(i) = isempty(tagindex{i});
end
emptytags = find(emptytags);
if ~isempty(emptytags)
  %uiwait(msgbox(['These tags are empty: ',tags{emptytags}],'','warn'));
  msgbox(['These tags are empty: ',tags{emptytags}],'','warn');
end
tagindex(emptytags) = [];
tagnum(emptytags) = [];
tags(emptytags) = [];
ntags = length(tagindex);
if (ntags == 0)
  return;
end

% Figure out whether we're doing channels or cells, and compute the
% channel or cell index for the appropriate data type
indexchancell = [];
if (supportlegacyppnumber & isfield(plotparams,'number'))
  indexchancell = plotparams.number;
else
  if strmatch(plotparams.fieldtoplot,{'sniptimes','snippets','wave','envelope'},'exact')
    if isfield(ephys,'channels')
      indexchancell = findainb(plotparams.channelnumber,ephys(1).channels);
    else % assume first - often has been recreated anyways
      indexchancell = 1;
    end
  elseif strmatch(plotparams.fieldtoplot,{'celltimes'},'exact')
    indexchancell = findainb(plotparams.cellnumber,ephys(1).cellnums);
  end
end

% Get colors: either a consistent pattern for each tag,
% or the default color order
if (isfield(plotparams,'alltags') & isfield(plotparams,'tagcolors'))
  co = plotparams.tagcolors(tagnum,:);
else
  co = get(hax,'ColorOrder');
end
ncol = size(co,1);
% Clean up from previous axis setting
set(hax,'YLimMode','auto','YTickMode','auto','YTickLabelMode','auto', ...
  'TickLength',[0.02 0.03],'TickDir','out','Tag','ephysplotax');
% Add any user defined axis requests
if plotparams.YTickOff
  set(hax,'YTick',[]);
end

% Check that intervals are consistent
[duration,toff] = ephysvalidatetimerange(ephys([tagindex{:}]));

% Labels
labels = tags;
if isfield(plotparams,'tagstext')
  % Use user-specified labels instead
  % Have to make the match
  for i = 1:length(tags)
    labelindex = strmatch(tags{i},plotparams.alltags,'exact');
    if ~(length(labelindex) == 1)
      error('Mismatch between plotparams.alltags and plotparams.tags');
    end
    labels{i} = plotparams.tagstext{labelindex};
  end
  if ischar(labels)
    labels = {labels};
  end
end

% Vertical spacing preprocessing and boundary-line plotting:
% useful for raster, spikepeaks, envelope, snipreconstruct, wave plots
% Leaves a gap for the different groups
% This does nothing for other kinds of plots
vertp = ephysplotvertspacing(ephys,plotparams,labels,tagindex,hax,[0 duration]+toff);
if (length(vertp) > 0)
  yc = vertp{1};
  ylimset = vertp{2};
end

% Do the plot
% Note: all plotting commands will be parented to given hax,
% in case it's not the current axis.
% The only bad thing here is that some of the plotting commands
% (stairs and errobar) do not allow one to set the parent at the time,
% but only after the object has been made. If the target axis'
% NextPlot property is replace (or replacechildren), then this
% will tank their contents. So, have to protect the current axis
% by setting it's NextPlot property to add.
if (gca ~= hax)
  curraxnp = get(gca,'NextPlot');
  set(gca,'NextPlot','add');
end
switch plotparams.fieldtoplot
  case 'stimulus',
    if (~isfield(plotparams,'type') | ~strcmp(plotparams.type,'bar'))
      hline = zeros(1,ntags);
      cnp = get(hax,'NextPlot');
      set(hax,'NextPlot','add');
      for i = 1:ntags
        colindx = mod(i-1,ncol)+1;
        if (isfield(plotparams,'type') & strcmp(plotparams.type,'first'))
          ntoplot = 1;
        else
          ntoplot = length(tagindex{i});
        end
        repeats = 0;
        if (isfield(plotparams,'type') & strcmp(plotparams.type,'repeats'))
          repeats = 1;
        end
        for j = 1:ntoplot
          cephys = ephys(tagindex{i}(j));
          t = (cephys.stimulus(2,:) - cephys.scanrange(1))/cephys.scanrate + toff;
          if repeats
            maxs = max(cephys.stimulus(1,:));
            hline(i) = stairs(t,0.8*cephys.stimulus(1,:)/maxs + yc{i}(j));
          else
            hline(i) = stairs(t,cephys.stimulus(1,:));
          end
          % (I did hline(i) to save only one line/color for legend command below)
          set(hline(i),'Color',co(colindx,:),'Parent',hax);
          if isfield(plotparams,'objectproperties')
            set(hline(i),plotparams.objectproperties{:});
          end
        end
      end
      if repeats
        set(hax,'YLim',ylimset);
      else
        ylim = get(hax,'YLim');
        set(hax,'YLim',ylim+[-.1 .1]);
        if (isfield(plotparams,'showtags') & plotparams.showtags)
          legend(hax,hline,labels{:});
        end
      end
      set(hax,'NextPlot',cnp);
    else
      % Do a bars plot of the stimulus instead
      hpatch = zeros(1,ntags);
      spacing = 1;
      if isfield(plotparams,'spacing')
        spacing = plotparams.spacing;
      end
      for i = 1:ntags
        colindx = mod(i-1,ncol)+1;
        cephys = ephys(tagindex{i}(1));
        t = (cephys.stimulus(2,:) - cephys.scanrange(1))/cephys.scanrate + ...
          toff + plotparams.offsettime;
        % Find the index where valve turned on
        if plotparams.stimOnAlwaysSecond % means data's been transformed past where pre-stim times are negative, and should just use second index as "on"...
          dtmin = t(2);
          ion = 2;
        else % means normal method should work
          [dtmin,ion] = min(abs(t));
        end
        hpatch(i) = patch(t([ion ion+1 ion+1 ion ion]), [0 0 1 1 0]*.7 / ...
          spacing + ntags-i, co(colindx,:), 'Parent', hax,'EdgeColor','none',...
          'Tag','bar');
      end
      if isfield(plotparams,'objectproperties')
        set(hpatch(i),plotparams.objectproperties{:});
      end
      set(hax,'YLim',[-.1 ntags+.1]);
      if (isfield(plotparams,'showtags') & plotparams.showtags)
        allx = get(hpatch,'XData');
        if iscell(allx)
          allx = [allx{:}];
        end
        maxx = max(max(allx));
        for i = 1:ntags
          ypos = unique(get(hpatch(i),'YData'));
          axes(hax)
          ht = text(1.1*maxx,mean(ypos),labels{i},'Tag','bartext');
          if ~isnan(plotparams.tagFontSize)
            set(ht,'FontSize',plotparams.tagFontSize);
          end
        end
      end
    end

  case {'sniptimes','celltimes'}
    type = 'raster';
    if isfield(plotparams,'type')
      type = plotparams.type;
    end
    check_celltimes = strcmp(plotparams.fieldtoplot,'celltimes') && ...
      isfield(ephys,'cellscantimes');
    switch type
      case 'raster',
        showmax = false;
        if isfield(plotparams,'raster_rmax_times')
          showmax = true;
        end
        set(hax,'YLim',ylimset);
        % Plot the rasters
        for i = 1:ntags
          for j = 1:length(tagindex{i})
            cephys = ephys(tagindex{i}(j));
            %spikes = getfield(cephys,plotparams.fieldtoplot,{plotparams.number});
            if indexchancell==1 & ~iscell(cephys.(plotparams.fieldtoplot))
              spikes = cephys.(plotparams.fieldtoplot)';
            else
              spikes = cephys.(plotparams.fieldtoplot){indexchancell}';
            end
            if size(spikes,1)>size(spikes,2)
              spikes = spikes';
            end
            %spikes = spikes{1}(:)';
            nspikes = length(spikes);
            thisy = [ones(1,nspikes)+0.2;ones(1,nspikes)-0.2] + yc{i}(j)-1;
            x = repmat((spikes-cephys.scanrange(1))/cephys.scanrate + toff,2,1);
            colindx = mod(i-1,ncol)+1;
            if check_celltimes && ...
                ~isequal(cephys.scanrange,cephys.cellscanrange{indexchancell})
              % Draw a black line indicating that this repeat is not valid
              x = cephys.scanrange;
              x = (x-x(1))/cephys.scanrate + toff;
              thisy = [0 0] + yc{i}(j);
              line(x,thisy,'Color','k');
            else
              hline = line(x,thisy,'Color',co(colindx,:),...
                'LineWidth',0.25,'Tag','Rast',...
                'Parent',hax);
              % objectproperties? But what if linewidth is a parameter?
              if isfield(plotparams,'objectproperties')
                set(hline,plotparams.objectproperties{:});
              end
            end
            if showmax && ~isempty(x)
%               [rmax,tmax] = deltarmax_calculation(x(1,:) - plotparams.raster_rmax_times(1),plotparams.raster_rmax_times(2));
%               line([1 1]*tmax,thisy(:,1)+[-1;1]*0.1*diff(thisy(:,1)),'Color','r','LineWidth',2);
            end
          end
        end
      case {'PSTH', 'PSTH w/ sem'},
        ebars = 0;
        if strcmp('PSTH w/ sem',plotparams.type)
          ebars = 1;
        end
        binwidth = 1;
        if isfield(plotparams,'binwidth')
          binwidth = plotparams.binwidth;
        end
        cnp = get(hax,'NextPlot');  % errorbar will execute newplot
        set(hax,'NextPlot','add');
        hline = zeros(1,ntags);
        herr = zeros(1,ntags);
        maxRplusErr = NaN;
        for i = 1:ntags
          colindx = mod(i-1,ncol)+1;        % Wrap-around color indices
          [r,t,err] = ephyspsth(ephys(tagindex{i}),binwidth,indexchancell, ...
            plotparams.fieldtoplot);
          maxRplusErr = max([maxRplusErr max(r+err)]);
          if ~ebars
            hline(i) = plot(t,r,'Color',co(colindx,:),'Parent',hax);
          else
            if(str2num(version('-release')) <=13 )
              htemp = errorbar(t,r,err);
            else
              htemp = errorbar('v6',t,r,err);
            end
            set(htemp,'Color',co(colindx,:),'Parent',hax);
            hline(i) = htemp(2);
            herr(i) = htemp(1);
          end
          out.psth_line_by_itag{i} = r;
          out.psth_errors_by_itag{i} = err;
          out.psth_timebins_by_itag{i} = t;
        end
        set(hax,'NextPlot',cnp);
        if isfield(plotparams,'objectproperties')
          set(hline,plotparams.objectproperties{:});
        end
        if (isfield(plotparams,'showtags') & plotparams.showtags)
          hleg = legend(hax,hline,labels{:});
          if isfield(plotparams,'legtextfs')
            set(hleg,'FontSize',plotparams.legtextfs)
          end
        end
        % Set tags on object for future manipulation
        set(hline,'Tag','PSTH');
        if ebars
          set(herr,'Tag','PSTHebars');
        end
        % Make the lower y limit = 0
        ylim = get(hax,'YLim');
        set(hax,'YLim',[0 ylim(2)]);
        % If plot params "autoYaxis" is on, make the Yaxis some number times the
        % max r+err value
        if plotparams.autoYaxis
          set(hax,'YLim',[0 plotparams.autoYaxis*maxRplusErr]);
        end
      otherwise,
        error('Unrecognized type')
    end % switch type
  case 'wave'
    % Axis preliminaries
    set(hax,'YLim',ylimset);
    set(hax,'NextPlot','add');
    scalefac = 1;
    if isfield(plotparams,'tomicrovolts')
      scalefac = plotparams.tomicrovolts;
    end
    % Plot the waveforms
    for i = 1:ntags
      for j = 1:length(tagindex{i})
        cephys = ephys(tagindex{i}(j));
        indx = indexchancell;
        thisy = cephys.wave(indx,:)*scalefac + yc{i}(j);
        x = linspace(0,duration,size(thisy,2)) + toff;
        colindx = mod(i-1,ncol)+1;
        h = plot(x,thisy,'Parent',hax);
        set(h,'Color',co(colindx,:),'Tag','wave');
        if isfield(plotparams,'objectproperties')
          set(h,plotparams.objectproperties{:});
        end
      end
    end
    set(hax,'NextPlot','replacechildren');
  case 'envelope'
    % Axis preliminaries
    set(hax,'YLim',ylimset);
    set(hax,'NextPlot','add');
    scalefac = 1;
    if isfield(plotparams,'tomicrovolts')
      scalefac = plotparams.tomicrovolts;
    end
    % Plot the envelopes
    for i = 1:ntags
      for j = 1:length(tagindex{i})
        cephys = ephys(tagindex{i}(j));
        indx = 2*indexchancell-1;
        thisy = cephys.envelope([indx indx+1],:)*scalefac + yc{i}(j);
        x = linspace(0,duration,size(thisy,2)) + toff;
        colindx = mod(i-1,ncol)+1;
        h = fillmm2(thisy(1,:),thisy(2,:),x,'Parent',hax);
        set(h,'FaceColor',co(colindx,:),'EdgeColor',co(colindx,:),'Tag','envelope');
        if isfield(plotparams,'objectproperties')
          set(h,plotparams.objectproperties{:});
        end
      end
    end
    set(hax,'NextPlot','replacechildren');
  case 'snippets'
    % Note this will need to be generalized if we're plotting
    % for cells rather than for channels
    if (isfield(plotparams,'type') & strcmp(plotparams.type,'cell'))
      error('Not implemented yet');
    end
    firstindx = tagindex{1}(1);
    sniprange = ephys(firstindx).sniprange;
    scanrate = ephys(firstindx).scanrate;
    scalefac = 1;
    if isfield(plotparams,'tomicrovolts')
      scalefac = plotparams.tomicrovolts;
    end
    % Time in milliseconds
    t = (sniprange(1):sniprange(2))/scanrate * 1000;
    % Use tag colors if there's more than one tag
    colorbytag = (ntags > 1);
    % If reconstructing, adjust
    if (isfield(plotparams,'type') & strcmp(plotparams.type,'reconstruct'))
      t0 = t/1000;   % save for offsets
      colorbytag = 1;
    end
    hlineleg = zeros(1,ntags);
    % (t is in milliseconds)
    for i = 1:ntags
      for j = 1:length(tagindex{i})
        cephys = ephys(tagindex{i}(j));
        if ~isempty(cephys.snippets{indexchancell})
          if (isfield(plotparams,'type') & strcmp(plotparams.type,'reconstruct'))
            [sizesnip,nsnips] = size(cephys.snippets{indexchancell});
            ttemp = cephys.sniptimes{indexchancell}(:)';
            ttemp = (ttemp - cephys.scanrange(1))/scanrate + toff;
            t = repmat(ttemp,sizesnip,1) + repmat(t0',1,nsnips);
            thisy = cephys.snippets{indexchancell}*scalefac + yc{i}(j);
            if isfield(cephys,'snipthresh')
              thisthresh = cephys.snipthresh(:,indexchancell);
              patch([0 duration duration 0 0]+toff,...
                thisthresh([1 1 2 2 1])*scalefac + yc{i}(j),...
                [0.9 0.9 0.9],'Parent',hax,'EdgeColor','none','tag','threshbar');
            end
          else
            thisy = scalefac*cephys.snippets{indexchancell};
          end
          hline = line(t,thisy,'Parent',hax,'Tag','snippet');
          hlineleg(i) = hline(1);
          if (colorbytag)
            colindx = mod(i-1,ncol)+1;
            set(hline,'Color',co(colindx,:));
          end
          % ? objectproperties?
        end
      end
    end
    if (~isfield(plotparams,'type') | ~strcmp(plotparams.type,'reconstruct'))
      if (isfield(plotparams,'showtags') & plotparams.showtags)
        legend(hax,hlineleg,labels{:});
      end
      % A hack so the xlim gets set properly later
      toff = t(1);
      duration = t(end)-toff;
    else
      set(hax,'YLim',ylimset);
      if (ntags == 1)
        set(hax,'YTickMode','auto','YTickLabelMode','auto');
      end
    end
  otherwise
    error('Unrecognized fieldtoplot')
end

% General axis properties
set(hax,'XLim',[0 duration] + toff);
if isfield(plotparams,'axisproperties')
  set(hax,plotparams.axisproperties{:});
end
if isfield(plotparams,'titlenumber')
  if plotparams.titlenumber
    if ischar(plotparams.titlenumber)
      str = plotparams.titlenumber;
    elseif (any(findstr('snip',plotparams.fieldtoplot)) | ~isfield(ephys,'celltimes'))
      % @note: here comes the limitation:
      %    we assume all files specified in ephys have same channels field.
      %    So ephys(1) is used to convert channel idx to channel.
      num = plotparams.channelnumber;
      str = 'Channel';
    else
      str = 'Cell';
      num = plotparams.cellnumber;
    end
    htitle = title([str,' ',num2str(num)]);
    set(htitle,'Parent',hax)
  end
end
% Now undo the protection of the current axis
if (gca ~= hax)
  set(gca,'NextPlot',curraxnp);
end

if exist('hline','var')
  out.hline = hline; % useful for doing a legend with more control
  out.labels = labels;
end


  function outdata = ephysplotvertspacing(ephys,plotparams,labels,tagindex,hax,xlim)
    % Do the grungy work of deciding how to space repeats vertically
    % outputs: {yc,ylim}
    outdata = {};
    % Define the "fieldtoplot" which may need this type of analysis
    vertppfield = {'stimulus','sniptimes','celltimes','envelope','snippets','wave'};
    % For each field, define the plot types which require this analysis
    vertpptype = {{'repeats'},{'raster'},{'raster'},{},{'reconstruct'},{}};
    % If the plot type is not defined, should we do this analysis? Specify
    % 'yes' by setting the corresponding entry in vertppbydefault to 1.
    vertppbydefault = [0 1 1 1 0 1];
    % Define the plot types which need spacing corresponding to physical units
    vertpprescale = [0 0 0 1 1 1];
    ppindx = strmatch(plotparams.fieldtoplot,vertppfield,'exact');
    % Get the channel/cell index
    indexchancell = [];
    if (supportlegacyppnumber & isfield(plotparams,'number'))
      indexchancell = plotparams.number;
    else
      if strmatch(plotparams.fieldtoplot,{'sniptimes','snippets','wave','envelope'},'exact')
        if isfield(ephys,'channels')
          indexchancell = findainb(plotparams.channelnumber,ephys(1).channels);
        else % assume 1st, probably recreated...
          indexchancell = 1;
        end
      elseif strmatch(plotparams.fieldtoplot,{'celltimes'},'exact')
        indexchancell = findainb(plotparams.cellnumber,ephys(1).cellnums);
      end
    end
    if (ppindx)
      if ( (vertppbydefault(ppindx) & ~isfield(plotparams,'type')) | ...
          ( isfield(plotparams,'type') & (isempty(vertpptype{ppindx}) | ...
          strmatch(plotparams.type,vertpptype{ppindx},'exact')) ) )
        % Figure out how many repeats in each group
        ntags = length(tagindex);
        nrpts = zeros(1,ntags);
        for j = 1:ntags
          nrpts(j) = length(tagindex{j});
        end
        cnrpts = cumsum(nrpts);
        nrptstot = cnrpts(end);
        bindx = [0,cnrpts];
        % Leave a gap of one unit between groups
        bshift = [0,cnrpts + (1:ntags)];
        yc = cell(1,ntags);
        % Adjust the spacing
        spacing = 1;
        if isfield(plotparams,'spacing')
          spacing = plotparams.spacing;
        end
        %spacing = spacing*1.2;
        % Unit scale is OK for several plot types, but
        % not for envelopes or snipreconstruct
        if vertpprescale(ppindx)
          scalefac = 1;
          if isfield(plotparams,'tomicrovolts')
            scalefac = plotparams.tomicrovolts;
          end
          if ~isfield(plotparams,'fixedspacing')
            alltagindx = [tagindex{:}];
            fmin = zeros(1,length(alltagindx));
            fmax = zeros(1,length(alltagindx));
            for i = 1:length(alltagindx)
              if strcmp(plotparams.fieldtoplot,'envelope')
                indx = 2*indexchancell-1;
                fmin(i) = min(ephys(alltagindx(i)).envelope(indx,:));
                fmax(i) = max(ephys(alltagindx(i)).envelope(indx+1,:));
              elseif strcmp(plotparams.fieldtoplot,'wave')
                indx = indexchancell;
                fmin(i) = min(ephys(alltagindx(i)).wave(indx,:));
                fmax(i) = max(ephys(alltagindx(i)).wave(indx,:));
              elseif strcmp(plotparams.fieldtoplot,'snippets')
                indx = indexchancell;
                tmp = ephys(alltagindx(i)).snippets{indx};
                if ~isempty(tmp)
                  fmin(i) = min(min(tmp));
                  fmax(i) = max(max(tmp));
                else
                  fmin(i) = nan;
                  fmax(i) = nan;
                end
              end
            end
            fmin = scalefac*fmin;
            fmax = scalefac*fmax;
            if (scalefac < 0)
              tmp = fmin;
              fmin = fmax;
              fmax = tmp;
            end
          else
            fmin = plotparams.fixedspacing(1);
            fmax = plotparams.fixedspacing(2);
          end
          range = max(fmax-fmin);
          spacing = spacing*range;
          bottomborder = fmin(end);%fmin(1);  % The current version works for scalefac < 0
          topborder = fmax(1);%fmax(end);
        else
          bottomborder = -spacing/2;
          topborder = spacing/2;
        end
        % Set the actual positions
        for i = 1:ntags
          yc{i} = (bshift(i) + (1:nrpts(i)))*spacing;
        end
        % Now reverse direction, so the first one is plotted at
        % the top of the axis
        % (better than setting the axis direction to reverse, because
        %  some will have vertical scale ticks)
        ymax = max([yc{:}]);
        for i = 1:ntags
          yc{i} = ymax - yc{i};
        end
        y = [yc{:}];
        % Prepare for labelling
        ym = zeros(1,ntags);
        for i = 1:ntags
          ym(i) = mean(yc{i});
        end
        labelord = length(ym):-1:1;
        % Plot the boundary lines
        for i = 1:ntags-1
          line(xlim,[0 0] + (yc{i}(end)+yc{i+1}(1))/2,'Color','k','LineStyle',':',...
            'Tag','Boundary','Parent',hax);
        end
        % Label the groups, if desired
        if (isfield(plotparams,'showtags') & plotparams.showtags)
          set(hax,'YTick',ym(labelord),'YTickLabel',labels(labelord),'TickLength',[0 0]);
        end
        % Y limits
        ylim = [min(y) + bottomborder, max(y) + topborder];
        if strcmp(plotparams.fieldtoplot,'stimulus')
          ylim(2) = ylim(2)+spacing;
        end
        if any(isnan(ylim))
          ylim = [0 eps];
        end
        % Return remaining parameters
        outdata{1} = yc;
        outdata{2} = ylim;
      end
    end
