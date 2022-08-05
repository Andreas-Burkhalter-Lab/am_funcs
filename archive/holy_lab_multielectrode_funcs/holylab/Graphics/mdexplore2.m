function [indx, aborted] = mdexplore2(xproj,xorig,options)
% MDEXPLORE: explore multidimensional data interactively
% This has two modes of operation:
%   explore: plots points in 2d, and clicking on a point displays the
%     corresponding multidimensional data
%   cluster: plot points in 2d, and user can draw a polygon around the
%     points; the function returns the index of the encircled points.
%     (left-drag to start polygon drawing, left click to continue, middle
%     click to finish.)
%
% Syntax for 'explore' mode (does not block command line):
%   mdexplore2(xproj,xorig)
%   mdexplore2(xproj,xorig,options)
%
% Syntax for 'cluster' mode (blocks command line until ctrl+q or ctrl+a is pressed):
%   [indx, aborted] = mdexplore2(xproj,xorig)
%   [indx, aborted] = mdexplore2(xproj,xorig,options)
%
% PRE:
%   xproj is a 2-by-N matrix containing projections into two dimensions;
%   xorig is a cell array containing the original multidimensional data;
%   options is a structure with the following fields:
%     plotfunc: a function handle which knows how to plot the
%       multidimensional data (default: @plot).
%     markersize: size of points to plot. Default depends on the mode,
%       larger points (size 9) for explore mode and smaller (size 4) for
%       cluster mode.
%     fignum: which figure to use (default: create a new figure)
%     subplotnum: if should be put in a subplot, what numbers to use
%     axisinfo: allows you to set axis limits to match another figure
%     separationPositions: the x positions to plot separation lines. 
% POST:
%   aborted: 1 when user closed the figure;
% 
% See also: MDEXPLORE_EPHYS.
  
  if (nargin < 3)
    options = struct;
  end
  mode = 'cluster';
  if (nargout == 0)
    mode = 'explore';
  end
  if ~isfield(options,'plotfunc')
    options.plotfunc = @plot;
  end
  if ~isfield(options,'markersize')
    switch mode
      case 'cluster'
        options.markersize = 4;
      case 'explore'
        options.markersize = 4;
    end
  end
  options = default(options,'win_locations',nan);
  
  default_options('separationPositions', []);
  separationPositions=options.separationPositions;

  [d,N] = size(xproj);
  if (d ~= 2)
    error('Input projection data must be a 2-by-N array');
  end

  if isfield(options,'fignum')
    figure(options.fignum)
  else
    figure
  end
  fig=gcf;
  if isfield(options,'subplotnum')
    t = options.subplotnum;
    subplot(t(1),t(2),t(3))
  end
  hax = gca;
  hpoints = line(xproj([1 1],:),xproj([2 2],:),...
    'MarkerFaceColor','k',...
    'MarkerEdgeColor','k',...
    'Marker','o',...
    'MarkerSize',options.markersize);
  if isfield(options,'axisinfo')
    axis(options.axisinfo);
  else
     xlim=get(hax, 'xlim');
     if(xlim(2)<0) xlim(2)=0; end
     if(xlim(1)>0) xlim(1)=0; end
     margin=diff(xlim)*0.05;
     set(hax, 'xlim', xlim+[-1 1]*margin);
     ylim=get(hax, 'ylim');
     if(ylim(2)<0) ylim(2)=0; end
     if(ylim(1)>0) ylim(1)=0; end
     margin=diff(ylim)*0.05;
     set(hax, 'ylim', ylim+[-1 1]*margin);
     
     % draw a light + at (0,0)
     hCrossLines=line([xlim', [0;0]], [[0;0], ylim']);
     set(hCrossLines, 'color', 0.8*[1 1 1]);
     set(hCrossLines, 'hitTest', 'off'); % to avoid these lines to cover points
  end
  
  % If wanting a return index, set up polygon drawing and use uiwait to
  % block until the user draws a polygon; then return the indices of the
  % points within the polygon

  % set up a callback so that when user clicks on a point, it
  % draws the corresponding md data
  setappdata(hax,'mdplotfcn',options.plotfunc);
  % Set up the callback cell array
  cbarray = cell(N,1);
  for i = 1:N
     cbarray{i} = {@mdexplore_plotcb,xorig{i}};
  end
  set(hpoints,{'ButtonDownFcn'},cbarray);
  %  set(hax,'HitTest','off');
  
  if (nargout > 0)
    % set(hpoints,'HitTest','off');
    bind_shortcut(fig, 'ctrl+q', @onDone);
    bind_shortcut(fig, 'ctrl+a', @onDone);
    helpMsg=sprintf('Press ctrl+q or ctrl+a to go to next raw cluster; \nClick "X" to exit');
    bind_shortcut(fig, 'f1', {@show_help, 'clustering', helpMsg});
    
    % add context menu
    cmenu = uicontextmenu;
    set(hax,'UIContextMenu',cmenu);
    uimenu(cmenu,'Label','delete clusters', ...
       'Callback',{@fDelClusters, struct('axes', hax)});
    uimenu(cmenu,'Label','show snippets', ...
       'Callback',{@fShowSnippets, struct('axes', hax)});
    
    % show tip about how to exit
    addTip(hax);
    
    set(fig, 'CloseRequestFcn', @cleanup);
    set(hax,'ButtonDownFcn',@mdexplore_polygon);
    setappdata(fig, 'xproj', xproj);
    setappdata(fig, 'xorig', xorig);
    setappdata(fig, 'hpoints', hpoints);
    setappdata(fig, 'clusters', {});
    setappdata(fig, 'separationPositions', separationPositions);
    setappdata(fig, 'win_locations',options.win_locations);
    set(fig, 'UserData', 0);
    waitfor(fig,'UserData');
    if(ishandle(hax))
       clusters=getappdata(fig, 'clusters');
       indx=clusters;
       aborted=0;
    else
       % user closed the figure, which means "abort"
       indx={};
       aborted=1;
    end
  end

function addTip(hax)
   hText=text(mean(get(gca, 'xlim')), mean(get(gca, 'ylim')), 'press F1 for help');
   set(hText, 'color', 'red');
   set(hText, 'fontSize', 20);
   pos=get(hText, 'extent');
   set(hText, 'position', [pos(1)-pos(3)/2 pos(2)]);
   t = timer('TimerFcn', {@tipTimerCallback, struct('hText', hText)}, ...
      'ExecutionMode', 'singleshot', 'startdelay', 3);
   start(t);

function tipTimerCallback(sender, event_data, args)
   hText=args.hText;
   if(ishandle(hText)) % make sure handle is still valid (figure isn't closed yet)
      delete(hText);
   end

function fDelClusters(sender, event_data, args)
   curAxes=args.axes;
   fig=get_parent_fig(curAxes);
   clusterColors=getappdata(fig, 'clusterColors');
   strings=cell(1, length(clusterColors));
   for idx=1:length(strings)
      strings{idx}=['cluster ' num2str(idx)];
   end
   selectedClusters=selection_dialog(strings, clusterColors);
   clusters=getappdata(fig, 'clusters');
   ptsIndices=[clusters{selectedClusters}];
   clusters(selectedClusters)=[];
   setappdata(fig, 'clusters', clusters);
   clusterColors(selectedClusters)=[];
   setappdata(fig, 'clusterColors', clusterColors);
   
   hpoints=getappdata(fig, 'hpoints');
   color='black';
   set(hpoints(ptsIndices), 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
   
function fShowSnippets(sender, event_data, args)
   curAxes=args.axes;
   fig=get_parent_fig(curAxes);
   clusterColors=getappdata(fig, 'clusterColors');
   strings=cell(1, length(clusterColors));
   for idx=1:length(strings)
      strings{idx}=['cluster ' num2str(idx)];
   end
   if length(clusterColors)>1
       selectedClusters=selection_dialog(strings, clusterColors); % pick clusters to view snippets
   elseif length(clusterColors)==1
       selectedClusters=1;
   else
       return
   end
   if(isempty(selectedClusters))
      return;
   end
   
   clusters=getappdata(fig, 'clusters');
   ptsIndices=[clusters{selectedClusters}];
   figSnip=figure;
   window_arrangement_profile = getappdata(fig,'window_arrangement_profile');
   if ~strcmp(window_arrangement_profile,'none')
      win_locations = getappdata(fig,'win_locations');
      if isstruct(win_locations)
        set(figSnip,'Position',win_locations.snip)
      end
   end
   xorig=getappdata(fig, 'xorig');
   snip=cat(2, xorig{:});
   hLines=plot(snip(:, ptsIndices));
   
   % if more than one clusters, use clusters' colors as lines' colors
   if(length(selectedClusters)>1)
      tNumOfLines=0;
      for tIdx=1:length(selectedClusters)
         pickedLines=hLines(tNumOfLines+[1:length(clusters{selectedClusters(tIdx)})]);
         set(pickedLines, 'color', clusterColors{selectedClusters(tIdx)});
         tNumOfLines=tNumOfLines+length(pickedLines);
      end % for,
   end % if, more than 1 cluster picked
   
   title(['Cluster(s): ' num2str(selectedClusters)]);
   
   set(hLines, 'ButtonDownFcn', @onClickOnSnip);
   axesSnip=gca;
   set(axesSnip, 'ButtonDownFcn', @onSelectMultipleSnip);
   
   % plot dot lines to separate channels
   separationPositions=getappdata(fig, 'separationPositions');
   plotSeparationLines(axesSnip, separationPositions);
   
%    nSamplesTotal=size(snip,1);
%    nSamplesPerChannel=64; % TODO: hardcoded number 64
%    ylim=get(axesSnip, 'ylim');
%    hold(axesSnip, 'on');
%    if(~isempty(nSamplesPerChannel:nSamplesPerChannel:nSamplesTotal-1))
%       hDotLines=plot(repmat(nSamplesPerChannel:nSamplesPerChannel:nSamplesTotal-1, 2,1), ylim');
%       set(hDotLines, 'lineStyle', ':', 'color', 'black', 'hitTest', 'off');
%    end
   
   % set appdata then waitfor ctrl+q
   setappdata(figSnip, 'ptsIndices', ptsIndices);
   setappdata(figSnip, 'hLines', hLines);
   setappdata(figSnip, 'parentFig', fig);
   bind_shortcut(figSnip, 'ctrl+q', @onDone);
   bind_shortcut(figSnip, 'ctrl+a', @onDone); 
    % since ctrl+q, if missassigned to a different matlab window will kill
    % the session
   bind_shortcut(figSnip, 'delete', @onRemoveSnipFromClusters);
   bind_shortcut(figSnip, 'ctrl+d', @onRemoveSnipFromClusters);
   helpMsg=sprintf('press ctrl+q or ctrl+a to close the figure; \nPress DEL or ctrl+d to remove snippets from the cluster(s)');
   bind_shortcut(figSnip, 'f1', {@show_help, 'clustering', helpMsg});
   set(figSnip, 'UserData', 0);
   waitfor(figSnip,'UserData');
   free(figSnip);

   

function onSelectMultipleSnip(sender, event)
   axesSnip=sender;
   figSnip=get_parent_fig(axesSnip);
   rect = GetSelRect;
   % loop over all lines and test if a line is in the rect
   hLines=getappdata(figSnip, 'hLines'); % OR: findobj() can be used.
   for idxLine=1:length(hLines)
      hLine=hLines(idxLine);
      if(is_line_in_rect(hLine, rect))
         toggleSelection(hLine);
      end
   end % for, each line

   
function onRemoveSnipFromClusters(sender, event)
   figSnip=sender;
   hLines=getappdata(figSnip, 'hLines'); % the line objs of selected clusters
   ptsIndices=getappdata(figSnip, 'ptsIndices'); % the indices of points that forms the selected clusters
   parentFig=getappdata(figSnip, 'parentFig');
   % NOTE: first you select some clusters, then you select a few snippets
   %       from these selected clusters 
   hSelectedLines=findobj(figSnip, 'type', 'line', 'Marker', '*'); % find selected lines in selected clusters
   if(isempty(hSelectedLines))
      return;
   end
   selectedPtsIndices=ptsIndices(indexainb(hLines, hSelectedLines));
   hpoints=getappdata(parentFig, 'hpoints'); % the handles of all points in pca axes
   hSelectedPoints=hpoints(selectedPtsIndices); % the handles of selected points in selected clusters
   
   % now got all info to delete the selected snippets from selected
   % clusters, so we can do deletion and updating
   
   % On snippet view axes, remove the line objs:
   delete(hSelectedLines);
   % update appdata of figSnip
   lineIdxToDel=indexainb(hLines, hSelectedLines);
   hLines(lineIdxToDel)=[];
   ptsIndices(lineIdxToDel)=[]; % NOTE: in this way, ptsIndices(i) is still aligned w/ hLines(i)
   setappdata(figSnip, 'hLines', hLines);
   setappdata(figSnip, 'ptsIndices', ptsIndices);
   
   % on pca axes, show selected points as non-clustered
   color='black';
   set(hSelectedPoints, 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'marker', 'o');
   % update appdata of parentFig (pca figure): clusters and clusterColors
   clusters=getappdata(parentFig, 'clusters');
   clusterColors=getappdata(parentFig, 'clusterColors');
   emptyClusterIndices=[];
   for idxCluster=1:length(clusters)
      ptsIndicesForCurCluster=clusters{idxCluster};
      ptsIndicesForCurCluster=setdiff(ptsIndicesForCurCluster, selectedPtsIndices);
      if(isempty(ptsIndicesForCurCluster))
         emptyClusterIndices(end+1)=idxCluster;
      else
         clusters{idxCluster}=ptsIndicesForCurCluster;
      end
   end
   clusters(emptyClusterIndices)=[];
   clusterColors(emptyClusterIndices)=[];
   setappdata(parentFig, 'clusters', clusters);
   setappdata(parentFig, 'clusterColors', clusterColors);

   
function toggleSelection(hLine)   
   figSnip=get_parent_fig(hLine);
   hLines=getappdata(figSnip, 'hLines');
   ptsIndices=getappdata(figSnip, 'ptsIndices');
   parentFig=getappdata(figSnip, 'parentFig');
   ptIndex=ptsIndices(hLines==hLine);
   hpoints=getappdata(parentFig, 'hpoints');
   hPoint=hpoints(ptIndex);
   
   if(strcmp(get(hLine, 'marker'), '*'))
      set(hLine, 'marker', 'none');
      set(hPoint, 'marker', 'o');
   else
      set(hLine, 'marker', '*');
      set(hPoint, 'marker', '+');
   end
   

function onClickOnSnip(sender, event)
   hLine=sender;
   toggleSelection(hLine);
   
   
function onDone(sender, event) 
   fig=sender;
   set(fig, 'UserData', 1);   
      
function mdexplore_polygon(hax,eventdata)
  [pvx,pvy] = GetSelPolygon('go','b');
  pv=[pvx; pvy];
  if(numel(pv)==0)
     return;
  end
  
  fig=get_parent_fig(hax);
  xproj=getappdata(fig, 'xproj');
  ptsIndices = find(inpolygon(xproj(1,:),xproj(2,:),pv(1,:),pv(2,:)));
  if(~isempty(ptsIndices))
     clusters=getappdata(fig, 'clusters');
     ptsIndices=setdiff(ptsIndices, cat(2, clusters{:}));
     if(~isempty(ptsIndices))
        clusters{end+1}=ptsIndices;
        setappdata(fig, 'clusters', clusters);
        hpoints=getappdata(fig, 'hpoints');
        color=pick_next_color(hax);
        set(hpoints(ptsIndices), 'MarkerFaceColor', color, 'MarkerEdgeColor', color);
        clusterColors=getappdata(fig, 'clusterColors');
        clusterColors{end+1}=color;
        setappdata(fig, 'clusterColors', clusterColors);
     end
  end

function mdexplore_plotcb(src,eventdata,plotdata)
   hax = get(src,'Parent');
   plotfunc = getappdata(hax,'mdplotfcn');
   if(isequal(plotfunc, @plot))
      childFig=figure;
      fig=get_parent_fig(hax);
      hChildFigs=getappdata(fig, 'hChildFigs');
      hChildFigs(end+1)=childFig;
      setappdata(fig, 'hChildFigs', hChildFigs);
   end % if, use the default plot func, use new figure too
   plotfunc(plotdata);
   axis tight


function cleanup(sender, event)
   fig=sender;
   hChildFigs=getappdata(fig, 'hChildFigs');
   free(hChildFigs);
   
   delete(fig);
   
      