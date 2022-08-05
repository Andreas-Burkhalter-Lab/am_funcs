function varargout = cass(varargin)
% CASS: Computer Aided Spike Sorting
% Type 'cass' on the command line to launch the application. You won't be
% able to do anything until you've made a selection from the File
% menu. When you start a new project, the input file you select is the
% output file written by AUTOSORT.
%
% You can affect the default behavior of CASS by supplying property/value
% pairs:
%   cass('property1',value1,...)
% where properties are from among the following:
%   dirname: the base name of the "*_sort_overview.mat" file you want to
%            process (supply only the part indicated by the *);
%   channel: a scalar indicating the channel number to open first
%
% History: 
%   ?           (TH)    wrote it
%   2007-04-12  (RCH)   added (1) ability to control where figures pop up;
%                       (2) ability to turn on "autoname" (an option that
%                       saves under the requested filename automatically
%                       every time you hit "Next Channel"; (3) the ability
%                       to turn on "auto t corr" (an option that pops up
%                       the t correlation plot automatically whenever you
%                       hit "Next Channel"
%
% See also: AUTOSORT.
  
% Copyright 2005 by Timothy E. Holy and Jason Guo
  
% Big changes afoot. The whole interaction with saved files, plus the
% possibility of the user selecting projection directions, changes many
% things.
% Here are some things I've learned:
%   -- Doing kmeans on 20000 snippets projected down to 4 dimensions is
%      still too slow (~6-8s on my laptop).  So, landmarks might need to
%      be set up before we get to the user-interactive part.  It's also
%      fair to say that it might not be great if the landmarks keep
%      changing on the user... (the user might want to be able to "back
%      up")
%   -- Any form of bootstrapping PCA is out of the question, this too
%      must be done in advance.  However, more direct PCA might be
%      possible. (Prob. not PLS-PCA.)  Update: true for Rebecca's data,
%      but with smaller d (Francesco's data), projection calculation in
%      real time is quite feasible.
%   -- mindist is _plenty_ fast to do interactively, at least in low
%      dimensions.  In fact, it's not bad even in high dimensions (0.6
%      sec for d = 100, 20000 "waveforms", 40 landmarks). It is, in fact,
%      about the same speed as projection; and after projection takes
%      place, mindist is an ignorable cost.
%      Therefore, landmarkIndex and R could be recalculated as needed.
%   -- meanshift with default options takes about half the time of
%      kmeans_hard.  This might just be do-able interactively, though not
%      likely over a range of alphaR.
% So here are changes to cass:
%   -- projectionDirections: one default for the whole file? One default
%      for the channel.  And the current pDs, which the user might
%      help define.  If nothing else, the user can cut the number of
%      dimensions down to size...
%   -- Need 2 landmarkIndices, the ones set up in advance, and other to
%      be the "current" version (perhaps re-assigned based on the
%      current projection directions)
%   -- However, landmark waveforms are fixed forever...perhaps those
%      should be saved by autosort?  If we are allowed to re-assign
%      landmarkIndices, then this appears to be necessary, since only one
%      set of landmarkIndices is saved to the file (per replicate)---by
%      saving, loading, saving, loading, the user would otherwise end up
%      doing iterations of kmeans!
%   -- In fact, perhaps the spike times and landmarkIndex simply should
%      not be saved to the sort_subset file?  That might simplify matters
%      a great deal!
%   -- The idea of replicates might eventually go away? Just use plenty
%      of landmarks. But let's not eliminate that possibility in the file
%      format, even if the cass GUI ends up having it eliminated...
%   -- Maintain the following data structures:
%         file_sort_info: all the sort_info in a given file
%         current_sort_info: the current choices of waveforms/clusters,
%           etc. Note this will have an additional field: landmarkDelete,
%           which contains landmark waveforms which should be "really
%           trashed", e.g. artifacts.

%CASS M-file for cass.fig
%      CASS, by itself, creates a new CASS or raises the existing
%      singleton*.
%
%      H = CASS returns the handle to a new CASS or the handle to
%      the existing singleton*.
%
%      CASS('Property','Value',...) creates a new CASS using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to cass_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CASS('CALLBACK') and CASS('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CASS.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES


% General bug history: 
%     * bug in file name selection in choose_sort_info_file fixed (RCH)
%     * hack job used to allow empty cells in peakHeight to be transformed
%       into peakHeightall without cat choking (RCH, 2005-10-30)
%     * (1) if you try to add start/stop cell timeMarkers when sorting
%       multielectrode data will just treat like start/stop all instead of
%       choking; (2) turns on the "modified" flag when a timeMarker is
%       added so will save if that's all you've changed when you advance
%       channels (RCH, 2007-03-28)

% Last Modified by GUIDE v2.5 30-Nov-2007 15:52:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cass_OpeningFcn, ...
                   'gui_OutputFcn',  @cass_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cass is made visible.
function cass_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

  % Choose default command line output for cass
  handles.output = hObject;

  % Set default values, letting the user supply command line
  % property/value pairs
  options = cass_options;
  optionsFieldnames = fieldnames(options);
  if ~isempty(varargin)
    for i = 1:2:length(varargin)
      if isempty(strmatch(varargin{i},optionsFieldnames,'exact'))
        error(['Property ' varargin{i} ' not recognized'])
      end
      options.(varargin{i}) = varargin{i+1};
    end
  end
  
  % Set up convenient access to cluster display handles
  objectTags = fieldnames(handles);   % determine how many there are
  axesW_index = strmatch('axesW',objectTags);
  nClusterDisplays = length(axesW_index);
  axes_IDs=1:nClusterDisplays;
  handle_matrix=nan(3, nClusterDisplays); % #rows: 2 row axes + 1 row text
  for axes_ID=axes_IDs
    handle_matrix(1,axes_ID)=handles.(['axesW' num2str(axes_ID)]); % Waveform
    handle_matrix(2,axes_ID)=handles.(['axesPh' num2str(axes_ID)]); % Peak Histogram
    handle_matrix(3,axes_ID)=handles.(['textClust' num2str(axes_ID)]);
  end
   
  handles.handle_matrix=handle_matrix;
  handles.isSplitting = false;    % State of splitting mode
  
  % Update handles structure
  guidata(hObject, handles);
  
  % Set up data to determine cluster selection and display order
  % We have only 7 axes to use for display, but because we can scroll we
  % can have a virtual number as high as we want. If we want to let the
  % user change the order in which clusters are displayed, then we need a
  % lookup, virtual_axis(clusternumber).  Because we'll use
  % clusternumber=0 for noise, let's actually do
  % virtual_axis(clusternumber+1). And, let's keep track of selections in
  % the same structure.
  axes_map=struct('virtual_axis',  num2cell([1]), ...
                  'selected', {'off'} ...
                  ); % 0 for noise, a positive number for each cluster
  % Because of scrolling, we have to keep track of current position.
  % This is in terms of the virtual_axis #, not the cluster #.
  axes_map_start_from=1;
  setappdata(handles.figCass, 'axes_map', axes_map);
  setappdata(handles.figCass, 'axes_map_start_from', axes_map_start_from);

  setappdata(handles.figCass, 'options', options);

  cass_start_project(handles);
  
  % Install mouse handler for cluster axes
  for axes_ID=axes_IDs
    tHandle=handles.handle_matrix(1, axes_ID);
    setappdata(tHandle, 'isDragging', 0);
    install_mouse_event_handler(tHandle, 'up', @axes_mouse_up);      
    install_mouse_event_handler(tHandle, 'down', @axes_mouse_down);      
    install_mouse_event_handler(tHandle, 'move', @axes_mouse_move);
    % Install context menu on cluster plot
    hcontext = uicontextmenu('parent',get(tHandle,'parent'));
    uimenu(hcontext,'Label','Mark for further processing',...
      'Callback',{@set_further_processing,handles,tHandle,'mark_cluster'});
    set(tHandle,'UIContextMenu',hcontext);
  end
  % Install handler for zoom on time axis & waveform projection axis
  zoomHandle = [handles.axesTime handles.axesVisualizeShape];
  for tHandle = zoomHandle
    install_mouse_event_handler(tHandle, 'up', @zoom_mouse_up);
    install_mouse_event_handler(tHandle, 'down', @zoom_mouse_down);
    install_mouse_event_handler(tHandle, 'move', @zoom_mouse_move);
    setappdata(tHandle, 'isDragging', 0);
  end
  % Install context menu on time axis
  hcontext = uicontextmenu('Parent',handles.figCass);
  set(handles.axesTime,'UIContextMenu',hcontext);
  uimenu(hcontext,'Label','Mark cell',...
         'Callback',{@set_time_marker,handles,'mark_cell'});
  uimenu(hcontext,'Label','Mark all',...
         'Callback',{@set_time_marker,handles,'mark_all'});
  uimenu(hcontext,'Label','Start cell',...
         'Callback',{@set_time_marker,handles,'start_cell'});
  uimenu(hcontext,'Label','Stop cell',...
         'Callback',{@set_time_marker,handles,'stop_cell'});
  uimenu(hcontext,'Label','Start all',...
         'Callback',{@set_time_marker,handles,'start_all'});
  uimenu(hcontext,'Label','Stop all',...
         'Callback',{@set_time_marker,handles,'stop_all'});
  %handles.uicontextTime = hcontext;
  %guidata(hObject,handles);
  % Install handler on figure for unselect
  install_mouse_event_handler(handles.figCass, 'up', @figure_mouse_up);

  % UIWAIT makes cass wait for user response (see UIRESUME)
  % uiwait(handles.figure1);


%
% Functions for manipulating cluster axes
%
function result=axes_mouse_down(sender, event_args)
   setappdata(sender, 'isMouseDown', 1);
   result=0;
   
function result=axes_mouse_move(sender, event_args)
   if(getappdata(sender, 'isMouseDown'))
      setappdata(sender, 'isDragging', 1);
   end
   result=0;
   
function result=axes_mouse_up(sender, event_args)
  if(~getappdata(sender, 'isDragging'))
    handles=guidata(sender);
    axes_map=getappdata(handles.figCass, 'axes_map');
    axes_map_start_from=getappdata(handles.figCass, 'axes_map_start_from');
    virtual_axis = [axes_map.virtual_axis];
    n_visible_axes = size(handles.handle_matrix,2);
    clusterIndex = find(virtual_axis >= axes_map_start_from & ...
                        virtual_axis < axes_map_start_from + ...
                        n_visible_axes);
    thisAxisIndex = find(handles.handle_matrix(1,:) == sender);
    clIndex = clusterIndex(thisAxisIndex); % Here's the cluster in axes_map
    update_selections(handles,clIndex,sender);
%   else % is dragging: that is, re-order bins
%     handles=guidata(sender);
%     axesDroppedOn=[];
%     axes_IDs=1:size(handles.handle_matrix,2);
%     for axes_ID=axes_IDs
%       tHandle=handles.handle_matrix(1,axes_ID);
%       if(is_mouse_over(tHandle))
%         axesDroppedOn=tHandle; break;
%       end
%     end
%     if(~isempty(axesDroppedOn))
%       tagFrom=get(sender, 'tag');
%       tagTo=get(axesDroppedOn, 'tag');
%       idFrom=str2num(tagFrom(end)); idTo=str2num(tagTo(end));
%       tDiff=idTo-idFrom;
%       axes_map  =getappdata(handles.figCass, 'axes_map');
%       axes_map_start_from =getappdata(handles.figCass, 'axes_map_start_from');
%       start_idx=idFrom+axes_map_start_from-1;
%       end_idx=idTo+axes_map_start_from-1;
%       tt=axes_map(start_idx);
%       axes_map(start_idx)=[];
%       axes_map=[axes_map(1:end_idx-1) tt axes_map(end_idx:end)];
%       setappdata(handles.figCass, 'axes_map', axes_map);
%     end
%   end
  end
  setappdata(sender, 'isMouseDown', 0);
  setappdata(sender, 'isDragging', 0);
   
  result=1;

  
%
% Functions for zooming on time & waveform projection axes; also for
% placing markers on the time axis, and selecting clusters on the
% projection axis.
%
% On the projection axis:
%   Left mouse drag: zoom in
%   Left mouse click: zoom all the way out
%   Right mouse: select (hold shift for adding---same as a middle click)
% On the time axis:
%   Left mouse drag: zoom in
%   Right mouse: add breakpoints
%
function result=zoom_mouse_down(sender, event_args)
  setappdata(sender, 'isMouseDown', 1);
  currentPoint = get(sender,'CurrentPoint');
  setappdata(sender,'startPoint',currentPoint(1,[1 2]));
  result=0;
   
function result=zoom_mouse_move(sender, event_args)
  handles=guidata(sender);
  if(getappdata(sender, 'isMouseDown'))
    setappdata(sender, 'isDragging', 1);
    currentPoint = get(sender,'CurrentPoint');
    currentPoint = currentPoint(1,[1 2]);
    startPoint = getappdata(sender,'startPoint');
    if handles.isSplitting
      xdata = [startPoint(1) currentPoint(1)];
      ydata = [startPoint(2) currentPoint(2)];
    else
      xdata = [startPoint(1) currentPoint([1 1]) startPoint([1 1])];
      ydata = [startPoint([2 2]) currentPoint([2 2]) startPoint(2)];
    end
    if ~isappdata(sender,'hrbbox')
      hrbbox = line(xdata,ydata,...
                    'LineStyle','--',...
                    'Color','k',...
                    'EraseMode','xor',...
                    'Tag','rbbox');
      setappdata(sender,'hrbbox',hrbbox);
    else
      hrbbox = getappdata(sender,'hrbbox');
      set(hrbbox,'XData',xdata,'YData',ydata);
    end
  end
  result=0;
   
function result=zoom_mouse_up(sender, event_args)
  handles=guidata(sender);
  if(~getappdata(sender, 'isDragging') && sender == handles.axesVisualizeShape)
    % Projection axis: this is either a zoom out (left click) or a
    % cluster selection (right click or shift-click/middle click)
    % Time axis: shouldn't be here (right click all handled through
    % context menu, and zooming operates in a new figure)
    startPoint = getappdata(sender,'startPoint');
    selection_type = lower(get(get_parent_fig(sender), 'SelectionType'));
    % Clean up mouse stuff now, to avoid weird behavior if this gets
    % interrupted
    setappdata(sender, 'isMouseDown', 0);
    setappdata(sender, 'isDragging', 0);
    result=1;
    if strcmp(selection_type,'normal')
      oxlim = getappdata(sender,'OXLim');
      oylim = getappdata(sender,'OYLim');
      set(sender,'XLim',oxlim,'YLim',oylim);
    elseif ~strcmp(selection_type,'open')
      % This is a cluster selection; find the closest point
      snipProj = getappdata(handles.figCass,'snipProj');
      x_coord = get(handles.popupXVis,'Value');
      y_coord = get(handles.popupYVis,'Value');
      [dist,closestIndex] = mindist(startPoint([1 2])',snipProj([x_coord ...
                          y_coord],:));
      % Determine the cluster identity of the closest point
      landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
      current_sort_info = getappdata(handles.figCass,'current_sort_info');
      clIndex = ...
          current_sort_info.landmarkClust(landmarkIndex(closestIndex))+1;
      % Determine whether this cluster is showing on the cluster axes,
      % and, if so, get the axis handle.
      axes_map = getappdata(handles.figCass,'axes_map');
      axes_map_start_from = getappdata(handles.figCass, ...
                                   'axes_map_start_from');
      virtual_axis = [axes_map.virtual_axis];
      n_visible_axes = size(handles.handle_matrix,2);
      clusterIndex = find(virtual_axis >= axes_map_start_from & ...
                          virtual_axis < axes_map_start_from + n_visible_axes);
      handleIndex = find(clusterIndex == clIndex);
      axHandle = [];
      if ~isempty(handleIndex)
        axHandle = handles.handle_matrix(1,handleIndex);
      end
      update_selections(handles,clIndex,axHandle);
    end      
% $$$     % This is either a zoom in around a point (left click), or a zoom out
% $$$     % (right click), or a "restore to original settings" (double click),
% $$$     % or it's a middle click on the time axis.
% $$$     handles=guidata(sender);
% $$$     startPoint = getappdata(sender,'startPoint');
% $$$     selection_type = lower(get(get_parent_fig(sender), 'SelectionType'))
% $$$     % Clean up mouse stuff now, to avoid weird behavior if this gets
% $$$     % interrupted
% $$$     setappdata(sender, 'isMouseDown', 0);
% $$$     setappdata(sender, 'isDragging', 0);
% $$$     result=1;
% $$$     if (sender == handles.axesVisualizeShape)
% $$$       xlim = get(sender,'XLim');
% $$$       ylim = get(sender,'YLim');
% $$$       dx = diff(xlim);
% $$$       dy = diff(ylim);
% $$$       oxlim = getappdata(sender,'OXLim');
% $$$       oylim = getappdata(sender,'OYLim');
% $$$       switch selection_type
% $$$        case 'normal'
% $$$         % zoom in
% $$$         xlim = startPoint(1) + dx*[-0.25 0.25];
% $$$         ylim = startPoint(2) + dy*[-0.25 0.25];
% $$$         set(sender,'XLim',xlim,'YLim',ylim);
% $$$        case 'alt'
% $$$         % zoom out
% $$$         xlim = startPoint(1) + dx*[-1 1];
% $$$         ylim = startPoint(2) + dy*[-1 1];
% $$$         xlim = IntersectIntervals(xlim,oxlim);
% $$$         ylim = IntersectIntervals(ylim,oylim);
% $$$         set(sender,'XLim',xlim,'YLim',ylim);
% $$$        case 'open'
% $$$         % restore to original settings
% $$$         set(sender,'XLim',oxlim,'YLim',oylim);
% $$$       end
% $$$     else
% $$$       % We're on the time axis, so set a marker
% $$$       if strcmp(selection_type,'alternate')
% $$$         set(handles.uicontextTime,'Visible','on');
% $$$       end
% $$$     end
  else
    % We're dragging. This either defines a zoom region, or a splitting
    % line. For zoom, if this is the waveform 
    % projection axis, we can take the limits literally, but if it's the
    % time axis then we create a new figure for reconstruction
    currentPoint = get(sender,'CurrentPoint');
    currentPoint = currentPoint(1,[1 2]);
    startPoint = getappdata(sender,'startPoint');
    % Clean up the mouse stuff, since some of the following code takes a
    % while to execute, and pauses can introduce some weird behavior
    setappdata(sender, 'isMouseDown', 0);
    setappdata(sender, 'isDragging', 0);
    if isappdata(sender,'hrbbox')
      hrbbox = getappdata(sender,'hrbbox');
      delete(hrbbox);
      rmappdata(sender,'hrbbox');
    end
    result=1;
    if (sender == handles.axesVisualizeShape)
      if handles.isSplitting
        cass_split_selected_cluster(handles,startPoint,currentPoint);
      else
        xlim = sort([startPoint(1) currentPoint(1)]);
        ylim = sort([startPoint(2) currentPoint(2)]);
        oxlim = getappdata(sender,'OXLim');
        oylim = getappdata(sender,'OYLim');
        xlim = IntersectIntervals(xlim,oxlim);
        ylim = IntersectIntervals(ylim,oylim);
        if (ylim(1)~=ylim(2)) & (xlim(1)~=xlim(2))
          set(sender,'XLim',xlim,'YLim',ylim);
        end
      end
    else
      % We're on the time axis
      selection_type = lower(get(get_parent_fig(sender), 'SelectionType'));
      if strcmp(selection_type,'normal')
        tlim = sort([startPoint(1) currentPoint(1)]);
        handles = guidata(sender);
        cass_waveform_reconstruct(handles,tlim)
      end
    end
  end
  
%
% Functions for manipulating time markers
%
function marker_dragged(src,eventData,info)
  xpos = get(gco,'XData');
  timeMarker = getappdata(info.dataHandle,'timeMarker');
  index = find([timeMarker.relativeTime] == info.startpos);
  timeMarker(index).relativeTime = xpos(1);
  setappdata(info.dataHandle,'timeMarker',timeMarker);
  % Update the startpos info so that it can match if it gets dragged again
  bdf = get(gco,'ButtonDownFcn');
  bdf{2}.doneargs{1}.startpos = xpos(1);
  set(gco,'ButtonDownFcn',bdf);
  if ~isequal(gcbf,info.dataHandle)
    handles = guidata(info.dataHandle);
    update_time_axis_display(handles)
  end
  
function delete_time_marker(src,eventData)
  xpos = get(gco,'XData');
  timeMarker = getappdata(gcbf,'timeMarker');
  index = find([timeMarker.relativeTime] == xpos(1));
  timeMarker(index) = [];
  setappdata(gcbf,'timeMarker',timeMarker)
  delete(gco)

function edit_time_marker_comment(src,eventData)
  xpos = get(gco,'XData');
  hfig = get_parent_fig(gco);
  timeMarker = getappdata(hfig,'timeMarker');
  index = find([timeMarker.relativeTime] == xpos(1));
  if (~isfield(timeMarker,'comment') || isempty(timeMarker(index).comment))
    timeMarker(index).comment = {''};
  end
  timeMarker(index).comment = inputdlg('Describe this marker',...
                                       'Time marker',3,...
                                       timeMarker(index).comment);
  setappdata(hfig,'timeMarker',timeMarker)
  
  
function timeMarker = convert_timeMarker_to_relativeTime(timeMarker,handles)
  tstart = getappdata(handles.figCass,'absoluteFileStartTime');
  tstart = tstart-min(tstart);
  sorthead = getappdata(handles.figCass,'sorthead');
  if ~isfield(timeMarker,'relativeTime')
    for i = 1:length(timeMarker)
      fileIndex = timeMarker(i).fileIndex;
      timeMarker(i).relativeTime = tstart(fileIndex) ...
        + timeMarker(i).fileTime/sorthead(fileIndex).scanrate;
    end
  end
  
function timeMarker = convert_timeMarker_to_fileTime(timeMarker,handles)
  tstart = getappdata(handles.figCass,'absoluteFileStartTime');
  tstart = tstart - min(tstart);
  tlength = getappdata(handles.figCass,'fileDuration');
  sorthead = getappdata(handles.figCass,'sorthead');
  for i = 1:length(timeMarker)
    if ~isfield(timeMarker,'fileIndex') | isempty(timeMarker(i).fileIndex)
      fileIndex = find(tstart < timeMarker(i).relativeTime,1,'last');
      timeMarker(i).fileIndex = fileIndex;
      timeMarker(i).fileTime = ...
        round((timeMarker(i).relativeTime-tstart(fileIndex))...
        *sorthead(fileIndex).scanrate);
    end
  end
      
%
% Function for unselecting cluster axes (click in figure)
%
function result=figure_mouse_up(sender, event_args)
  handles = guidata(sender);
  set(handles.handle_matrix,'selected','off')
  axes_map = getappdata(handles.figCass,'axes_map');
  for i = 1:length(axes_map)
    axes_map(i).selected = 'off';
  end
  setappdata(handles.figCass,'axes_map',axes_map);
  % Update marker sizes, too
  plH = [axes_map.projLineHandle];
  markersize = get(plH,'MarkerSize');
  if iscell(markersize)
    markersize = cell2mat(markersize);
  end
  indexFat = find(markersize > 1);
  set(plH(indexFat),'MarkerSize',1);
  result = 0;


%
% Functions called once per sorting project
%
function cass_parse_overview(handles)
  options = getappdata(handles.figCass,'options');
  if isempty(options.dirname)
    return
  end
  clear_overview_data(handles); % eliminate any previous settings, if this
                                % function fails (e.g., user clicks cancel)
  sortoverview = load([options.dirname filesep 'overview']);
  
  % check if the autosort result is fake or not
  if(isfield(sortoverview, 'isFake'))
     isFake=sortoverview.isFake;
  else
     isFake=0;
  end
  setappdata(handles.figCass, 'isFake', isFake);
  
  % load fake info if isFake
  if(isFake)
     % TODO:   
     
  end % if, isFake
  
  sorthead = sortoverview.sorthead;
  options_autosort = sortoverview.options;
  options.channels = sorted_get_channels(options.dirname);
  startIndex = find(options.channels == options.channel);
  if isempty(startIndex)
    startIndex = 1;
  end
  set(handles.popupChannel,...
      'String',num2str(options.channels'),...
      'Value',startIndex);
  % Make sure we can find the individual snippet files; prompt user for
  % help if they can't be found
  [sorthead,pathChanged] = cass_fix_path(sorthead,struct('first_try_pathName',[options.dirname filesep]));
  setappdata(handles.figCass,'sorthead',sorthead); % store path changes
  % If we had to ask the user for assistance, let's give him/her a chance
  % to store the path change in the sort_overview file, so that we don't
  % have to ask for help if cass is run again later
  if pathChanged
    button = questdlg(['Do you want me to save the path changes in the ' ...
                       options.dirname '_sort_overview file?'],...
                      'Path changes','Yes','No','Yes');
    if strcmp(button,'Yes')
      savestruct.sorthead = sorthead;
      savestruct.options = options_autosort;
      save([options.dirname 'overview'],'-struct',savestruct);
    end
  end
  % Store data about file timing
  tstart = sortheader_absolute_starttime(sorthead); % File start time
  tlen = [sorthead.nscans]./[sorthead.scanrate];
  setappdata(handles.figCass,'absoluteFileStartTime',tstart);
  setappdata(handles.figCass,'fileDuration',tlen);
  setappdata(handles.figCass,'options',options);


function clear_overview_data(handles)
  fields_to_clear = {'sorthead','absoluteFileStartTime','fileDuration'};
  for i = 1:length(fields_to_clear)
    if isappdata(handles.figCass,fields_to_clear{i})
      rmappdata(handles.figCass,fields_to_clear{i})
    end
  end


  
%
% Functions called when a channel is started, replicate/alpha chosen, or
% reset
%
function filename = cass_choose_sort_info_file(handles)
  options = getappdata(handles.figCass,'options');
  if isempty(options.dirname)
    filename = '';
    return;
  end
  currentChannel = getappdata(handles.figCass,'currentChannel');
  dirchan = [options.dirname filesep 'chan' num2str(currentChannel) ...
             filesep];
  filename = choose_sort_info_file(dirchan,options.force_autosort_info);
  %filename = [dirchan filename];

  
function cass_load_sort_info(handles,filename)
  if isempty(filename)
    return
  end
  options = getappdata(handles.figCass,'options');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  dirchan = [options.dirname filesep 'chan' num2str(currentChannel) ...
             filesep];
  sort_info = load_sort_info([dirchan filename]);
  if (length(sort_info) > 1)
    % Enable replicates
    set(handles.popupReplicate,...
        'String',num2str((1:length(sort_info))'),...
        'Value',1,...
        'Enable','on');
  else
    set(handles.popupReplicate,'Enable','off');
  end
  if isfield(sort_info,'Rfactor') && ~isequal(sort_info(1).Rfactor,nan)
    % Enable alphaR (& choose alphaR closest to default)
    [tmp,alphaRIndex] = min(abs(sort_info(1).Rfactor - ...
                              options.RfactorDefault));
    set(handles.popupAlphaR,...
        'String',num2str(sort_info(1).Rfactor'),...
        'Value',alphaRIndex,...
        'Enable','on');
  else
    set(handles.popupAlphaR,'String','--','Enable','off');
  end
  setappdata(handles.figCass,'sort_info',sort_info);
  

function [current_sort_info,timeMarker] = select_current_sort_info(handles)
% Select sorting results, perhaps based on popup values.
  current_sort_info = getappdata(handles.figCass,'sort_info');
  if isempty(current_sort_info)
    return
  end
  if (length(current_sort_info) > 1)
    replicate = get(handles.popupReplicate,'Value');
    current_sort_info = current_sort_info(replicate);
  end
  if isfield(current_sort_info,'Rfactor') && length(current_sort_info.Rfactor) > 1
    alphaRIndex = get(handles.popupAlphaR,'Value');
    current_sort_info.Rfactor = current_sort_info.Rfactor(alphaRIndex);
    current_sort_info.landmarkClust = ...
        current_sort_info.landmarkClust(alphaRIndex,:);
  end
  if isfield(current_sort_info,'timeMarker')
    % Have to convert to relativetime
    timeMarker = current_sort_info.timeMarker;
    timeMarker = convert_timeMarker_to_relativeTime(timeMarker,handles);
  else
    timeMarker = [];
  end


function store_current_sort_info(handles,current_sort_info,timeMarker)
  if isempty(current_sort_info)
    return
  end
  setappdata(handles.figCass,'current_sort_info',current_sort_info);
  setappdata(handles.figCass,'timeMarker', timeMarker);
  if isfield(current_sort_info,'comment')
    set(handles.editComment,'String',current_sort_info.comment);
  end
  setappdata(handles.figCass,'modified',0);
  

function set_current_sort_info(handles)
% Set/reset sorting results, perhaps based on popup values.
  [current_sort_info,timeMarker] = select_current_sort_info(handles);
  store_current_sort_info(handles,current_sort_info,timeMarker);
  




function cass_load_snippets(handles)
  if ~isappdata(handles.figCass,'sorthead')
    return
  end
  % Get the file data
  sorthead = getappdata(handles.figCass,'sorthead');
  options = getappdata(handles.figCass,'options');
  % Find the selected channel & import it
  currentChannel = getappdata(handles.figCass,'currentChannel');
  shc = sortheader_importchan(sorthead,currentChannel);
  sort_info = getappdata(handles.figCass,'sort_info');
  % Pick the snippets we're going to load (i.e., choose snipIndex)
  nFiles = length(shc);
  nsnips_per_file = [shc.numofsnips];
  nsnips_total = sum(nsnips_per_file);
  sniprange = reshape([shc.sniprange],2,length(shc));
  sniprange = unique(sniprange','rows');
  if (size(sniprange,1) > 1)
    error(['Snippets do not have the same number of samples in each ' ...
           'file']);
  end
  % The calculation of the # of snippets to sort changed once we went to
  % avoiding storing whole snippets in memory
  if ~sort_info(1).use_projection
    storage_per_snip = diff(sniprange) + 6;
  else
    storage_per_snip = size(sort_info(1).projectDirections,2) + 6;
    % the +6 is an attempt to account for sniptime, landmarkIndex, minmax,
    % etc.
  end
  max_snips_in_mem = min(round(options.max_snip_memsize/storage_per_snip),...
    options.max_snips_to_cluster);
  subset_fraction = min(1,max_snips_in_mem/nsnips_total);
  stride = ceil(1/subset_fraction);
  snipIndex = cell(1,nFiles);
  nsnips_loaded = nan(1,nFiles);
  for fileIndex = 1:nFiles
    % Select evenly-spaced spikes
    snipIndex{fileIndex} = 1:stride:nsnips_per_file(fileIndex);
    nsnips_loaded(fileIndex) = length(snipIndex{fileIndex});
  end
  % Update the %spikes display
  percentSpikes = round(100*sum(nsnips_loaded)/nsnips_total);
  set(handles.txtPercentSpikes,'String',[num2str(percentSpikes) '%']);
  % Load the data
  tstart = getappdata(handles.figCass,'absoluteFileStartTime');
  peakHeight = cell(1,length(sorthead));
  absoluteSnipTime = cell(1,length(sorthead));
  if isfield(shc,'fc_rank')
    snipRank = cell(1,length(sorthead));
    maxRank = 0;
    for i = 1:nFiles
      maxRank = max(maxRank,max(shc(i).fc_rank));
    end
    set(handles.popupEventRank,...
      'Visible','on',...
      'String',num2str((1:maxRank)'),...
      'Value',maxRank);
    set(handles.textEventRank,...
      'Visible','on');
  else
    set([handles.popupEventRank handles.textEventRank],...
      'Visible','off');
  end
  %snip = cell(1,length(sorthead));
  for i = 1:nFiles
    % Grab peak heights
    peakHeight{i} = double(shc(i).peakHeight(snipIndex{i})');
    % Convert peakHeight to volts
    peakHeight{i} = peakHeight{i}*sorthead(i).scalemult + ...
        sorthead(i).scaleoff;
    % Get the snip times
    % For user feedback, convert sniptimes to absolute timescale
    sniptimes = double(shc(i).sniptimes(snipIndex{i})');
    absoluteSnipTime{i} = sniptimes/sorthead(i).scanrate ...
        + tstart(i) - min(tstart);
    if isfield(shc,'fc_rank')
      snipRank{i} = shc(i).fc_rank(snipIndex{i});
    end
    % Read in snippets
    %snip(i) = sortheader_readsnips(shc(i),snipIndex(i));
  end
  % Check to see if min/max and projections have been pre-calculated
  projfilename = [options.dirname filesep 'chan' num2str(currentChannel) ...
             filesep 'proj.mat'];
  precalc = false;
  if exist(projfilename,'file')
    precalc = true;
    stmp = load(projfilename);
    for fileIndex = 1:nFiles
      snipminmax{fileIndex} = stmp.snipmm{fileIndex}(:,snipIndex{fileIndex});
      snipProj{fileIndex} = stmp.snipProj{fileIndex}(:,snipIndex{fileIndex});
    end
  else
    % Read in snippet min/max pairs
    snipminmax = sortheader_snipminmax(shc,snipIndex);
  end
  % Calculate binning parameters for peak height
  % We only have to do this once when a channel is first loaded in; these
  % parameters are independent of clustering
  % Check that polarity is consistent across files
  polarity = unique([sorthead.polarity]);
  if (length(polarity) > 1)
    errordlg('Polarity not consistent across files');
  end
  peak_height_data.polarity = polarity;
  % Determine whether the recorded voltage range is consistent across
  % files
  if isfield(sorthead,'voltageMax')
    voltageMax = unique([sorthead.voltageMax]);
    voltageMin = unique([sorthead.voltageMin]);
    if (length(voltageMax) > 1 || length(voltageMin) > 1)
      warndlg('Voltage range was not consistent across files');
    end
    peak_height_data.voltageRange = [max(voltageMin) min(voltageMax)];
  end
  % Calculate thresholds for this channel, making sure they're consistent
  % across files
  for i = 1:length(sorthead)
    channelIndex = find(sorthead(i).channels == currentChannel);
    thresh(:,i) = sorthead(i).thresh(:,channelIndex)*sorthead(i).scalemult ...
        + sorthead(i).scaleoff;
  end
  uthresh = unique(thresh','rows');
  if (size(uthresh,1) > 1)
    warndlg('Different snippet files use different voltage thresholds!')
    % Choose the most extreme threshold(s)
    switch polarity
     case -1
      uthresh = min(uthresh);
     case 1
      uthresh = max(uthresh);
     case 0
      uthresh = [min(uthresh(:,1)) max(uthresh(:,2))];
    end
  end
  peak_height_data.thresh = uthresh;
  % Calculate binning parameters for plotting peakHeight distribution
  n = length(peakHeight);
  for nth = 1:n
    if size(peakHeight{nth},1) == 0
      peakHeight{nth} = (peakHeight{nth})';
    end
  end
  peakHeightall = cat(2,peakHeight{:});
  minpH = min(peakHeightall);
  maxpH = max(peakHeightall);
  nBins = min(ceil(sqrt(length(peakHeightall))),50);
  peak_height_data.binEdges = linspace(minpH,maxpH*(1+eps),nBins);
  % Save data in figure
  setappdata(handles.figCass,'peak_height_data',peak_height_data);
  setappdata(handles.figCass,'peakHeight',peakHeight);
  setappdata(handles.figCass,'snipIndex',snipIndex);
  setappdata(handles.figCass,'snipIndexStride',stride);
  setappdata(handles.figCass,'absoluteSnipTime',absoluteSnipTime);
  %setappdata(handles.figCass,'snip',snip);
  setappdata(handles.figCass,'nsnips_per_file',nsnips_per_file);
  setappdata(handles.figCass,'snipminmax',snipminmax);
  setappdata(handles.figCass,'precalc',precalc);
  if precalc
    setappdata(handles.figCass,'snipProj0',snipProj);
  end
  if isfield(shc,'fc_rank')
    setappdata(handles.figCass,'snipRank',snipRank);
  end

function clear_channel_data(handles)
  fields_to_clear = {'snipIndex','peakHeight','absoluteSnipTime','snip',...
                     'snipminmax',...
                     'peak_height_data','snipIndexStride',...
                     'nsnips_per_file','landmarkIndex', ...
                     'sort_info','current_sort_info',...
                     'saveFile',...
                     'modified',...
                     'landmarkIndex_total', ...
                     'cluster_axes_data',...
                     'snipProj','landmarkWaveformProj',...
                     'alphaRIndexStable','stabilityValid'};
  for i = 1:length(fields_to_clear)
    if isappdata(handles.figCass,fields_to_clear{i})
      rmappdata(handles.figCass,fields_to_clear{i})
    end
  end
  % Put the cluster panner in the right state
  setappdata(handles.figCass,'axes_map_start_from',1)
  % Reset zoom limits
  if isappdata(handles.axesVisualizeShape,'OXLim')
    rmappdata(handles.axesVisualizeShape,'OXLim')
  end
  if isappdata(handles.axesVisualizeShape,'OYLim')
    rmappdata(handles.axesVisualizeShape,'OYLim')
  end
  % Clear the edit text
  set(handles.editComment,'String',{''});
  

  
%
% Functions called when projections are reset, or when a channel is started
%
function set_projection_pulldowns(handles)
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  options = getappdata(handles.figCass,'options');
  % Even if not using projection to sort the data, will use it for
  % visualization
  n_spatial_directions = size(current_sort_info.projectDirections,2);
  strPulldown = num2str((1:n_spatial_directions)');
  strPulldown(end+1,1) = 't';  % add time as a direction
  if isfield(options,'proj_default_xaxis')
    if ischar(options.proj_default_xaxis)
      if strmatch(options.proj_default_xaxis,'t')
        set(handles.popupXVis,'String',strPulldown,'Value', ...
                          size(strPulldown,1));
      end
    else
      set(handles.popupXVis,'String',strPulldown,'Value', ...
                        options.proj_default_axis);
    end
  else
    set(handles.popupXVis,'String',strPulldown,'Value',2);
  end
  set(handles.popupYVis,'String',strPulldown,'Value',1);
  if isappdata(handles.axesVisualizeShape,'OXLim')
    rmappdata(handles.axesVisualizeShape,'OXLim')
  end
  if isappdata(handles.axesVisualizeShape,'OYLim')
    rmappdata(handles.axesVisualizeShape,'OYLim')
  end
  
%
% Recalculation functions
%
function update_calculations(handles)
  if ~isappdata(handles.figCass,'sort_info')
    return
  end
  set_current_sort_info(handles)
  set_projection_pulldowns(handles)
  update_calculations_rest(handles)
  
function update_calculations_rest(handles)
  update_projections(handles)
  update_landmarkIndex(handles)
  update_cluster_axes_calculations(handles)
  update_alphaR_range_calculations(handles)
  

function update_projections(handles)
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  if isempty(current_sort_info)
    return;
  end
  % Get the snippets
  %snip = getappdata(handles.figCass,'snip');
  snipIndex = getappdata(handles.figCass,'snipIndex');
  sorthead = getappdata(handles.figCass,'sorthead');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  isFake=getappdata(handles.figCass, 'isFake');
  shc = sortheader_importchan(sorthead,currentChannel);
  precalc = getappdata(handles.figCass,'precalc');
  if ~current_sort_info.use_projection
    % If we're not using projection, then we have to load them all
    snip = sortheader_readsnips(shc,snipIndex);
    snipProj = cat(2,snip{:});
    if(isFake)
       error('not implemented for fake data');
    end
  elseif precalc
    % We don't have to do much, except calculate snipIndicesUsed
    n_files = length(snipIndex);
    for fileIndex = 1:n_files
      snipIndicesUsed{fileIndex} = ...
	  [repmat(fileIndex,1,length(snipIndex{fileIndex}));
	   snipIndex{fileIndex}];
    end
    snipIndicesUsed = cat(2,snipIndicesUsed{:});
    snipProj = getappdata(handles.figCass,'snipProj0');
    snipProj = cat(2,snipProj{:});
  else
    % We're using projection.  We should collect the snippets in blocks,
    % in case we don't have the memory to store them all at once
    blocksize = 5000;  % process in blocks of 5000 snippets
    snipProj = {};
    % Create subsets of snipIndex. It might be easier to do this directly
    % in terms of the disk representation, but then we lose our
    % abstraction through the sortheader_readsnips function.
    n_files = length(snipIndex);
    for fileIndex = 1:n_files
      n_snips_per_file(fileIndex) = length(snipIndex{fileIndex});
    end
    cum_n_snips_per_file = [0 cumsum(n_snips_per_file)];
    offset = 1;
    
    snipIndicesUsed=zeros(2,0); % the used snip indices, top row is file indices
    
    while (offset < cum_n_snips_per_file(end))
      snipIndexTmp = cell(1,n_files);
      for fileIndex = 1:n_files
        rangemin = max(1,offset-cum_n_snips_per_file(fileIndex));
        rangemax = min(n_snips_per_file(fileIndex),...
          offset+blocksize-1-cum_n_snips_per_file(fileIndex));
        snipIndexTmp{fileIndex} = snipIndex{fileIndex}(rangemin: ...
          rangemax);
       
        snipIndicesUsed=[snipIndicesUsed ...
           [fileIndex*ones(1, length(snipIndexTmp{fileIndex})); snipIndexTmp{fileIndex}] ];
      end
      sniptmp = sortheader_readsnips(shc,snipIndexTmp);
      sniptmp = cat(2,sniptmp{:}); % cat sniptmp from different files
      snipProj{end+1} = current_sort_info.projectDirections'*sniptmp;
      offset = offset+blocksize;
    end
    snipProj = cat(2,snipProj{:}); % cat snipProj from different blocks
  end
  
  % actually all snippets are used but snipIndex can't be used b/c snippets
  % are readed block by block
  if(isFake)
     setappdata(handles.figCass, 'snipIndicesUsed', snipIndicesUsed);
     setappdata(handles.figCass, 'shc', shc);
  end
  
  
  landmarkWaveform = current_sort_info.landmarkWaveform;
  if current_sort_info.use_projection
    landmarkWaveformProj = ...
      current_sort_info.projectDirections'*landmarkWaveform;
  end
  % Add time as a coordinate
  sniptime = getappdata(handles.figCass,'absoluteSnipTime');
  n = length(sniptime);
  for nth = 1:n
    if size(sniptime{nth},1) == 0
      sniptime{nth} = (sniptime{nth})';
    end
  end
  sniptime = cat(2,sniptime{:});
  t2V = current_sort_info.t2V;
  if ~t2V
    t2V = 1;  % This will allow user to see time, even if not used for sorting
  end
  snipProj = [snipProj; sniptime*t2V];
  landmarkWaveformProj = [landmarkWaveformProj; ...
    current_sort_info.landmarkT*t2V];
  setappdata(handles.figCass,'snipProj',snipProj);
  setappdata(handles.figCass,'landmarkWaveformProj',landmarkWaveformProj);
  
  
function update_landmarkIndex(handles,rmsflag)
  isFake=getappdata(handles.figCass, 'isFake');
  if(isempty(isFake))
     isFake=false; % not strictly necessary, but make code clear
  end

  if(isFake)
     currentChannel = getappdata(handles.figCass,'currentChannel');
     shc = getappdata(handles.figCass, 'shc');
     snipIndicesUsed=getappdata(handles.figCass, 'snipIndicesUsed');
     nfiles=length(shc);
     snipLabels=cell(1, nfiles);
     for fileIndex=1:nfiles
        ssnpFile=fullfile(shc(fileIndex).fh.abspathstr, shc(fileIndex).fh.filename);
        fakeChannelInfoFile=replace_extension(ssnpFile, '.fake_info');
        tt=load(fakeChannelInfoFile, '-mat');
        fakeInfo=tt.fakeInfo;
        fakeChannelIndex=find([fakeInfo.templateID]==currentChannel);
        snipLabels{fileIndex}=fakeInfo(fakeChannelIndex).snipLabels;
     end % for, each file
     landmarkIndex=zeros(1, size(snipIndicesUsed,2) );
     % TODO: next loop is  inefficient
     for col=1:length(landmarkIndex)
        fileIndex=snipIndicesUsed(1,col);
        snipIndex=snipIndicesUsed(2,col);
        landmarkIndex(col)=snipLabels{fileIndex}(snipIndex);
     end
  else
     current_sort_info = getappdata(handles.figCass,'current_sort_info');
     snipProj = getappdata(handles.figCass,'snipProj');
     landmarkWaveformProj = getappdata(handles.figCass,'landmarkWaveformProj');
     % Determine the closest landmark to each snippet
     [dist,landmarkIndex] = mindist(snipProj,landmarkWaveformProj);
  end
  
  setappdata(handles.figCass,'landmarkIndex',landmarkIndex);
  
%   if (nargin > 1 && rmsflag)
%     clabel = agglabel(landmarkIndex);
%     nLandmarks = size(landmarkWaveformProj,2);
%     rmsd = zeros(1,nLandmarks);
%     for i = 1:length(clabel)
%       if ~isempty(clabel{i})
%         rmsd(i) = sqrt(mean(dist(clabel{i})));
%       end
%     end
%     setappdata(handles.figCass,'landmarkRMSD',rmsd);
%   end
  
function status = update_landmarkClust(handles)
% Re-cluster landmarks
% status = 1 means completed successfully, 0 is cancel
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  % Decide on alphaR
  %alphaRIndex = get(handles.popupAlphaR,'Value');
  alphaR = current_sort_info.Rfactor;
  if isnan(alphaR)
    alphaR = inputdlg('Choose alphaR:');
    if isempty(alphaR)
      status = 0; % abort
      return
    end
    alphaR = str2num(alphaR);
    current_sort_info.Rfactor = alphaR;
  end
  snipProj = getappdata(handles.figCass,'snipProj');
  landmarkWaveformProj = getappdata(handles.figCass, ...
                                    'landmarkWaveformProj');
  %landmarkRMSD = getappdata(handles.figCass, 'landmarkRMSD');
  % Get rid of R = 0 cases
  %nzIndex = landmarkRMSD > 0;
  %landmarkRMSD(~nzIndex) = mean(landmarkRMSD(nzIndex));
  % Do the clustering
  options = getappdata(handles.figCass,'options');
  cops = struct('ploteach',get(handles.boxWatchClustering,'Value'));
  if cops.ploteach
    hfigwatch = figure;
    cops.ploteach = gca;
  end
  landmarkClust = options.cluster_func(snipProj,...
                                       landmarkWaveformProj,...
                                       alphaR,...
                                       cops);
  current_sort_info.landmarkClust = landmarkClust;
  setappdata(handles.figCass,'current_sort_info',current_sort_info);
  setappdata(handles.figCass,'stabilityValid',0);
  status = 1;

  
function update_cluster_axes_calculations(handles)
  % First zero everything in case of empty or error
  setappdata(handles.figCass,'cluster_axes_data',[]);
  if (~isappdata(handles.figCass,'snipIndex') || ...
      ~isappdata(handles.figCass,'current_sort_info'))
    return
  end
  % Fetch the needed data
  snipIndex = getappdata(handles.figCass,'snipIndex');
  landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  peak_height_data = getappdata(handles.figCass,'peak_height_data');
  peakHeight = getappdata(handles.figCass,'peakHeight');
  peakHeight = cat(2,peakHeight{:});
  % Collect waveforms associated with each landmark into clusters
  clabel = agglabel(current_sort_info.landmarkClust+2);  % +2 for
                                                         % delete=-1
  clabel(1) = []; % get rid of the delete pile
  indexToShow = cat(2,clabel{:});  % Make sure our limits don't include delete
  ylim = [min(min(current_sort_info.landmarkWaveform(:,indexToShow))) ...
    max(max(current_sort_info.landmarkWaveform(:,indexToShow)))];
  for i = 1:length(clabel)
    cluster_axes_data(i).landmarkWaveforms = ...
        current_sort_info.landmarkWaveform(:,clabel{i});
    cluster_axes_data(i).waveformXLim = ...
        [1 size(current_sort_info.landmarkWaveform,1)];
    cluster_axes_data(i).waveformYLim = ylim;
  end
  % Distribute the peak height data into the clusters
  snipClust = current_sort_info.landmarkClust(landmarkIndex);
  clabel = agglabel(snipClust+2);
  clabel(1) = [];     % Get rid of the delete pile
  nClusters = length(clabel);
  for i = 1:nClusters
    tClust = clabel{i};
    if ~isempty(tClust)
      tmp = histc(peakHeight(tClust),peak_height_data.binEdges);
      tmp(end-1) = tmp(end-1)+tmp(end); % Add in number of pts on edge
      cluster_axes_data(i).peakHeightBinned = tmp(1:end-1);
    else
      cluster_axes_data(i).peakHeightBinned = ...
        zeros(1,length(peak_height_data.binEdges)-1);
    end
    cluster_axes_data(i).num_in_cluster = length(tClust);
  end
  setappdata(handles.figCass,'cluster_axes_data',cluster_axes_data);
  % Now we have to update axes_map
  axes_map = struct('virtual_axis',  num2cell(1:length(clabel)), ...
                  'selected', repmat({'off'},1,length(clabel)) ...
                  ); % 0 for noise, a positive number for each cluster
  setappdata(handles.figCass,'axes_map',axes_map);


  
function update_alphaR_range_calculations(handles)
  % Determine the "range" of alphaRs over which a cluster is "stable"
  % Here's the definition of instability:
  % Consider two sets of cluster labels, "candidate" (which might be
  % from alphaR(chosen), or might include the result of merges) and
  % "automatically-identified" at alphaR(other).  If there exists an
  % automatically-identified cluster at alphaR(other) such that it
  % contains landmarks that have different cluster numbers in the candidate
  % clustering, then: both candidate clusters are unstable at
  % alphaR(other).
  setappdata(handles.figCass,'alphaRIndexStable',[]);
  stabilityValid = getappdata(handles.figCass,'stabilityValid');
  if ~isempty(stabilityValid)
    if ~stabilityValid
      return
    end
  end
  sort_info = getappdata(handles.figCass,'sort_info');
  if isempty(sort_info)
    return;
  end
  replicate = 1;
  if (length(sort_info) > 1)
    replicate = get(handles.popupReplicate,'Value');
  end
  if ~isfield(sort_info,'Rfactor') || (length(sort_info(replicate).Rfactor) == 1)
    % We only have 1 alphaR, so this is meaningless
    return
  end
  % We have multiple alphaR, so we can proceed
  alphaRIndex = get(handles.popupAlphaR,'Value');
  landmarkClust = sort_info(replicate).landmarkClust;
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  clandmarkClust = current_sort_info.landmarkClust;
  n_alphaR = size(landmarkClust,1);
  n_candidate_clust = max(clandmarkClust);
  alphaRIndexStable = nan(n_alphaR,n_candidate_clust);
  for i = 1:n_alphaR
    tClust = landmarkClust(i,:);
    clabel = agglabel(tClust+2);  % +2 for delete=-1
    clabel(1) = [];  % Eliminate the delete pile
    tStable = ones(1,n_candidate_clust);
    for j = 1:length(clabel)  % for each automatically-identified cluster
      candidate_labels = clandmarkClust(clabel{j});
      u_c_labels = unique(candidate_labels);  % get candidate cluster lbls
      u_c_labels = u_c_labels(u_c_labels > -1); % get rid of delete pile
      if (length(u_c_labels) > 1)  % if not unique...
        u_c_labels = u_c_labels(find(u_c_labels)); % get rid of noise
        tStable(u_c_labels) = 0;   % all participating candidates are unstable
      end
    end
    alphaRIndexStable(i,:) = tStable;
  end
  % Store the results
  setappdata(handles.figCass,'alphaRIndexStable',alphaRIndexStable);
  
  
%
% Display updating functions
%
function update_display(handles)
  if ~isappdata(handles.figCass,'sort_info')
    return
  end
  update_cluster_axes_display(handles)
  update_time_axis_display(handles)
  update_projection_axis_display(handles)
  update_stability_axis_display(handles)
  % NOTE: if this section changes substantially, should go back and
  % integrate any necessary changes into btnMerge_Callback below, which is
  % designed to carry out everything of this function except for
  % update_time_axis_disply
  
  
function update_cluster_axes_display(handles)
  cluster_axes_data = getappdata(handles.figCass,'cluster_axes_data');
  % If there are no data, clear everything
  n_visible_axes = size(handles.handle_matrix,2);
  if (1)
  %if isempty(cluster_axes_data)
    for i = 1:2
      for j = 1:n_visible_axes
        cla(handles.handle_matrix(i,j));
      end
    end
    set(handles.handle_matrix(3,1),'String','Noise:');
    for j = 2:n_visible_axes
      set(handles.handle_matrix(3,j),'String',[num2str(j-1) ':']);
    end
  end
  if isempty(cluster_axes_data)
    return
  end
  axes_map = getappdata(handles.figCass,'axes_map');
  axes_map_start_from = getappdata(handles.figCass, ...
                                   'axes_map_start_from');
  peak_height_data = getappdata(handles.figCass,'peak_height_data');
  
  % Determine which clusters are showing
  virtual_axis = [axes_map.virtual_axis];
  clusterIndex = find(virtual_axis >= axes_map_start_from & ...
                      virtual_axis < axes_map_start_from + n_visible_axes);
  % clusterIndex = clusterNumber+1, because of noise=0

  % Display each cluster
  for i = 1:length(clusterIndex)
    tIndex = clusterIndex(i);
    col = unique_color(tIndex,length(axes_map));
    % Show waveforms
    waveformsToShow = cluster_axes_data(tIndex).landmarkWaveforms;
    if ~isempty(waveformsToShow)
      plot(handles.handle_matrix(1,i),...
        waveformsToShow,...
        'Color',col);
    else
      cla(handles.handle_matrix(1,i));
    end
    if (i > 1)
      set(handles.handle_matrix(1,i),'YTick',[])
    end
    % Histogram of peak heights
    % Use stairs + patch instead of bar for performance reasons
    binEdges = peak_height_data.binEdges;
    nPerBin = cluster_axes_data(tIndex).peakHeightBinned;
    [x,y] = stairs(binEdges([1 1:end]),...
                   [0 nPerBin 0]);
    cla(handles.handle_matrix(2,i));
    hpatch = patch(x,y,col,'Parent',handles.handle_matrix(2,i),...
      'EdgeColor',col);
    % Show threshold line(s)
    yrange = [0 max(nPerBin)*1.05+1]';
    line(peak_height_data.thresh([1 1],:),...
         yrange(:,ones(1,length(peak_height_data.thresh))),...
         'Parent',handles.handle_matrix(2,i),...
         'LineStyle','--',...
         'Color','k');
    xlim = [min(min(peak_height_data.thresh),min(binEdges)) ...
      max(max(peak_height_data.thresh),max(binEdges))];
    xlim = xlim + [-.05 .05]*diff(xlim);
    % Display voltage min/max lines
    if (peak_height_data.voltageRange(1) > xlim(1))
      line(peak_height_data.voltageRange([1 1]),...
           yrange,...
           'Parent',handles.handle_matrix(2,i),...
           'LineStyle','-',...
           'Color','k');
    end
    if (peak_height_data.voltageRange(2) < xlim(2))
      line(peak_height_data.voltageRange([2 2]),...
           yrange,...
           'Parent',handles.handle_matrix(2,i),...
           'LineStyle','-',...
           'Color','k');
    end
    set(handles.handle_matrix(2,i),'XLim',xlim,'YLim',yrange);
    % Update the text strings
    if (tIndex == 1)
      clustString = ['Noise: ' ...
                     num2str(cluster_axes_data(tIndex).num_in_cluster)];
    else
      clustString = [num2str(tIndex-1) ': ' ...
                     num2str(cluster_axes_data(tIndex).num_in_cluster)];
    end
    set(handles.handle_matrix(3,i),'String',clustString);
  end
  % Set waveform axes to correct scale (in fact, common across all)
  set(handles.handle_matrix(1,:),...
      'XLim',cluster_axes_data(1).waveformXLim,...
      'YLim',cluster_axes_data(1).waveformYLim);
  % Prettify the plot
  set(handles.handle_matrix([1 2],:),'XTick',[],'YTick',[],'Box','off');
  set(handles.handle_matrix(1,1),'YTickMode','auto');
  % Set the selection handles
  selected = {axes_map(clusterIndex).selected};
  selected = [selected repmat({'off'},1,n_visible_axes-length(selected))];
  set(handles.handle_matrix(1,:),{'Selected'},selected');
  

function update_projection_axis_display(handles)
  % Before we do anything to modify the axis, we need to check and see if
  % we need to preserve axis limits (for zooming)
  xlim = [];
  if isappdata(handles.axesVisualizeShape,'OXLim')
    % If x data are valid from a previous draw, store the _current_ xlim
    xlim = get(handles.axesVisualizeShape,'XLim');
  else
    % If not valid, make sure the axis limits are set to default
    set(handles.axesVisualizeShape,'XLimMode','auto')
  end
  ylim = [];
  if isappdata(handles.axesVisualizeShape,'OYLim')
    ylim = get(handles.axesVisualizeShape,'YLim');
  else
    set(handles.axesVisualizeShape,'YLimMode','auto')
  end
  % Now draw everything from scratch
  cla(handles.axesVisualizeShape);
  snipProj = getappdata(handles.figCass,'snipProj');
  keepflag = true(1,size(snipProj,2));
  if isappdata(handles.figCass,'snipRank')
    snipRank = getappdata(handles.figCass,'snipRank');
    rankThreshold = get(handles.popupEventRank,'Value');
    keepflag = (cat(2,snipRank{:}) <= rankThreshold);
  end
  snipProj = snipProj(:,keepflag);
  landmarkWaveformProj = getappdata(handles.figCass,'landmarkWaveformProj');
  landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
  landmarkIndex = landmarkIndex(keepflag);
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  axes_map=getappdata(handles.figCass, 'axes_map');
  x_coord = get(handles.popupXVis,'Value');
  y_coord = get(handles.popupYVis,'Value');
  landmarkClust = current_sort_info.landmarkClust;
  clabelLandmark = agglabel(landmarkClust+2); % +2 for delete = -1
  clabelLandmark(1) = []; % get rid of delete group
  clabelSnip = agglabel(landmarkClust(landmarkIndex)+2);
  clabelSnip(1) = [];
  %set(handles.axesVisualizeShape,'XLimMode','auto','YLimMode','auto')
  landmarkVisibleText = {'off','on'};
  landmarkVisible = landmarkVisibleText{get(handles.boxShowLandmarks,...
                                            'Value')+1};      
  for i = 1:length(clabelSnip)
    if ~isempty(clabelSnip{i})
      col = unique_color(i,length(clabelSnip));
      hline = line(snipProj(x_coord,clabelSnip{i}),...
           snipProj(y_coord,clabelSnip{i}),...
           'Parent',handles.axesVisualizeShape,...
           'LineStyle','none',...
           'Marker','.',...
           'MarkerSize',1,...
           'Color',col);
      axes_map(i).projLineHandle = hline;
      if (strcmp(axes_map(i).selected,'on'))
        set(hline,'MarkerSize',9);
      end
      % Display the landmarks
      hline = line(landmarkWaveformProj(x_coord,clabelLandmark{i}),...
           landmarkWaveformProj(y_coord,clabelLandmark{i}),...
           'Parent',handles.axesVisualizeShape,...
           'LineStyle','none',...
           'Marker','o',...
           'MarkerSize',12,...
           'MarkerFaceColor',col,...
           'MarkerEdgeColor',col,...
           'Visible',landmarkVisible,...
           'Tag','landmarkProjection');
         axes_map(i).landmarkLineHandle = hline;
    end
  end
  set(handles.axesVisualizeShape,'XTick',[],'YTick',[])
  % Put a '+' at 0,0 to facilitate identification of noise cloud
  line(0,0,'Parent',handles.axesVisualizeShape,'Marker','+','Color','k');
  %axis(handles.axesVisualizeShape,'equal')
  if get(handles.boxLockAspectRatio,'Value')
    set(handles.axesVisualizeShape,'DataAspectRatio',[1 1 1]);
  else
    set(handles.axesVisualizeShape,'DataAspectRatioMode','auto');
  end
  % Retrieve/store axis limits (for zooming)
  if ~isempty(xlim)
    set(handles.axesVisualizeShape,'XLim',xlim);
  else
    setappdata(handles.axesVisualizeShape,'OXLim',...
      get(handles.axesVisualizeShape,'XLim'));
  end
  if ~isempty(ylim)
    set(handles.axesVisualizeShape,'YLim',ylim);
  else
    setappdata(handles.axesVisualizeShape,'OYLim',...
      get(handles.axesVisualizeShape,'YLim'));
  end
  setappdata(handles.figCass,'axes_map',axes_map);
  

function update_time_axis_display(handles)
  t = getappdata(handles.figCass,'absoluteSnipTime');
  if isempty(t)
    return
  end
  n = length(t);
  for nth = 1:n
    if size(t{nth},1) == 0
      t{nth} = (t{nth})';
    end
  end
  t = cat(2,t{:});
  landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  landmarkClust = current_sort_info.landmarkClust;
  clabel = agglabel(landmarkClust(landmarkIndex)+2);  % +1 for delete=-1
  clabel(1) = [];  % Get rid of the delete bin
  cla(handles.axesTime);
  % Draw the spike ticks
  for i = 1:length(clabel)
    if ~isempty(clabel{i})
      line(t(clabel{i}),...
           repmat(i,1,length(clabel{i})),...
           'Parent',handles.axesTime,...
           'LineStyle','none',...
           'Marker','.',...
           'MarkerSize',1,...
           'HitTest','off',...
           'Color',unique_color(i,length(clabel)));
    end
  end
  % Draw the landmark ticks
  landmarkVisibleText = {'off','on'};
  landmarkVisible = landmarkVisibleText{get(handles.boxShowLandmarks,...
                                            'Value')+1};      
  tickSize = 0.2;
  for i = 1:length(landmarkClust)
    if (landmarkClust(i) > -1)
      col_index = landmarkClust(i)+1;
      tick_yrange = landmarkClust(i)+1+[-1 1]*tickSize;
      line(current_sort_info.landmarkT([i i]),...
           tick_yrange,...
           'Parent',handles.axesTime,...
           'Color',unique_color(col_index,length(clabel)),...
           'Visible',landmarkVisible,...
           'HitTest','off',...
           'Tag','landmarkT');
    end
  end
  ylim = [0 length(clabel)+1];
  % Draw the time markers
  cass_draw_time_markers(handles,[0 max(t)],ylim,handles.axesTime);
  % Finish it off
  set(handles.axesTime,'YLim',ylim,...
    'XLim',[0 max(t)],...
    'TickDir','out',...
    'YTick',[]);
  xlabel(handles.axesTime,'Time (s)')
  

function update_stability_axis_display(handles)
  cla(handles.axesAlphaRRange)
  if ~isappdata(handles.figCass,'alphaRIndexStable')
    return;
  end
  alphaRIndexStable = getappdata(handles.figCass,'alphaRIndexStable');
  if isempty(alphaRIndexStable)
    return
  end
  n_candidate_clusters = size(alphaRIndexStable,2);
  sort_info = getappdata(handles.figCass,'sort_info');
  replicate = 1;
  if (length(sort_info) > 1)
    replicate = get(handles.popupReplicate,'Value');
  end
  alphaR = sort_info(replicate).Rfactor;
  alphaRIndex = get(handles.popupAlphaR,'Value');  %index of current alphaR
  set(handles.axesAlphaRRange,'XTickMode','auto')
  for i = 1:size(alphaRIndexStable,2)
    tStable = find(alphaRIndexStable(:,i)); % These are the stable alphaRs
    if ~isempty(tStable)
      tAlphaR = nan(size(alphaR)); % use NaNs to create line breaks
      tAlphaR(tStable) = alphaR(tStable);
      col = unique_color(i+1,n_candidate_clusters+1);  %+1 for noise
      line(repmat(i,1,length(tAlphaR)),...
           tAlphaR,...
           'Parent',handles.axesAlphaRRange,...
           'LineStyle','-',...
           'Marker','+',...
           'Color',col,...
           'MarkerFaceColor',col,...
           'MarkerEdgeColor',col);
    end
  end
  % Prettify the plot
  % Set up the YTicks: use up to 4 of them, including beginning, current
  % alphaR, and end
  ytickIndex = [1 alphaRIndex length(alphaR)];
  ytick = alphaR(ytickIndex);
  % Place the 4th in the largest gap (logarithmically)
  log_ytick = log(ytick);
  [tmp,gapIndex] = max(diff(log_ytick));
  ytick_add = mean(log_ytick([gapIndex gapIndex+1]));
  [tmp,ytick_add_index] = min(abs(log(alphaR) - ytick_add));
  ytickIndex = unique([ytickIndex ytick_add_index]);
  set(handles.axesAlphaRRange,...
      'XLim',[0 n_candidate_clusters+1],...
      'YLim',alphaR([1 end]),...
      'YScale','log',...
      'TickLength',[0.03 0.05],...
      'YTick',alphaR(ytickIndex),...
      'YTickLabel',num2str(alphaR(ytickIndex)'),...
      'YMinorTick','off');
  xtick = get(handles.axesAlphaRRange,'XTick');
  xtick = setdiff(xtick,[0 n_candidate_clusters+1]);
  xtick = xtick(find(xtick == round(xtick)));
  set(handles.axesAlphaRRange,'XTick',xtick);
  ylabel(handles.axesAlphaRRange,'a_R','Interpreter','none');
  
%
% Utility functions
%
function update_selections(handles,clIndex,clusterAxHandle)
  axes_map=getappdata(handles.figCass, 'axes_map');
  is_shift_down=strcmp(get(handles.figCass, 'SelectionType'), ...
                       'extend');
  current_state = axes_map(clIndex).selected;
  if(is_shift_down) % toggle current axes only
    if(strcmp(current_state, 'on'))
      axes_map(clIndex).selected = 'off';
      set(axes_map(clIndex).projLineHandle,'MarkerSize',1);
      if ishandle(clusterAxHandle)
        set(clusterAxHandle, 'selected', 'off');
      end
    else
      axes_map(clIndex).selected = 'on';
      set(axes_map(clIndex).projLineHandle,'MarkerSize',9);
      if ishandle(clusterAxHandle)
        set(clusterAxHandle, 'selected', 'on');
      end
    end
  else % set current axes object selected, all other unselected
    set(handles.handle_matrix(1,:), 'selected', 'off');
    [axes_map.selected] = deal('off');
    if ishandle(clusterAxHandle)
      set(clusterAxHandle, 'selected', 'on');
    end
    axes_map(clIndex).selected = 'on';
    % Update marker sizes, too
    plH = [axes_map.projLineHandle];
    markersize = get(plH,'MarkerSize');
    if iscell(markersize)
      markersize = cell2mat(markersize);
    end
    indexFat = find(markersize > 1);
    set(plH(indexFat),'MarkerSize',1); %just reset ones that need it (faster)
    set(axes_map(clIndex).projLineHandle,'MarkerSize',9);
  end
  setappdata(handles.figCass, 'axes_map', axes_map);
    
    
% $$$ function update_selections(handles)
% $$$   return  % Now this is handled in axes_mouse_up
% $$$   % Determine which clusters are showing
% $$$   axes_map=getappdata(handles.figCass, 'axes_map');
% $$$   axes_map_start_from=getappdata(handles.figCass, 'axes_map_start_from');
% $$$   virtual_axis = [axes_map.virtual_axis];
% $$$   n_visible_axes = size(handles.handle_matrix,2);
% $$$   clusterIndex = find(virtual_axis >= axes_map_start_from & ...
% $$$                       virtual_axis < axes_map_start_from + n_visible_axes);
% $$$   % Get their selection state & copy it over to virtual axes
% $$$   selected = get(handles.handle_matrix(1,:),'Selected');
% $$$   for i = 1:length(clusterIndex)
% $$$     axes_map(clusterIndex(i)).selected = selected{i};
% $$$   end
% $$$   setappdata(handles.figCass, 'axes_map',axes_map);
  

function invert_selections(handles)
  axes_map=getappdata(handles.figCass, 'axes_map');
  offIndex = strmatch('off',{axes_map.selected});
  onIndex = setdiff(1:length(axes_map),offIndex);
  [axes_map(offIndex).selected] = deal('on');
  [axes_map(onIndex).selected] = deal('off');
  onHandles = [axes_map(onIndex).projLineHandle];
  offHandles = [axes_map(offIndex).projLineHandle];
  set(onHandles,'MarkerSize',1);
  set(offHandles,'MarkerSize',9);
  setappdata(handles.figCass, 'axes_map', axes_map);
  % Figure out which axes are showing, and toggle their selection state
  % This way may be more laborious than needed, but it does insure
  % that "empty" axes don't get their selection state turned on
  axes_map_start_from=getappdata(handles.figCass, 'axes_map_start_from');
  virtual_axis = [axes_map.virtual_axis];
  n_visible_axes = size(handles.handle_matrix,2);
  clusterIndex = find(virtual_axis >= axes_map_start_from & ...
                      virtual_axis < axes_map_start_from + ...
                      n_visible_axes);
  [tmpcI,axOffIndex] = intersect(clusterIndex,offIndex);
  set(handles.handle_matrix(1,:),'Selected','off');
  set(handles.handle_matrix(1,axOffIndex),'Selected','on');


function cass_draw_time_markers(handles,tlim,ylim,hax)
  is_axesTime = (hax == handles.axesTime); % axesTime or zoom window?
  % Get the file boundaries within time interval
  tstart = getappdata(handles.figCass,'absoluteFileStartTime');
  tstart = tstart-min(tstart);
  tlength = getappdata(handles.figCass,'fileDuration');
  tend = tstart+tlength;
  selectIndex = find(tstart >= tlim(1) & tstart <= tlim(2));
  tstart = tstart(selectIndex);
  selectIndex = find(tend >= tlim(1) & tend <= tlim(2));
  tend = tend(selectIndex);
  % Figure out how many clusters we have
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  landmarkClust = current_sort_info.landmarkClust+1;  % +1 because noise=0
  maxclust = max(landmarkClust);
  % Draw marker lines at file start & end times
  for i = 1:length(tstart)
    line(tstart([i i])+tlength(i),...
         ylim,...
         'Parent',hax,...
         'LineStyle',':',...
         'Marker','none',...
         'HitTest','off',...
         'Color','k');
    line(tstart([i i]),...
         ylim,...
         'Parent',hax,...
         'LineStyle','--',...
         'Marker','none',...
         'HitTest','off',...
         'Color','k');
  end
  % Draw manually-set time markers
  timeMarker = getappdata(handles.figCass,'timeMarker');
  if (length(timeMarker) > 0)
    % Find the time markers that fall within the time limits
    indexKeep = find([timeMarker.relativeTime] >= tlim(1) & ...
      [timeMarker.relativeTime] <= tlim(2));
    timeMarker = timeMarker(indexKeep);
  end
  markerString = {'mark','start','stop'};
  markerStyle = {'.','o','x'};
  if (length(timeMarker) > 0)
    % Set up data for dragging time markers
    drag_marker_info.done = @marker_dragged;
    if is_axesTime
      % Set up context menu for deleting & editing time markers
      hcontext = uicontextmenu;
      uimenu(hcontext,'Label','Delete',...
        'Callback',@delete_time_marker);
      uimenu(hcontext,'Label','Edit comment',...
        'Callback',@edit_time_marker_comment);
    else
      hcontext = [];
    end
  end
  tickSize = 0.2;
  for i = 1:length(timeMarker)
    if isnan(timeMarker(i).clustNum) | (timeMarker(i).clustNum == -1) % means fake clust num from multielectode sorting
      tmy = ylim;
      col = 'k';
    else
      if is_axesTime
        tmy = timeMarker(i).clustNum+1 + [-1 1]*tickSize;
      else
        tmy = ylim;
      end
      %col = unique_color(timeMarker(i).clustNum+1,ylim(2)-1);
      col = unique_color(timeMarker(i).clustNum+1,maxclust);
    end
    tmstyle = markerStyle{strmatch(timeMarker(i).action(1:4), ...
                                   markerString)};
    drag_marker_info.doneargs = {struct('dataHandle',handles.figCass,...
                                       'startpos',timeMarker(i).relativeTime)};
    hmarkline = line(timeMarker(i).relativeTime([1 1]),tmy,...
                     'Color',col,...
                     'Marker',tmstyle,...
                     'EraseMode','xor',...
                     'ButtonDownFcn',{@moveline,drag_marker_info},...
                     'UIContextMenu',hcontext,...
                     'Parent',hax);
  end
  

function cass_waveform_reconstruct(handles,tlim)
  % Reconstruct the waveform from snippets
  sorthead = getappdata(handles.figCass,'sorthead');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  shc = sortheader_importchan(sorthead,currentChannel);
  snipminmax = getappdata(handles.figCass,'snipminmax');
  snipIndex = getappdata(handles.figCass,'snipIndex');
  % Learn the cluster assignment of each of these spikes
  landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  landmarkClust = current_sort_info.landmarkClust;
  clust = landmarkClust(landmarkIndex);
  % Split clust back out to cell array
  n_files = length(snipIndex);
  clust_cell = cell(1,n_files);
  offset = 0;
  for fileIndex = 1:n_files
    n_snips = length(snipIndex{fileIndex});
    clust_cell{fileIndex} = clust((1:n_snips)+offset);
    offset = offset+n_snips;
  end
  % Create a new figure for reconstruction plot, occupying the upper 20%
  % of the screen
  units_root = get(0,'Units');
  set(0,'Units','pixels');
  screensize = get(0,'ScreenSize');
  pos = screensize;
  pos(2) = 0.8*(screensize(2)+screensize(4));
  pos(4) = (screensize(2)+screensize(4)-pos(2));
  hfig = figure('Units','pixels','Position',round(pos));
  pos(2) = screensize(2); % sliderwindow at bottom
  hax_reconstruct = gca;
  % Call the standalone function
  display_reconstruction(hax_reconstruct,shc,snipminmax,clust_cell,snipIndex);
  ylim = get(hax_reconstruct,'YLim');
  % Now zoom in
  display_reconstruction(gca,'XLim',tlim);
  % Draw the time markers (file boundaries & manually-set markers)
  cass_draw_time_markers(handles,tlim,ylim,hax_reconstruct);
  % Create sliderwindow
  sliderwindow(hax_reconstruct,struct('axisupdatefcn',@display_reconstruction,'position',round(pos)))
  return
  
  % Old code
  % Find spikes that occured within time interval
  t = getappdata(handles.figCass,'absoluteSnipTime');
  t = cat(2,t{:});
  selectIndex = find(t >= tlim(1) & t <= tlim(2));
  t = t(selectIndex);
  snip = getappdata(handles.figCass,'snip');
  snip = cat(2,snip{:});
  snip = snip(:,selectIndex);
  % Learn the cluster assignment of each of these spikes
  landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  landmarkClust = current_sort_info.landmarkClust+1;  % +1 because noise=0
  clust = landmarkClust(landmarkIndex);
  clust = clust(selectIndex);
  % Get the colors
  maxclust = max(landmarkClust);
  for i = 1:maxclust
    colmatrix(i,:) = unique_color(i,maxclust);
  end
  % Determine the timing of each snippet
  sorthead = getappdata(handles.figCass,'sorthead');
  spiketiming = (sorthead(1).sniprange(1):sorthead(1).sniprange(2))...
      /sorthead(1).scanrate;
  % Create a new figure for reconstruction plot, occupying the upper 20%
  % of the screen
  units_root = get(0,'Units');
  set(0,'Units','pixels');
  screensize = get(0,'ScreenSize');
  pos = screensize;
  pos(2) = 0.8*(screensize(2)+screensize(4));
  pos(4) = (screensize(2)+screensize(4)-pos(2));
  hfig = figure('Units','pixels','Position',round(pos));
  pos(2) = screensize(2); % sliderwindow at bottom
  hax_reconstruct = gca;
  % Draw the spikes; first draw the noise spikes so that in case of
  % waveform overlap these get over-written by the "real" spikes
  noiseIndex = find(clust == 1);  % we've already added 1 to cluster number
  for i = noiseIndex
    line(spiketiming+t(i),snip(:,i),'Color',colmatrix(1,:));
  end
  restIndex = find(clust > 1);
  for i = restIndex
    line(spiketiming+t(i),snip(:,i),'Color',colmatrix(clust(i),:));
  end
  ylim = [min(snip(:)) max(snip(:))];
  ylim = ylim + [-0.1 0.1]*diff(ylim);
  % Draw the time markers (file boundaries & manually-set markers)
  cass_draw_time_markers(handles,tlim,ylim,hax_reconstruct);
  set(gca,'YLim',ylim)
  % Create a slider window (for zooming)
  sliderwindow(hax_reconstruct,struct('position',round(pos)))
  
  
    
  
function merge_selected(handles,label)
  axes_map = getappdata(handles.figCass,'axes_map');
  selected = {axes_map.selected};
  selectedIndex = strmatch('on',selected,'exact');
  if isempty(selectedIndex)
    return
  end
  if (nargin < 2)
    label = selectedIndex(1)-1;
  end
  if (length(selectedIndex) > 1 || label == -1)
    % Get the current data
    current_sort_info = getappdata(handles.figCass,'current_sort_info');
    timeMarker = getappdata(handles.figCass,'timeMarker');
    % Save the data for "undo" functionality
    setappdata(handles.figCass,'old_sort_info',current_sort_info);
    setappdata(handles.figCass,'old_timeMarker',timeMarker);
    set(handles.UndoBtn,'Enable','on')
    % Build a lookup table for resetting the cluster numbers of the
    % landmarks
    landmarkClust = current_sort_info.landmarkClust;
    newClustLabel = 0:max(landmarkClust);
    newClustLabel(selectedIndex) = label;
    newClustLabel = [-1 newClustLabel]; % for delete = -1
    [ulabel,tmp,ulabelmap] = unique(newClustLabel);
    newClustLabel = ulabelmap-1 + min(ulabel);
    landmarkClust = newClustLabel(landmarkClust+2); %+2 for delete=-1
    current_sort_info.landmarkClust = landmarkClust;
    % If there are time markers, we may have to merge those too
    for i = 1:length(timeMarker)
      if ~isnan(timeMarker(i).clustNum)
        timeMarker(i).clustNum = newClustLabel(timeMarker(i).clustNum+2);
      end
    end
    if (~any(selectedIndex == 1) && label ~= -1)
      % If we're doing anything other than putting spikes into the noise
      % bin, mark this sorting as manually tweaked
      current_sort_info.Rfactor = NaN;
    end
    setappdata(handles.figCass,'current_sort_info',current_sort_info)
    setappdata(handles.figCass,'timeMarker',timeMarker)
  end
  set(handles.handle_matrix,'Selected','off')
  update_cluster_axes_calculations(handles)
  setappdata(handles.figCass,'modified',1)
  

function cass_split_selected_cluster(handles,startPoint,currentPoint)
  % Find the equation of the line between the two clicks
  dxy = currentPoint - startPoint;
  n = [dxy(2), -dxy(1)];
  c = n*currentPoint(:);   % Eq. is n \cdot x = c
  % Get the currently selected cluster
  axes_map = getappdata(handles.figCass,'axes_map');
  clIndex = strmatch('on',{axes_map.selected},'exact');
  if (length(clIndex) ~= 1)
    error('Must have 1 selected cluster');
  end
  % Get the xy coordinates of the landmarks in this cluster
  x_coord = get(handles.popupXVis,'Value');
  y_coord = get(handles.popupYVis,'Value');
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  landmarkIndex = find(current_sort_info.landmarkClust == clIndex-1);
  landmarkWaveformProj = getappdata(handles.figCass,'landmarkWaveformProj');
  xy = landmarkWaveformProj([x_coord y_coord],landmarkIndex);
  % Partition into two groups on opposite sides of line
  proj = n * xy;
  current_sort_info.landmarkClust(landmarkIndex(proj < c)) = length(axes_map);
  setappdata(handles.figCass,'current_sort_info',current_sort_info);
  handles.isSplitting = false;
  guidata(handles.figCass,handles);
  update_cluster_axes_calculations(handles)
  update_display(handles)
  %setappdata(handles.figCass,'modified',1)
  %update_calculations_rest(handles)

function set_time_marker(hObject, eventdata, handles, action)
  % Set a time marker in the time axes
  % Get click position
  startPoint = getappdata(handles.axesTime,'startPoint');
  clustNum = round(startPoint(2))-1;
  if strcmp(action(end-2:end),'all') ...
          | (clustNum == -1) % means it's a fake cluster, from multielec. data
    clustNum = NaN;
  end
  timeMarker = getappdata(handles.figCass,'timeMarker');
  timeMarker(end+1).relativeTime = startPoint(1);
  timeMarker(end).clustNum = clustNum;
  timeMarker(end).action = action;
  if strcmp(action(1:4),'mark')
    % Get the user to type in a comment
    timeMarker(end).comment = inputdlg('Describe this marker',...
                                       'Time marker',3);
  end
  setappdata(handles.figCass,'timeMarker',timeMarker);
  setappdata(handles.figCass,'modified',1);
  update_time_axis_display(handles)
  
  
% --- Outputs from this function are returned to the command line.
function varargout = cass_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editComment_Callback(hObject, eventdata, handles)
% hObject    handle to editComment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editComment as text
%        str2double(get(hObject,'String')) returns contents of editComment as a double
  setappdata(handles.figCass,'modified',1);
  

% --- Executes during object creation, after setting all properties.
function editComment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editComment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnScrollLeft.
function btnScrollLeft_Callback(hObject, eventdata, handles, scroll_step)
% hObject    handle to btnScrollLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if(nargin < 4) scroll_step=1; end
   % handles=guidata(hObject);
   %axes_map=getappdata(handles.figCass, 'axes_map');
   % make sure not over-scroll:
   %if(axes_map{2}-scroll_step<1) scroll_step=axes_map{2}-1; end
   %axes_map=[axes_map(1), num2cell(cat(2, axes_map{2:end})-scroll_step)];
   %% todo: this has bug
   %setappdata(handles.figCass, 'axes_map', axes_map);
   %update_selections(handles)
   axes_map_start_from=getappdata(handles.figCass, 'axes_map_start_from');
   axes_map_start_from=max(1,axes_map_start_from-scroll_step);
   setappdata(handles.figCass, 'axes_map_start_from',axes_map_start_from);
   update_cluster_axes_display(handles)



% --- Executes on button press in btnScrollRight.
function btnScrollRight_Callback(hObject, eventdata, handles, scroll_step)
% hObject    handle to btnScrollRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   if(nargin < 4) scroll_step=1; end
   % handles=guidata(hObject);
   %update_selections(handles)
   axes_map=getappdata(handles.figCass, 'axes_map');
   axes_map_start_from=getappdata(handles.figCass, 'axes_map_start_from');
   %axes_map=[axes_map(1), num2cell(cat(2, axes_map{2:end})+scroll_step)];
   %% todo: this has bug
   axes_map_start_from = min(axes_map_start_from+scroll_step,...
     max(1,length(axes_map)-size(handles.handle_matrix,2)+1));
   %setappdata(handles.figCass, 'axes_map', axes_map);
   setappdata(handles.figCass, 'axes_map_start_from', axes_map_start_from);
   
   update_cluster_axes_display(handles);


% --- Executes on button press in btJumpRight.
function btJumpRight_Callback(hObject, eventdata, handles)
% hObject    handle to btJumpRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  btnScrollRight_Callback(hObject,eventdata,handles,size(handles.handle_matrix,2));

% --- Executes on button press in btnJumpLeft.
function btnJumpLeft_Callback(hObject, eventdata, handles)
% hObject    handle to btnJumpLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  btnScrollLeft_Callback(hObject,eventdata,handles,size(handles.handle_matrix,2));


% --- Executes on button press in btnGoFirst.
function btnGoFirst_Callback(hObject, eventdata, handles)
% hObject    handle to btnGoFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   setappdata(handles.figCass, 'axes_map_start_from',1);
   update_cluster_axes_display(handles)

% --- Executes on button press in btnGoLast.
function btnGoLast_Callback(hObject, eventdata, handles)
% hObject    handle to btnGoLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  %axes_map = getappdata(handles.figCass,'axes_map');
  %axes_map_start_from = length(axes_map) - size(handles.handle_matrix,
  % setappdata(handles.figCass, 'axes_map_start_from',1);
  % update_cluster_axes_display(handles)
  btnScrollRight_Callback(hObject, eventdata, handles, inf);



% --------------------------------------------------------------------
function menuitemOpenProject_Callback(hObject, eventdata, handles)
% hObject    handle to menuitemOpenProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  options = getappdata(handles.figCass,'options');
  pathok = 0;
  while ~pathok
    pathin = uigetdir(pwd,'Choose project directory');
    if ~pathin
      return   % user hit cancel, quit trying
    end
    project_overview_name = [pathin filesep 'overview.mat'];
    if (exist(project_overview_name) ~= 2)
      errordlg(['Project overview file does not exist, please pick '...
                'a different directory']);
    else
      pathok = 1;
    end
  end
  options.dirname = pathin;
  setappdata(handles.figCass,'options',options);
  cass_start_project(handles)
  
  
function cass_start_project(handles)
  % Load the project overview file
  cass_parse_overview(handles)
  % Set the current channel
  options = getappdata(handles.figCass,'options');
  if isfield(options,'channels')
    startIndex = [];
    if isfield(options,'channel')
      startIndex = find(options.channel == options.channels);
    end
    if isempty(startIndex)
      startIndex = 1;
    end
    setappdata(handles.figCass,'currentChannel',options.channels(startIndex))
  else
    return
  end
  % Clear all the old channel-specific data
  clear_channel_data(handles)
  % Load the first channel's data
  % must load sort info before snippets, so we know how big projections are
  filename = cass_choose_sort_info_file(handles);
  cass_load_sort_info(handles,filename)
  cass_load_snippets(handles)
  % Get going!
  update_calculations(handles);
  update_display(handles);


% --------------------------------------------------------------------
function menuitemExit_Callback(hObject, eventdata, handles)
% hObject    handle to menuitemExit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  status = cass_save_check(handles);
  if status
    close(handles.figCass)
  end

% --------------------------------------------------------------------
function menuItemAbout_Callback(hObject, eventdata, handles)
% hObject    handle to menuItemAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuitemFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuitemFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuitemHelp_Callback(hObject, eventdata, handles)
% hObject    handle to menuitemHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupWaveform.
function popupWaveform_Callback(hObject, eventdata, handles)
% hObject    handle to popupWaveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupWaveform contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupWaveform


% --- Executes on button press in boxShowLandmarks.
function boxShowLandmarks_Callback(hObject, eventdata, handles)
% hObject    handle to boxShowLandmarks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boxShowLandmarks
  landmarkVisibleText = {'off','on'};
  landmarkVisible = landmarkVisibleText{get(hObject,'Value')+1};
  hlandmarks = findobj(handles.axesVisualizeShape,'Tag', ...
                       'landmarkProjection');
  set(hlandmarks,'Visible',landmarkVisible);
  hlandmarks = findobj(handles.axesTime,'Tag','landmarkT');
  set(hlandmarks,'Visible',landmarkVisible);


% --- Executes on selection change in popupAlphaR.
function popupAlphaR_Callback(hObject, eventdata, handles)
% hObject    handle to popupAlphaR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupAlphaR contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupAlphaR
  update_calculations(handles)
  update_display(handles)
  

% --- Executes on selection change in popupReplicate.
function popupReplicate_Callback(hObject, eventdata, handles)
% hObject    handle to popupReplicate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupReplicate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupReplicate
  if isappdata(handles.figCass,'landmarkIndex_total')
    rmappdata(handles.figCass,'landmarkIndex_total')
  end
  update_calculations(handles)
  update_display(handles)

  
% --- Executes on button press in ResetButton.
function ResetButton_Callback(hObject, eventdata, handles)
% hObject    handle to ResetButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  update_calculations(handles)
  update_display(handles)


% --- Executes on button press in btnCorrelation.
function btnCorrelation_Callback(hObject, eventdata, handles)
% hObject    handle to btnCorrelation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  sorthead = getappdata(handles.figCass,'sorthead');
  options = getappdata(handles.figCass,'options');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  shc = sortheader_importchan(sorthead,currentChannel);
  % The main issue is to go from the subset of spikes to the totality of
  % spikes.
  % Have we had to calculate assignments previously?
  if isappdata(handles.figCass,'landmarkIndex_total')
    landmarkIndex_total = getappdata(handles.figCass,...
                                     'landmarkIndex_total');
  else
    % Nope, this is our first time
    % Determine whether our "subset" of spikes is in fact the whole thing
    snipIndexStride = getappdata(handles.figCass,'snipIndexStride');
    have_all = (snipIndexStride == 1);
    if ~have_all
      % We have to load from disk & apply sorting
      current_sort_info = getappdata(handles.figCass,'current_sort_info');
      landmarkIndex_total = landmark_sort_apply(current_sort_info,shc);
      setappdata(handles.figCass,'landmarkIndex_total',landmarkIndex_total);
    else
      % or, our "subset" is indeed the whole thing
      landmarkIndex = getappdata(handles.figCass,'landmarkIndex');
      % but we have to split it back to individual files
      nsnips_per_file = getappdata(handles.figCass,'nsnips_per_file');
      cum_snips_per_file = [0 cumsum(nsnips_per_file)];
      for i = 1:length(nsnips_per_file)
        landmarkIndex_total{i} = ...
          landmarkIndex(cum_snips_per_file(i)+1:cum_snips_per_file(i+1));
      end
    end
  end
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  landmarkClust = current_sort_info.landmarkClust;
  % Determine the cluster number assigned to each spike
  for i = 1:length(landmarkIndex_total)
    spikeClust{i} = landmarkClust(landmarkIndex_total{i});
  end
  % Set options
  if ~isfield(options,'corr_tmin')
    options.corr_tmin = 0;  % will choose a single scan
  end
  % One option: study only the selected clusters
  axes_map = getappdata(handles.figCass,'axes_map');
  selected = {axes_map.selected};
  selectedIndex = strmatch('on',selected,'exact');
  if isempty(selectedIndex)
    selectedIndex = 2:length(selected);  % Never select delete bin
  end
  options.selectedIndex = selectedIndex;
  % Get the time markers
  if isfield(current_sort_info,'timeMarker')
    timeMarker = current_sort_info.timeMarker;
    timeMarker = convert_timeMarker_to_fileTime(timeMarker,handles);
    options.timeMarker = timeMarker;
  end
  % Now do the correlation computations and display the result
  cluster_correlations(shc,spikeClust,options);
  

% --- Executes on selection change in popupChannel.
function popupChannel_Callback(hObject, eventdata, handles)
% hObject    handle to popupChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupChannel
  status = cass_save_check(handles);
  if status
    cass_new_channel(handles);
  end
  
  
% --- Executes on button press in btnNextChannel.
function btnNextChannel_Callback(hObject, eventdata, handles)
% hObject    handle to btnNextChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  status = cass_save_check(handles);
  if status
    % Increment the channel index
    channelIndex = get(handles.popupChannel,'Value');
    channelIndex = min(channelIndex+1,...
      size(get(handles.popupChannel,'String'),1));
    set(handles.popupChannel,'Value',channelIndex);
    % Do everything needed for the new channel
    cass_new_channel(handles);
  end
  
  % if the "auto t corr" box is checked, automatically run it
  if getappdata(handles.autoTcorr,'autoTcorrOn')
      h = NaN;
      e = NaN;
      btnCorrelation_Callback(h, e, handles)
  end
  

  
function cass_new_channel(handles)  
  % Set the channel number
  % This has to be stored separately from the popup value, as when the
  % popup value changes there are some things we have to do to save data,
  % etc, before going on to the chosen channel.
  options = getappdata(handles.figCass,'options');
  currentChannel = options.channels(get(handles.popupChannel, ...
    'Value'));
  setappdata(handles.figCass,'currentChannel',currentChannel);
  % Clear all the old channel-specific data
  clear_channel_data(handles)
  % Do the needed processing to start a new channel
  filename = cass_choose_sort_info_file(handles);
  cass_load_sort_info(handles,filename)
  cass_load_snippets(handles)
  update_calculations(handles)
  % Show the data
  update_display(handles)


% --- Executes on button press in btnMerge.
function btnMerge_Callback(hObject, eventdata, handles)
% hObject    handle to btnMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  %update_selections(handles);
  merge_selected(handles)
  update_alphaR_range_calculations(handles)
  if getappdata(handles.auto_update_timeaxis,'autoUpdateOn')
    update_display(handles)
  else % cheap way to do all of update_display except for updating the time axis; 
       % probably better to set a flag, but lower risk of breaking things
       % this way...!
    update_cluster_axes_display(handles)
    update_projection_axis_display(handles)
    update_stability_axis_display(handles)
  end
  
  
  
% --- Executes on button press in btnSplit.
function btnSplit_Callback(hObject, eventdata, handles)
% hObject    handle to btnSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show the landmarks for just this cluster
axes_map = getappdata(handles.figCass,'axes_map');
selected = strmatch('on',{axes_map.selected},'exact');
if (length(selected) ~= 1)
  errordlg('You must select a single cluster to split')
end
%hlandmarks = findobj(handles.axesVisualizeShape,'Tag', ...
%		     'landmarkProjection');
set([axes_map.landmarkLineHandle],'Visible','off');
set(axes_map(selected).landmarkLineHandle,'Visible','on');
set(axes_map(selected).projLineHandle,'MarkerSize',1);
% Turn on splitting mode
handles.isSplitting = true;
guidata(hObject,handles);


% --- Executes on key press over figCass with no controls selected.
function figCass_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figCass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  currentChar = get(hObject,'CurrentCharacter');
  %double(currentChar)
  if strcmp(currentChar,char(8))
    % Merge selected axes with the noise
    %update_selections(handles);
    axes_map = getappdata(handles.figCass,'axes_map');
    axes_map(1).selected = 'on';
    setappdata(handles.figCass,'axes_map',axes_map);
    merge_selected(handles)
    update_alphaR_range_calculations(handles)
    update_display(handles)
%  elseif strcmp(currentChar,char(127))
%    % Mark the waveforms as needing deletion
%    update_selections(handles);
%    merge_selected(handles,-1);  % Put in the delete bin
%    update_alphaR_range_calculations(handles)
%    update_display(handles)
  end

  
% --- Executes on button press in btnExpand.
function btnExpand_Callback(hObject, eventdata, handles)
% hObject    handle to btnExpand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  %update_selections(handles)
  axes_map = getappdata(handles.figCass,'axes_map');
  clusterIndex = strmatch('on',{axes_map.selected},'exact')';
  if ~isempty(clusterIndex)
    % Open a new figure
    figure
    hfig_exp = gcf;
    options = getappdata(handles.figCass,'options');
    if isfield(options,'fig_positions')
     if isstruct(options.fig_positions)
       set(hfig_exp,'Position',options.fig_positions.expand);
     end
    end
    hold on
    cluster_axes_data = getappdata(handles.figCass,'cluster_axes_data');
    for tIndex = clusterIndex
      col = unique_color(tIndex,length(axes_map));
      % Show waveforms
      waveformsToShow = cluster_axes_data(tIndex).landmarkWaveforms;
      if ~isempty(waveformsToShow)
        plot(waveformsToShow,...
             'Color',col);
      end
    end
    hold off
    axis tight
  end
      


% --- Executes on button press in btnOpenSortInfo.
function btnOpenSortInfo_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenSortInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  options = getappdata(handles.figCass,'options');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  dirchan = [options.dirname filesep 'chan' num2str(currentChannel) ...
             filesep];
  filename = uigetfile([dirchan '*.mat'],'Pick a sort_info file:');
  if ~isequal(filename,0)
    % Clear all the old channel-specific data
    clear_channel_data(handles)
    % Start the new sort
    cass_load_sort_info(handles,filename)
    cass_load_snippets(handles)
    update_calculations(handles);
    update_display(handles);
  end
    

function cass_save_sort_info(handles)
  saveFile = getappdata(handles.figCass,'saveFile');
  if isempty(saveFile)
    return
  end
  if (strncmp(saveFile,'autosort_info',13) || strcmp(saveFile,'.mat'))
    errordlg('Save filename is messed up; going into debugger')
    keyboard
  end
  options = getappdata(handles.figCass,'options');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  dirchan = [options.dirname filesep 'chan' num2str(currentChannel) ...
             filesep];
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  % Put the text that's in the comment box into this structure
  current_sort_info.comment = get(handles.editComment,'String');
  % Put info about time markers into this structure
  timeMarkerTmp = getappdata(handles.figCass,'timeMarker');
  if ~isempty(timeMarkerTmp)
    % Convert relative time to file Index, scan #
    timeMarkerTmp = convert_timeMarker_to_fileTime(timeMarkerTmp, ...
						   handles);
    timeMarker = rmfield(timeMarkerTmp,'relativeTime');
    for i = 1:length(timeMarker)
      if ~isfield(timeMarker(i),'comment')
        timeMarker(i).comment = '';
      end
    end
    current_sort_info.timeMarker = timeMarker;
  end
  save_sort_info(current_sort_info,[dirchan saveFile]);
  setappdata(handles.figCass,'modified',0); % clear the modification flag
  
  
function cass_pick_save_file(handles)
  options = getappdata(handles.figCass,'options');
  currentChannel = getappdata(handles.figCass,'currentChannel');
  dirchan = [options.dirname filesep 'chan' num2str(currentChannel) ...
             filesep];
  filename = uiputfile([dirchan '*.mat'],'Pick a filename for saving:');
  if ~isequal(filename,0)
    setappdata(handles.figCass,'saveFile',filename);
  end


function status = cass_save_check(handles)
  % status output is 1 if we should proceed with going on to the next
  % channel perhaps after saving; status = 0 if user hits cancel.
  modified = getappdata(handles.figCass,'modified');
  
  % new wrinkle: if autoname is checked, ASSUME user wants to (a) save, and
  % (b) do so using whatever text is in the autoname text box; otherwise,
  % do what have always done....
  if modified
      if getappdata(handles.autoname,'useautosaveflag')
          status = 1;
          savename = get(handles.autoname_text,'String');
          if length(savename) < 1 % means nothing entered, so just go ahead and ask
              saveFile = getappdata(handles.figCass,'saveFile');
              if isempty(saveFile)
                  cass_pick_save_file(handles)
              end
          else
              setappdata(handles.figCass,'saveFile',savename);
          end
          cass_save_sort_info(handles)
      else
          button = questdlg('Do you want to save sorting results?');
          switch button
              case 'Yes'
                  saveFile = getappdata(handles.figCass,'saveFile');
                  if isempty(saveFile)
                      cass_pick_save_file(handles)
                  end
                  cass_save_sort_info(handles)
                  status = 1;
              case 'No'
                  status = 1;
              case 'Cancel'
                  status = 0;
          end
      end
  else
    status = 1;
  end
    
  
% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  saveFile = getappdata(handles.figCass,'saveFile');
  if isempty(saveFile)
    cass_pick_save_file(handles);
  end
  cass_save_sort_info(handles)

  
% --- Executes on button press in btnSaveAs.
function btnSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  cass_pick_save_file(handles)
  cass_save_sort_info(handles)


% --- Executes on selection change in popupXVis.
function popupXVis_Callback(hObject, eventdata, handles)
% hObject    handle to popupXVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupXVis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupXVis
  if isappdata(handles.axesVisualizeShape,'OXLim')
    rmappdata(handles.axesVisualizeShape,'OXLim')
  end
  update_projection_axis_display(handles)


% --- Executes on selection change in popupYVis.
function popupYVis_Callback(hObject, eventdata, handles)
% hObject    handle to popupYVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupYVis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupYVis
  if isappdata(handles.axesVisualizeShape,'OYLim')
    rmappdata(handles.axesVisualizeShape,'OYLim')
  end
  update_projection_axis_display(handles)


% --- Executes on button press in btnProjReset.
function btnProjReset_Callback(hObject, eventdata, handles)
% hObject    handle to btnProjReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  tmp_sort_info = select_current_sort_info(handles);
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  current_sort_info.projectDirections = tmp_sort_info.projectDirections;
  current_sort_info.t2V = tmp_sort_info.t2V;
  setappdata(handles.figCass,'current_sort_info',current_sort_info);
  % Recluster & redisplay
  set_projection_pulldowns(handles)
  update_projections(handles)
  update_landmarkIndex(handles)
  status = update_landmarkClust(handles);
  if status
    update_cluster_axes_calculations(handles)
    update_alphaR_range_calculations(handles)
  end
  update_display(handles)

  
% --- Executes on button press in btnProjChoose.
function btnProjChoose_Callback(hObject, eventdata, handles)
% hObject    handle to btnProjChoose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  nDirs = size(current_sort_info.projectDirections,2);
  % Give the user a choice of directions
  liststring = num2str((1:nDirs)');
  liststring(end+1,1) = 't';
  liststring = cellstr(liststring);
  selection = listdlg('ListString',liststring,...
                      'SelectionMode','multiple',...
                      'PromptString','Choose among projections:');
  if isempty(selection)
    return
  end
  % Set the chosen directions in current_sort_info
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  if (max(selection) <= nDirs)
    % User didn't choose time
    current_sort_info.t2V = 0;
  else
    % User chose time, pull this one off the selection list
    selection = setdiff(selection,nDirs+1);
  end
  current_sort_info.projectDirections = ...
      current_sort_info.projectDirections(:,selection);
  setappdata(handles.figCass,'current_sort_info',current_sort_info)
  set_projection_pulldowns(handles)
  update_projections(handles)
  update_landmarkIndex(handles)
  status = update_landmarkClust(handles);
  if status
    update_cluster_axes_calculations(handles)
    update_alphaR_range_calculations(handles)
  end
  update_display(handles)


% --- Executes on button press in btnViewProj.
function btnViewProj_Callback(hObject, eventdata, handles)
% hObject    handle to btnViewProj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  nDirs = size(current_sort_info.projectDirections,2);
  % Plot the projection directions
  figure
  dimx = ceil(sqrt(nDirs));
  dimy = ceil(nDirs/dimx);
  for i = 1:nDirs
    subplot(dimy,dimx,i)
    plot(current_sort_info.projectDirections(:,i));
    axis tight
  end


% --- Executes on button press in boxWatchClustering.
function boxWatchClustering_Callback(hObject, eventdata, handles)
% hObject    handle to boxWatchClustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boxWatchClustering


% --- Executes on button press in boxLockAspectRatio.
function boxLockAspectRatio_Callback(hObject, eventdata, handles)
% hObject    handle to boxLockAspectRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boxLockAspectRatio
  if get(handles.boxLockAspectRatio,'Value')
    set(handles.axesVisualizeShape,'DataAspectRatio',[1 1 1]);
  else
    set(handles.axesVisualizeShape,'DataAspectRatioMode','auto');
  end


% --- Executes on button press in invertSelectionsBtn.
function invertSelectionsBtn_Callback(hObject, eventdata, handles)
% hObject    handle to invertSelectionsBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  invert_selections(handles);


% --- Executes on button press in UndoBtn.
function UndoBtn_Callback(hObject, eventdata, handles)
% hObject    handle to UndoBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  old_sort_info = getappdata(handles.figCass,'old_sort_info');
  old_timeMarker = getappdata(handles.figCass,'old_timeMarker');
  store_current_sort_info(handles,old_sort_info,old_timeMarker);
  setappdata(handles.figCass,'old_sort_info',[]);
  setappdata(handles.figCass,'old_timeMarker',[]);
  set(handles.UndoBtn,'Enable','off')
  update_calculations_rest(handles)
  update_display(handles)

function set_further_processing(hObject, eventdata, handles, hax, action)
  % Determine the cluster that was clicked
  current_sort_info = getappdata(handles.figCass,'current_sort_info');
  axes_map=getappdata(handles.figCass, 'axes_map');
  axes_map_start_from=getappdata(handles.figCass, 'axes_map_start_from');
  virtual_axis = [axes_map.virtual_axis];
  n_visible_axes = size(handles.handle_matrix,2);
  clusterIndex = find(virtual_axis >= axes_map_start_from & ...
    virtual_axis < axes_map_start_from + ...
    n_visible_axes);
  thisAxisIndex = find(handles.handle_matrix(1,:) == hax);
  clIndex = clusterIndex(thisAxisIndex)-1; % Here's the cluster in axes_map
  % Set the flag for further processing
  current_sort_info.further_processing(clIndex) = true;
  setappdata(handles.figCass,'current_sort_info',current_sort_info);
  
  

  


% --- Executes on button press in autoname.
function autoname_Callback(hObject, eventdata, handles)
% hObject    handle to autoname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoname
setappdata(handles.autoname,'useautosaveflag',get(hObject,'Value'));
testpoint = 1;


function autoname_text_Callback(hObject, eventdata, handles)
% hObject    handle to autoname_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of autoname_text as text
%        str2double(get(hObject,'String')) returns contents of autoname_text as a double


% --- Executes during object creation, after setting all properties.
function autoname_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to autoname_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set(handles.autoname_text,'String','1')
%setappdata(handles.autoname,'useautosaveflag',0)
% GET OUT LATER WITH: get(handles.autoname_text,'String')


% --- Executes on button press in autoTcorr.
function autoTcorr_Callback(hObject, eventdata, handles)
% hObject    handle to autoTcorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoTcorr
setappdata(handles.autoTcorr,'autoTcorrOn',get(hObject,'Value'));


% --- Executes on selection change in popupEventRank.
function popupEventRank_Callback(hObject, eventdata, handles)
% hObject    handle to popupEventRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupEventRank contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupEventRank
  update_projection_axis_display(handles);

% --- Executes during object creation, after setting all properties.
function popupEventRank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupEventRank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_update_timeaxis.
function auto_update_timeaxis_Callback(hObject, eventdata, handles)
% hObject    handle to auto_update_timeaxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of auto_update_timeaxis
setappdata(handles.auto_update_timeaxis,'autoUpdateOn',get(hObject,'Value'));

% --- Executes on button press in update_time_asix_plot.
function update_time_asix_plot_Callback(hObject, eventdata, handles)
% hObject    handle to update_time_asix_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_time_axis_display(handles)

