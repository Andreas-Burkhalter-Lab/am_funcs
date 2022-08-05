function varargout = mdcluster(varargin)
% MDCLUSTER: multidimensional clustering GUI
% Syntax:
%   labels = mdcluster(data)
%   labels = mdcluster(data,options)
% where
%   data is either a d-by-N matrix, where each column is one data point in
%   d dimensions, or a structure with the following fields:
%     points: a d-by-N matrix
%     labels (optional): a 1-by-N input vector of initial cluster numbers
%     time (optional): a 1-by-N vector, each data point is labeled with a
%       "time" of occurrence (useful if the characteristics of your clusters
%       might change over time or any other 1-dimensional parameter).
%     timerange (optional): a 2-vector, supplied in case you want to use a
%       consistent range for the time axis.
%     corrTMax (optional): the maximum time (in the same units as the time
%       vector) used in plotting correlations.
%   options is an (optional) structure of the format described in
%     MDCLUSTER_OPTIONS. Perhaps the best plan is to execute
%        options = mdcluster_options
%     and then to examine/add to/modify the options structure.
%
%     In addition to the settings described in MDCLUSTER_OPTIONS, you may
%     specify the default value of any GUI element as additional fields,
%     where the "tag" of the GUI element is the name of the field, and its
%     value is assigned to the "Value" property of that GUI element on
%     startup.  GUI element tags can be determined by running GUIDE on
%     mdcluster and examining tags with the Property Browser.
%
% On output,
%   labels is a 1-by-N vector giving the cluster number assigned to each
%     data point.
%
% This function benefits from having the following auxillary software installed:
%    Dimensionality Reduction Toolbox
%    FastICA
%
% Example:
%   % Create a data set with 3 clusters that can be seen by inspecting
%   % coordinates 3 and 4 (but it's a slightly tricky data set because
%   % coordinates 1 and 2 have the largest variance, so PCA and k-means
%   % will fail)
%   x = randn(4,1200); g2 = 401:800; g3 = 801:1200; x(1:2,:) = 20*x(1:2,:); x(3,g2) = x(3,g2)+30; x(4,g3) = x(4,g3)+30;
%   clust = mdcluster(x);
%
% See also: MDCLUSTER_OPTIONS.
  
% Copyright 2008 by Timothy E. Holy
  
%        
% MDCLUSTER M-file for mdcluster.fig
%      MDCLUSTER, by itself, creates a new MDCLUSTER or raises the existing
%      singleton*.
%
%      H = MDCLUSTER returns the handle to a new MDCLUSTER or the handle to
%      the existing singleton*.
%
%      MDCLUSTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MDCLUSTER.M with the given input arguments.
%
%      MDCLUSTER('Property','Value',...) creates a new MDCLUSTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mdcluster_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mdcluster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mdcluster

% Last Modified by GUIDE v2.5 16-Jul-2008 12:06:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mdcluster_OpeningFcn, ...
                   'gui_OutputFcn',  @mdcluster_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before mdcluster is made visible.
function mdcluster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mdcluster (see VARARGIN)

% Choose default command line output for mdcluster
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if (length(varargin) < 1)
  error('You have to at least supply data points to start mdclsuter');
end
data = varargin{1};
if ~isstruct(data)
  tmp = data;
  clear data;
  data.points = tmp;
end
if (length(varargin) > 1)
  options = varargin{2};
else
  options = mdcluster_options;
end
setappdata(hObject,'data',data);
sz = size(data.points);
textInfo = sprintf('%d points, %d dimensions',sz(2),sz(1));
setappdata(hObject,'textInfoBase',textInfo);

if isfield(data,'labels')
  labels = data.labels;
  set(handles.btnRevert,'Enable','on','Visible','on');
else
  labels = [];
end
setappdata(hObject,'labels',labels);
if ~isfield(data,'time')
  set(handles.checkboxShowTime,'Visible','off');
  set(handles.checkboxUseTime,'Visible','off');
  set(handles.btnTCorrelations,'Visible','off');
  set(handles.editt2V,'Visible','off');
  set(handles.textt2V,'Visible','off');
end

set(handles.popupmenuProjectionMethods,'String',{''});
set(handles.popupmenuClusterMethods,'String',{''});

setappdata(hObject,'options',options);

% Define the usable methods
handles = mdcluster_initialize(handles,options);

set(handles.popupmenuX,'Value',2);
setappdata(handles.figure_mdcluster,'markerSizes',options.markerSizes);

% Set the projection method
indx = strmatch(options.defaultProjectionMethod,get(handles.popupmenuProjectionMethods,'String'),'exact');
if ~isempty(indx) && indx ~= get(handles.popupmenuProjectionMethods,'Value')
  set(handles.popupmenuProjectionMethods,'Value',indx);
  handles = mdcluster_render_parameters(handles,options,'projectionMethods');
end

guidata(hObject, handles);
mdcluster_project(handles);

% Allow clicking to select clusters
set(handles.axesProjection,'ButtonDownFcn',@mdcluster_click);

% UIWAIT makes mdcluster wait for user response (see UIRESUME)
if isscalar(options.blocking)
  if options.blocking
    uiwait(handles.figure_mdcluster);
  end
end


% --- Outputs from this function are returned to the command line.
function varargout = mdcluster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ishandle(hObject)
  varargout{1} = getappdata(hObject,'labels');
  delete(hObject);
else
  varargout{1} = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Display code                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mdcluster_update_display(handles)
  proj = getappdata(handles.figure_mdcluster,'proj');
  options = getappdata(handles.figure_mdcluster,'options');
  coordIndex = [get(handles.popupmenuX,'Value') get(handles.popupmenuY,'Value')];
  coordIndex = min([coordIndex; [1 1]*size(proj,1)],[],1);
  set(handles.popupmenuX,'Value',coordIndex(1));
  set(handles.popupmenuY,'Value',coordIndex(2));
  y = proj(coordIndex(2),:);
  showingTime = false;
  if get(handles.checkboxShowTime,'Value')
    data = getappdata(handles.figure_mdcluster,'data');
    x = data.time;
    showingTime = true;
  else
    x = proj(coordIndex(1),:);
  end
  labels = getappdata(handles.figure_mdcluster,'labels');
  if isempty(labels)
    clabel = {1:length(x)};
  else
    clabel = agglabel(labels);
  end
  n_clusters = length(clabel);
  hax = handles.axesProjection;
  cla(hax);
  oneDMode = false;
  if isPlottingDensity(handles)
    % Density mode
    sliderVal = get(handles.sliderDensity,'Value');
    nbins = 30*(1-sliderVal) + 300*sliderVal;
    nbins = round(nbins);
    if (size(proj,1) > 1 || showingTime)
      % Two-dimensional display
      rect = [min(x) max(x) min(y) max(y)];
      w = 1.2;
      rect(1:2) = widenrange(rect(1:2),w);
      rect(3:4) = widenrange(rect(3:4),w);
      [n,xC,yC] = hist2d(x,y,rect,nbins,nbins);
      imagesc(xC,yC,log(n+1)','Parent',hax);
      set(hax,'YDir','normal');
      colormap(1-gray);
      setappdata(handles.figure_mdcluster,'clusterLineH',[]);
    else
      % One-dimensional mode
      xr = [min(proj) max(proj)];
      w = 1.2;
      xr = widenrange(xr,w);
      xc = linspace(xr(1),xr(2),nbins);
      for i = 1:n_clusters
        col = unique_color(i,n_clusters);
        n = hist(proj(clabel{i}),xc);
        hb = bar(xc,n,1);
        hold on
        set(hb,'FaceColor',col,'EdgeColor','none');
      end
      oneDMode = true;
      set(hax,'XLim',xr);
    end
  else
    % Points mode
    markerSizes = getappdata(handles.figure_mdcluster,'markerSizes');
    hline = nan(1,n_clusters);
    for i = 1:n_clusters
      col = unique_color(i,n_clusters);
      hline(i) = line(x(clabel{i}),y(clabel{i}),'LineStyle','none','Marker','.','Color',col,'Parent',hax,'HitTest','off','SelectionHighlight','off','MarkerSize',markerSizes(1));
    end
    setappdata(handles.figure_mdcluster,'clusterLineH',hline);
    set(hax,'XLimMode','auto','YLimMode','auto');
    if options.axis_equal
      axis(handles.axesProjection,'equal');
    end
  end

  if showingTime
    xlabel('Time')
    if isfield(data,'timerange')
      set(hax,'XLim',data.timerange)
    end
  else
    xlabel(['Dim. ' num2str(coordIndex(1))])
  end
  if oneDMode
    ylabel('Counts');
  else
    ylabel(['Dim. ' num2str(coordIndex(2))]);
  end
  
  textInfo = getappdata(handles.figure_mdcluster,'textInfoBase');
  if ~isempty(labels)
    textInfo = [textInfo sprintf('; %d clusters',n_clusters)];
  end
  set(handles.textInfo,'String',textInfo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Method execution code                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = mdcluster_prepare_s(parameterList)
% A utility function to evaluate parameter lists and prepare the
% parameter/option structure
  s = struct;   % in case parameterList is empty
  for i = 1:length(parameterList)
    name = parameterList(i).name;
    % Change the name to be a valid fieldname
    name = strrep(name,'#','n');
    name = strrep(name,' ','_');
    name = strrep(name,'.','');
    value = parameterList(i).value;
%     if ~isempty(parameterList(i).castfcn)
%       value = parameterList(i).castfcn(parameterList(i).value);
%     else
%       value = parameterList(i).value;
%     end
    s.(name) = value;
  end

function output = mdcluster_method_executor(handles,fieldstruct)
% Gather data and run the "execute" function for a method
  options = getappdata(handles.figure_mdcluster,'options');
  % Get updated labels (required by a few methods, like LDA)
  labels = getappdata(handles.figure_mdcluster,'labels');
  data = getappdata(handles.figure_mdcluster,'data');
  % For cluster methods, we have to decide whether to use projected data
  % or raw data, and we might have to append time
  if strcmp(fieldstruct(1:7),'cluster')
    if (get(handles.checkboxRawData,'Value'))
      x = data.points;
    else
      x = getappdata(handles.figure_mdcluster,'proj');
    end
    if (get(handles.checkboxUseTime,'Value') && isfield(data,'time'))
      t2V = str2num(get(handles.editt2V,'String'));
      x = [x; t2V*data.time];
    end
  else
    % For projection, use the raw points
    x = data.points;
  end
  methodIndex = mdcluster_get_methodIndex(handles,options,fieldstruct);
  s = mdcluster_prepare_s(options.(fieldstruct)(methodIndex).parameterList);
  s.labels = labels;
  hline = getappdata(handles.figure_mdcluster,'clusterLineH');
  s.selected = mdcluster_get_selected(hline);
  set(handles.textStatus,'String','Busy')
  drawnow
  output = options.(fieldstruct)(methodIndex).execute(x,s);
  set(handles.textStatus,'String','')


function mdcluster_project(handles)
% Calculate projections of data points
  proj = mdcluster_method_executor(handles,'projectionMethods');
  setappdata(handles.figure_mdcluster,'proj',proj);
  ndims = size(proj,1);
  str = int2str((1:ndims)');
  cstr = cellstr(str)';
  set(handles.popupmenuY,'String',cstr);
  set(handles.popupmenuX,'String',cstr);
  % If the old setting was 1-dimensional, and the new one isn't, then by
  % default go back to 2d
  if (ndims > 1)
    if (get(handles.popupmenuY,'Value') == 1 && get(handles.popupmenuX,'Value') == 1)
      set(handles.popupmenuX,'Value',2);
    end
  else
    % Set on density display no matter what
    set(handles.checkboxDensityHistogram,'Value',1);
    set(handles.sliderDensity,'Visible','on');
  end
  mdcluster_update_display(handles);


function mdcluster_cluster(handles)
% Cluster the data points
  oldlabels = getappdata(handles.figure_mdcluster,'labels');
  labels = mdcluster_method_executor(handles,'clusterMethods');
  setappdata(handles.figure_mdcluster,'labels',labels);
  mdcluster_update_display(handles);
  if (isempty(oldlabels) || max(oldlabels) == 1)
    % Since we only now have clusters, some new projection methods might be
    % available
    mdcluster_create_methods(handles,'projectionMethods');
  end

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Utilities                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [methodIndex,fieldgui] = mdcluster_get_methodIndex(handles,options,fieldstruct)
% Utility function: determine the index into the methods structure array
% of the currently-selected method
  % Mangle the names to get the tag of the object
  fieldgui = fieldstruct;
  fieldgui(1) = upper(fieldgui(1));
  fieldgui = ['popupmenu' fieldgui];
  % Determine the current setting and the corresponding entry in options
  str = get(handles.(fieldgui),'String');
  setting = str{get(handles.(fieldgui),'Value')};
  methodIndex = strmatch(setting,{options.(fieldstruct).name},'exact');
    
function isselected = mdcluster_get_selected(hline)
  selected = get(hline,'Selected');
  if ~iscell(selected)
    selected = {selected};
  end
  isselected = cellfun(@(s) strcmp(s,'on'),selected);
  
function mdcluster_set_selected(hfig,hline,isselected)
  markerSizes = getappdata(hfig,'markerSizes');
  set(hline(isselected),'MarkerSize',markerSizes(2));
  set(hline(~isselected),'MarkerSize',markerSizes(1));
  selected(isselected) = {'on'};
  selected(~isselected) = {'off'};
  set(hline,{'Selected'},selected(:));
  if (sum(isselected) == 1)
    % Make sure selected cluster is on top
    hax = get(hline(1),'Parent');
    set(hax,'Children',[hline(isselected) hline(~isselected)]);
  end

  
function handled = mdcluster_click(sender,event)
  hax = sender;
  hfig = get_parent_fig(hax);
  cp = get(hax,'CurrentPoint');
  cp = cp(1,1:2);
  [pos,hline_click] = findpoint(cp,hax);
  selectionType = get(hfig,'SelectionType');
  hline = getappdata(hfig,'clusterLineH');
  if isempty(hline)
    return  % Can't do anything if don't have line (e.g., in density mode)
  end
  % Check to see if the distance is so large that we should unselect all
  % Make this measurement in normalized units
  [xnclick,ynclick] = data2norm(cp(1),cp(2),hax);
  [xnclosest,ynclosest] = data2norm(pos(1),pos(2),hax);
  if ((xnclick-xnclosest)^2 + (ynclick-ynclosest)^2 > 0.02^2 && strcmp(selectionType,'normal'))
    isselected = false(length(hline),1);
  else
    % Close enough. Either select it or toggle the state
    clickIndex = find(hline == hline_click);
    isselected = mdcluster_get_selected(hline);
    switch selectionType
      case 'normal'
        isselected = false(length(hline),1);
        isselected(clickIndex) = true;
      case 'extend'
        isselected(clickIndex) = ~isselected(clickIndex);
    end
  end
  mdcluster_set_selected(hfig,hline,isselected);  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Initializiation & GUI state code                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = mdcluster_initialize(handles,options)
  handles = mdcluster_create_methods(handles,'projectionMethods');
  handles = mdcluster_create_methods(handles,'clusterMethods');
  mdcluster_set_gui_elements(handles,options);
  
function mdcluster_set_gui_elements(handles,options)
% Sets the states of GUI elements from the input options structure
  fn_options = fieldnames(options);
  fn_handles = fieldnames(handles);
  c = intersect(fn_options,fn_handles);
  for i = 1:length(c)
    set(handles.(c{i}),'Value',options.(c{i}));
  end
  if ~get(handles.checkboxDensityHistogram,'Value')
    set(handles.sliderDensity,'Visible','off');
  end
  states = {'off','on'};
  state = get(handles.checkboxShowTime,'Value');
  set(handles.textTimeX,'Visible',states{state+1});
  set(handles.popupmenuX,'Visible',states{2-state});

function handles = mdcluster_create_methods(handles,fieldstruct)
% Determines which methods are available and fills the popupmenus
% accordingly
% This is called on initialization but also any time an important
% dependency changes, e.g., if labels become available then methods that
% specify haveClusters suddenly become available
%
  % Mangle the names to get the tag of the object
  fieldgui = fieldstruct;
  fieldgui(1) = upper(fieldgui(1));
  fieldgui = ['popupmenu' fieldgui];
  % Store the old setting of the popup menu
  str = get(handles.(fieldgui),'String');
  oldsetting = str{get(handles.(fieldgui),'Value')};
  options = getappdata(handles.figure_mdcluster,'options');
  % Determine which methods are available & pass the valid test
  labels = getappdata(handles.figure_mdcluster,'labels');
  haveClusters = false;
  if (~isempty(labels) && max(labels) > 1)
    haveClusters = true;
  end
  n_methods = length(options.(fieldstruct));
  method_names = {options.(fieldstruct).name};
  valid = false(1,n_methods);
  for methodIndex = 1:n_methods
    if ischar(options.(fieldstruct)(methodIndex).valid)
      if strcmp(options.(fieldstruct)(methodIndex).valid,'haveClusters')
        valid(methodIndex) = haveClusters;
      end
    else
      valid(methodIndex) = ...
        options.(fieldstruct)(methodIndex).valid;
    end
  end
  method_names = method_names(valid);
  % Set the list of available methods, and update the popupmenu's "Value"
  % to point to the current one
  set(handles.(fieldgui),'String',method_names);
  value = strmatch(oldsetting,method_names,'exact');
  reset_parameters = false;
  if isempty(value)
    % The previous setting is not available, pick the first one
    value = 1;
    reset_parameters = true;
  end
  set(handles.(fieldgui),'Value',value);
  if reset_parameters
    handles = mdcluster_render_parameters(handles,options,fieldstruct);
  end
  
function handles = mdcluster_render_parameters(handles,options,fieldstruct)
% This function sets up the text and edit boxes corresponding to
% parameterList for each method. It also defines a callback function that
% executes to update the options structure to store changes.
  fieldhandles = [fieldstruct 'Handles'];
  [methodIndex,fieldgui] = mdcluster_get_methodIndex(handles,options,fieldstruct);
  thisMethod = options.(fieldstruct)(methodIndex);
  % Delete any current parameter GUI elements
  if isfield(handles,fieldhandles)
    for i = 1:length(handles.(fieldhandles))
      if ishandle(handles.(fieldhandles)(i))
        delete(handles.(fieldhandles)(i));
      end
    end
  end
  % Calculate the placement
  units = get(handles.(fieldgui),'Units');
  set(handles.(fieldgui),'Units','pixels');
  pos = get(handles.(fieldgui),'Position');
  set(handles.(fieldgui),'Units',units);
  pos(2) = pos(2)-30;
  % Draw the GUI objects
  handles.(fieldhandles) = [];
  for i = 1:length(thisMethod.parameterList)
    celem = thisMethod.parameterList(i);
    htxt = uicontrol('Units','pixels','Position',pos,'Style','text','String',celem.name);
    extnt = get(htxt,'Extent');
    pos(3:4) = extnt(3:4);
    set(htxt,'Position',pos);
    pos(1) = pos(1)+pos(3)+10;
    pos(3) = 40;
    cbdata = {fieldstruct,methodIndex,i};
    hedit = uicontrol('Units','pixels','Position',pos,'Style','edit', ...
		      'String',celem.value,'BackgroundColor','white',...
          'Callback',{@mdcluster_set_parameter,cbdata});
    pos(1) = pos(1)+pos(3)+20;
    handles.(fieldhandles)(2*i-1:2*i) = [htxt hedit];  % Save the handles
  end
  set(handles.(fieldhandles),'Units','normalized');  % make resizable

function mdcluster_set_parameter(hObject,eventData,cbdata)
% This is the callback function that executes to update the options
% structure to store changes (set up in mdcluster_render_parameters)
  handles = guidata(hObject);
  options = getappdata(handles.figure_mdcluster,'options');
  value = get(hObject,'String');
  % Validate the input
  castfcn = ...
      options.(cbdata{1})(cbdata{2}).parameterList(cbdata{3}).castfcn;
  if ~isempty(castfcn)
    value = castfcn(value);
  end
  % Save the result
  options.(cbdata{1})(cbdata{2}).parameterList(cbdata{3}).value = value;
  setappdata(handles.figure_mdcluster,'options',options);
  % Put the validated version into the edit box
  if ~ischar(value)
    value = num2str(value);
  end
  set(hObject,'String',value);

function xy = widenrange(xyi,w)
xy = w/2*diff(xyi)*[-1 1]+mean(xyi);

function iPD = isPlottingDensity(handles)
iPD = get(handles.checkboxDensityHistogram,'Value');
  
  
% --- Executes on selection change in popupmenuProjectionMethods.
function popupmenuProjectionMethods_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectionMethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuProjectionMethods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuProjectionMethods
options = getappdata(handles.figure_mdcluster,'options');
handles = mdcluster_render_parameters(handles,options,'projectionMethods');
guidata(handles.figure_mdcluster,handles);

% --- Executes during object creation, after setting all properties.
function popupmenuProjectionMethods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuProjectionMethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxShowTime.
function checkboxShowTime_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Toggle the visibility of the X-coordinate popupmenu
states = {'on','off'};
state = get(hObject,'Value');
set(handles.popupmenuX,'Visible',states{state+1});
set(handles.textTimeX,'Visible',states{2-state});
% Redisplay
mdcluster_update_display(handles);


function editt2V_Callback(hObject, eventdata, handles)
% hObject    handle to editt2V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editt2V as text
%        str2double(get(hObject,'String')) returns contents of editt2V as a double


% --- Executes during object creation, after setting all properties.
function editt2V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editt2V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxDensityHistogram.
function checkboxDensityHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDensityHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDensityHistogram
if get(hObject,'Value')
  set(handles.sliderDensity,'Visible','on')
else
  set(handles.sliderDensity,'Visible','off')
end
mdcluster_update_display(handles);
if ~isPlottingDensity(handles)
  set(handles.axesProjection,'ButtonDownFcn',@mdcluster_click);
end

% --- Executes on button press in checkboxUseTime.
function checkboxUseTime_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUseTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUseTime


% --- Executes on selection change in popupmenuClusterMethods.
function popupmenuClusterMethods_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuClusterMethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the text on the Cluster button appropriately
val = get(handles.popupmenuClusterMethods,'Value');
method = get(handles.popupmenuClusterMethods,'String');
method = method{val};
if strcmp(method,'Manual merge')
  set(handles.buttonCluster,'String','Merge');
  set(handles.buttonCycle,'Visible','on');
else
  set(handles.buttonCluster,'String','Cluster');
  set(handles.buttonCycle,'Visible','off');
end
options = getappdata(handles.figure_mdcluster,'options');
handles = mdcluster_render_parameters(handles,options,'clusterMethods');
guidata(handles.figure_mdcluster,handles);


% --- Executes during object creation, after setting all properties.
function popupmenuClusterMethods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuClusterMethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxRawData.
function checkboxRawData_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRawData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxRawData


% --- Executes on button press in buttonCancel.
function buttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to buttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure_mdcluster);

% --- Executes on button press in buttonDone.
function buttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to buttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure_mdcluster);

% --- Executes on selection change in popupmenuY.
function popupmenuY_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuY contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuY
mdcluster_update_display(handles);


% --- Executes during object creation, after setting all properties.
function popupmenuY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuX.
function popupmenuX_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenuX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuX
mdcluster_update_display(handles);

% --- Executes during object creation, after setting all properties.
function popupmenuX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonProject.
function buttonProject_Callback(hObject, eventdata, handles)
% hObject    handle to buttonProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  mdcluster_project(handles);
  
% --- Executes on button press in buttonCluster.
function buttonCluster_Callback(hObject, eventdata, handles)
% hObject    handle to buttonCluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  mdcluster_cluster(handles);

% --- Executes on slider movement.
function sliderDensity_Callback(hObject, eventdata, handles)
% hObject    handle to sliderDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
mdcluster_update_display(handles);


% --- Executes during object creation, after setting all properties.
function sliderDensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btnRevert.
function btnRevert_Callback(hObject, eventdata, handles)
% hObject    handle to btnRevert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(handles.figure_mdcluster,'data');
labels = data.labels;
setappdata(handles.figure_mdcluster,'labels',labels);
mdcluster_update_display(handles);


% --- Executes on button press in buttonCycle.
function buttonCycle_Callback(hObject, eventdata, handles)
% hObject    handle to buttonCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hline = getappdata(handles.figure_mdcluster,'clusterLineH');
isselected = mdcluster_get_selected(hline);
if (sum(isselected) > 1)
  isselected(:) = false;
end
if (sum(isselected) == 0)
  isselected(1) = true;
else
  selectedIndex = find(isselected);
  isselected(selectedIndex) = false;
  selectedIndex = selectedIndex+1;
  if (selectedIndex > length(isselected))
    selectedIndex = 1;
  end
  isselected(selectedIndex) = true;
end
mdcluster_set_selected(handles.figure_mdcluster,hline,isselected);

% --- Executes on button press in btnTCorrelations.
function btnTCorrelations_Callback(hObject, eventdata, handles)
% hObject    handle to btnTCorrelations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = getappdata(handles.figure_mdcluster,'data');
labels = getappdata(handles.figure_mdcluster,'labels');
if isempty(labels)
  clabel = {1:length(data.time)};
else
  clabel = agglabel(labels);
end
ops.trange = data.corrTMax;
n_clusters = length(clabel);
for i = 1:n_clusters
  ops.col(i,:) = unique_color(i,n_clusters);
end
if ~isPlottingDensity(handles)
  hline = getappdata(handles.figure_mdcluster,'clusterLineH');
  selected = mdcluster_get_selected(hline);
  if any(selected)
    clabel = clabel(selected);
    ops.col = ops.col(selected,:);
  end
end
n_clusters = length(clabel);
t = cell(1,n_clusters);
for i = 1:n_clusters
  t{i} = data.time(clabel{i});
  if ~issorted(t{i})
    t{i} = sort(t{i});
  end
end
set(handles.textStatus,'String','Busy')
drawnow
plot_all_corr(t,ops);
set(handles.textStatus,'String','')
