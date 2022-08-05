%function add_landmarks(sortfile_in,sortfile_out,replicate)
function varargout = add_landmarks(varargin)
% ADD_LANDMARKS M-file for add_landmarks.fig
%      ADD_LANDMARKS, by itself, creates a new ADD_LANDMARKS or raises the existing
%      singleton*.
%
%      H = ADD_LANDMARKS returns the handle to a new ADD_LANDMARKS or the handle to
%      the existing singleton*.
%
%      ADD_LANDMARKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADD_LANDMARKS.M with the given input arguments.
%
%      ADD_LANDMARKS('Property','Value',...) creates a new ADD_LANDMARKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before add_landmarks_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to add_landmarks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help add_landmarks

% Last Modified by GUIDE v2.5 20-Feb-2006 22:16:09
% 2006-04-17 RCH: altered save function so it leaves you in the same
%      directory you started in

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @add_landmarks_OpeningFcn, ...
                   'gui_OutputFcn',  @add_landmarks_OutputFcn, ...
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


% --- Executes just before add_landmarks is made visible.
function add_landmarks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to add_landmarks (see VARARGIN)

% Choose default command line output for add_landmarks
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

install_mouse_event_handler(handles.axesVisualizeShape,'down',...
                            @zoom_mouse_down)
%install_mouse_event_handler(handles.axesVisualizeShape,'up',...
%                            @zoom_mouse_up)
%install_mouse_event_handler(handles.axesVisualizeShape,'move',...
%                            @zoom_mouse_move)


% UIWAIT makes add_landmarks wait for user response (see UIRESUME)
% uiwait(handles.figAL);


% --- Outputs from this function are returned to the command line.
function varargout = add_landmarks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function add_landmarks_open(handles,path_in,sortfile_in,replicate)
  if (nargin < 4)
    replicate = 1;
  end
  load([path_in sortfile_in])
  current_sort_info = sort_info(replicate);
  load([path_in filesep '../overview'])
  shc = sortheader_importchan(sorthead,sort_info(1).channel);
  nFiles = length(shc);
  nsnips_per_file = [shc.numofsnips];
  nsnips_total = sum(nsnips_per_file);
  sniprange = reshape([shc.sniprange],2,length(shc));
  sniprange = unique(sniprange','rows');
  if (size(sniprange,1) > 1)
    error(['Snippets do not have the same number of samples in each ' ...
           'file']);
  end
  snipLength = diff(sniprange);
  max_snips_in_mem = min(round(options.max_snip_memsize/snipLength),...
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
  tstart = sortheader_absolute_starttime(sorthead);
  absoluteSnipTime = cell(1,length(sorthead));
  snip = cell(1,length(sorthead));
  for i = 1:nFiles
    % For user feedback, convert sniptimes to absolute timescale
    sniptimes = double(shc(i).sniptimes(snipIndex{i})');
    absoluteSnipTime{i} = sniptimes/sorthead(i).scanrate ...
        + tstart(i) - min(tstart);
    if size(absoluteSnipTime{i},1) == 0
      absoluteSnipTime{i} = (absoluteSnipTime{i})';
    end
    % Read in snippets
    snip(i) = sortheader_readsnips(shc(i),snipIndex(i));
  end
  % Calculate projections
  snip = cat(2,snip{:});
  if current_sort_info.use_projection
    snipProj = current_sort_info.projectDirections'*snip;
    landmarkWaveformProj = ...
      current_sort_info.projectDirections'*current_sort_info.landmarkWaveform;
  end
  % Add time as a sorting coordinate
  sniptime = cat(2,absoluteSnipTime{:});
  snipProj = [snipProj; sniptime*current_sort_info.t2V];
  landmarkWaveformProj = [landmarkWaveformProj; ...
              current_sort_info.landmarkT*current_sort_info.t2V];
  setappdata(handles.figAL,'snipProj',snipProj)
  setappdata(handles.figAL,'landmarkWaveformProj',landmarkWaveformProj)
  setappdata(handles.figAL,'current_sort_info',current_sort_info)
  setappdata(handles.figAL,'saveDir',path_in)
  setappdata(handles.figAL,'snip',snip)
  setappdata(handles.figAL,'sniptime',sniptime)
  
  n_spatial_directions = size(current_sort_info.projectDirections,2);
  strPulldown = num2str((1:n_spatial_directions)');
  if current_sort_info.t2V
    strPulldown(end+1,1) = 't';  % add time as a direction
  end
  set(handles.popupXVis,'String',strPulldown,'Value',size(snipProj,1));
  set(handles.popupYVis,'String',strPulldown,'Value',1);

  update_landmarkIndex(handles)
  
  add_landmarks_plot(handles)

function result=zoom_mouse_down(sender, event_args)
  %setappdata(sender, 'isMouseDown', 1);
  currentPoint = get(sender,'CurrentPoint');
  currentPoint = currentPoint(1,[1 2]);
  result=1;
  %setappdata(sender,'startPoint',currentPoint(1,[1 2]));
  % Find the closest point to this point in the 2-d projection
  handles=guidata(sender);
  snipProj = getappdata(handles.figAL,'snipProj');
  x_coord = get(handles.popupXVis,'Value');
  y_coord = get(handles.popupYVis,'Value');
  [dist,index] = mindist(currentPoint',snipProj([x_coord y_coord],:));
  newProj = snipProj(:,index);
  % Now shift this new point by a kmeans-type algorithm in the space of
  % full dimensionality
  landmarkWaveformProj = getappdata(handles.figAL,'landmarkWaveformProj');
  n_landmarks = size(landmarkWaveformProj,2)+1;
  is_moved = 1;
  snipIndexOld = [];
  while is_moved
    [dist,index_closest] = mindist(snipProj,[landmarkWaveformProj newProj]);
    snipIndex = find(index_closest == n_landmarks);
    if (~isequal(snipIndex,snipIndexOld) && ~isempty(snipIndex))
      snipIndexOld = snipIndex;
      newProj = mean(snipProj(:,snipIndex),2);
    else
      is_moved = 0;
    end
  end
  if ~isempty(snipIndex)
    snip = getappdata(handles.figAL,'snip');
    newlandmarkWaveform = mean(snip(:,snipIndex),2);
    sniptime = getappdata(handles.figAL,'sniptime');
    newlandmarkT = mean(sniptime(:,snipIndex));
    current_sort_info = getappdata(handles.figAL,'current_sort_info');
    current_sort_info.landmarkWaveform = ...
        [current_sort_info.landmarkWaveform newlandmarkWaveform];
    current_sort_info.landmarkT = ...
        [current_sort_info.landmarkT newlandmarkT];
    lmc = max(current_sort_info.landmarkClust,[],2) + 1;
    current_sort_info.landmarkClust = ...
        [current_sort_info.landmarkClust, lmc];
    landmarkWaveformProj = ...
      [current_sort_info.projectDirections'*current_sort_info.landmarkWaveform;
      current_sort_info.landmarkT*current_sort_info.t2V];
    setappdata(handles.figAL,'current_sort_info',current_sort_info);
    setappdata(handles.figAL,'landmarkWaveformProj',landmarkWaveformProj)
    update_landmarkIndex(handles)
    add_landmarks_plot(handles)
  end
    
  
% function result=zoom_mouse_move(sender, event_args)
%   if(getappdata(sender, 'isMouseDown'))
%     setappdata(sender, 'isDragging', 1);
%     currentPoint = get(sender,'CurrentPoint');
%     currentPoint = currentPoint(1,[1 2]);
%     startPoint = getappdata(sender,'startPoint');
%     xdata = [startPoint(1) currentPoint([1 1]) startPoint([1 1])];
%     ydata = [startPoint([2 2]) currentPoint([2 2]) startPoint(2)];
%     if ~isappdata(sender,'hrbbox')
%       hrbbox = line(xdata,ydata,...
%                     'LineStyle','--',...
%                     'Color','k',...
%                     'EraseMode','xor',...
%                     'Tag','rbbox');
%       setappdata(sender,'hrbbox',hrbbox);
%     else
%       hrbbox = getappdata(sender,'hrbbox');
%       set(hrbbox,'XData',xdata,'YData',ydata);
%     end
%   end
%   result=0;
%    
% function result=zoom_mouse_up(sender, event_args)
%   handles=guidata(sender);
%   currentPoint = get(sender,'CurrentPoint');
%   currentPoint = currentPoint(1,[1 2]);
%   startPoint = getappdata(sender,'startPoint');
%   % Clean up the mouse stuff
%   setappdata(sender, 'isMouseDown', 0);
%   setappdata(sender, 'isDragging', 0);
%   if isappdata(sender,'hrbbox')
%     hrbbox = getappdata(sender,'hrbbox');
%     delete(hrbbox);
%     rmappdata(sender,'hrbbox');
%   end
%   result=1;
%   % Find the points inside the rectangle
%   p1 = min(startPoint,currentPoint);
%   offset = abs(startPoint - currentPoint);
%   x_coord = get(handles.popupXVis,'Value');
%   y_coord = get(handles.popupYVis,'Value');
%   snipProj = getappdata(handles.figAL,'snipProj');
%   snipIndex = find(snipProj(x_coord,:) > p1(1) & ...
%                    snipProj(x_coord,:) < p1(1)+offset(1) & ...
%                    snipProj(y_coord,:) > p1(2) & ...
%                    snipProj(y_coord,:) < p1(2)+offset(2));
%   if ~isempty(snipIndex)
%     % Introduce a new landmark at the mean of the snippets in the rectangle
%     snip = getappdata(handles.figAL,'snip');
%     newlandmarkWaveform = mean(snip(:,snipIndex),2);
%     sniptime = getappdata(handles.figAL,'sniptime');
%     newlandmarkT = mean(sniptime(:,snipIndex));
%     current_sort_info = getappdata(handles.figAL,'current_sort_info');
%     current_sort_info.landmarkWaveform = ...
%         [current_sort_info.landmarkWaveform newlandmarkWaveform];
%     current_sort_info.landmarkT = ...
%         [current_sort_info.landmarkT newlandmarkT];
%     lmc = max(current_sort_info.landmarkClust,[],2) + 1;
%     current_sort_info.landmarkClust = ...
%         [current_sort_info.landmarkClust, lmc];
%     landmarkWaveformProj = ...
%       [current_sort_info.projectDirections'*current_sort_info.landmarkWaveform;
%       current_sort_info.landmarkT*current_sort_info.t2V];
%     setappdata(handles.figAL,'current_sort_info',current_sort_info);
%     setappdata(handles.figAL,'landmarkWaveformProj',landmarkWaveformProj)
%     update_landmarkIndex(handles)
%     add_landmarks_plot(handles)
%   end

  
function update_landmarkIndex(handles)
  current_sort_info = getappdata(handles.figAL,'current_sort_info');
  snipProj = getappdata(handles.figAL,'snipProj');
  landmarkWaveformProj = getappdata(handles.figAL,'landmarkWaveformProj');
  % Determine the closest landmark to each snippet
  [dist,landmarkIndex] = mindist(snipProj,landmarkWaveformProj);
  setappdata(handles.figAL,'landmarkIndex',landmarkIndex);


function add_landmarks_plot(handles)
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
  snipProj = getappdata(handles.figAL,'snipProj');
  landmarkWaveformProj = getappdata(handles.figAL,'landmarkWaveformProj');
  landmarkIndex = getappdata(handles.figAL,'landmarkIndex');
  clabelSnip = agglabel(landmarkIndex);
  current_sort_info = getappdata(handles.figAL,'current_sort_info');
  x_coord = get(handles.popupXVis,'Value');
  y_coord = get(handles.popupYVis,'Value');
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
      % Display the landmarks
      line(landmarkWaveformProj(x_coord,i),...
           landmarkWaveformProj(y_coord,i),...
           'Parent',handles.axesVisualizeShape,...
           'LineStyle','none',...
           'Marker','o',...
           'MarkerSize',12,...
           'MarkerFaceColor',col,...
           'MarkerEdgeColor',col,...
           'Visible',landmarkVisible,...
           'Tag','landmarkProjection');
    end
  end
  set(handles.axesVisualizeShape,'XTick',[],'YTick',[])
  % Put a '+' at 0,0 to facilitate identification of noise cloud
  line(0,0,'Parent',handles.axesVisualizeShape,'Marker','+','Color','k');
  %axis(handles.axesVisualizeShape,'equal')
  %if get(handles.boxLockAspectRatio,'Value')
  %  set(handles.axesVisualizeShape,'DataAspectRatio',[1 1 1]);
  %else
  %  set(handles.axesVisualizeShape,'DataAspectRatioMode','auto');
  %end
  % Retrieve/store axis limits (for zooming)
  %if ~isempty(xlim)
  %  set(handles.axesVisualizeShape,'XLim',xlim);
  %else
  %  setappdata(handles.axesVisualizeShape,'OXLim',...
  %    get(handles.axesVisualizeShape,'XLim'));
  %end
  %if ~isempty(ylim)
  %  set(handles.axesVisualizeShape,'YLim',ylim);
  %else
  %  setappdata(handles.axesVisualizeShape,'OYLim',...
  %    get(handles.axesVisualizeShape,'YLim'));
  %end


% --- Executes on selection change in popupYVis.
function popupYVis_Callback(hObject, eventdata, handles)
% hObject    handle to popupYVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupYVis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupYVis
  add_landmarks_plot(handles)

% --- Executes during object creation, after setting all properties.
function popupYVis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupYVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupXVis.
function popupXVis_Callback(hObject, eventdata, handles)
% hObject    handle to popupXVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupXVis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupXVis
  add_landmarks_plot(handles)


% --- Executes during object creation, after setting all properties.
function popupXVis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupXVis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cancelBtn.
function cancelBtn_Callback(hObject, eventdata, handles)
% hObject    handle to cancelBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figAL)


% --- Executes on button press in saveBtn.
function saveBtn_Callback(hObject, eventdata, handles)
% hObject    handle to saveBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curDir = cd;
saveDir = getappdata(handles.figAL,'saveDir');
cd(saveDir)
[filename,pathname] = uiputfile('*.mat','Pick an output file');
if ~isequal(filename,0)
  sort_info = getappdata(handles.figAL,'current_sort_info');
  save([pathname filename],'sort_info');
end
cd(curDir)


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


% --------------------------------------------------------------------
function File_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Open_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [filename,pathname,filterindex] = uigetfile('*.mat',['Choose sorting ' ...
                      'file']);
  if ~isequal(filename,0)
    add_landmarks_open(handles,pathname,filename);
  end
  

