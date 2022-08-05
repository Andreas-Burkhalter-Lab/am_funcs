function varargout = drawmask_gui(varargin)
% DRAWMASK_GUI: specify the "interesting" region interactively.
% Syntax:
%   mask = drawmask_gui(stk)
% where
%   stk is a nx-by-ny-by-nz stack of images
% and
%   mask is a nx-by-ny-by-nz logical array, true in the pixels inside the
%     polygon that will be defined by the user.
  
  
% DRAWMASK_GUI M-file for drawmask_gui.fig
%      DRAWMASK_GUI, by itself, creates a new DRAWMASK_GUI or raises the existing
%      singleton*.
%
%      H = DRAWMASK_GUI returns the handle to a new DRAWMASK_GUI or the handle to
%      the existing singleton*.
%
%      DRAWMASK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRAWMASK_GUI.M with the given input arguments.
%
%      DRAWMASK_GUI('Property','Value',...) creates a new DRAWMASK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drawmask_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drawmask_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help drawmask_gui

% Last Modified by GUIDE v2.5 18-May-2010 14:21:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drawmask_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @drawmask_gui_OutputFcn, ...
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


% --- Executes just before drawmask_gui is made visible.
function drawmask_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drawmask_gui (see VARARGIN)

% Choose default command line output for drawmask_gui
handles.output = hObject;

% Draw the up & down buttons
icons = swicons;
set(handles.pushbuttonUp,'CData',icons.up);
set(handles.pushbuttonDown,'CData',icons.down);

% Parse inputs
stk = varargin{1};
nz = size(stk,3);
polyc = cell(1,nz);
setappdata(handles.figureMain,'stk',stk);
setappdata(handles.figureMain,'polyc',polyc);
setappdata(handles.figureMain,'nz',nz);
set(handles.textTotFrameNumber,'String',sprintf('/%d',nz));

% Update handles structure
guidata(hObject, handles);

update_display(handles)

% UIWAIT makes drawmask_gui wait for user response (see UIRESUME)
uiwait(handles.figureMain);


% --- Outputs from this function are returned to the command line.
function varargout = drawmask_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
stk = getappdata(handles.figureMain,'stk');
polyc = getappdata(handles.figureMain,'polyc');
mask = false(size(stk));
[nx ny nz] = size(stk);
x = 1:nx;
y = 1:ny;
[X,Y] = ndgrid(x,y);
for i = 1:nz
  thispoly = getpoly(polyc,i);
  mask(:,:,i) = inpolygon(X,Y,thispoly(:,2),thispoly(:,1));
end
varargout = {mask};
close(handles.figureMain)


function editFrameNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameNumber as text
%        str2double(get(hObject,'String')) returns contents of editFrameNumber as a double
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  nz = getappdata(handles.figureMain,'nz');
  framenumber = round(framenumber);
  if (framenumber < 1)
    framenumber = 1;
  elseif (framenumber > nz)
    framenumber = nz;
  end
  set(handles.editFrameNumber,'String',num2str(framenumber));
  update_display(handles)
  

% --- Executes during object creation, after setting all properties.
function editFrameNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonUp.
function pushbuttonUp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  framenumber = framenumber-1;
  if (framenumber < 1)
    framenumber = 1;
  end
  set(handles.editFrameNumber,'String',num2str(framenumber));
  update_display(handles)

% --- Executes on button press in pushbuttonDown.
function pushbuttonDown_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  nz = getappdata(handles.figureMain,'nz');
  framenumber = framenumber+1;
  if (framenumber > nz)
    framenumber = nz;
  end
  set(handles.editFrameNumber,'String',num2str(framenumber));
  update_display(handles)


% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  uiresume(handles.figureMain)

function update_display(handles)
  stk = getappdata(handles.figureMain,'stk');
  polyc = getappdata(handles.figureMain,'polyc');
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  himg = imshowsc(stk(:,:,framenumber),'Parent',handles.axesImage);
  set(himg,'HitTest','off');
  thispoly = getpoly(polyc,framenumber);
  if ~isempty(thispoly)
    thispoly = [thispoly; thispoly(1,:)];
    hline = line(thispoly(:,1),thispoly(:,2),'Parent',handles.axesImage,'Marker','o');
    set(hline,'Tag','poly');
    set(hline,'ButtonDownFcn',@start_dragging_vertex);
  else
    set(handles.axesImage,'ButtonDownFcn',@draw_poly);
  end
  update_buttons(handles)

function update_buttons(handles)
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  nz = getappdata(handles.figureMain,'nz');
  if (framenumber == 1)
    set(handles.pushbuttonUp,'Enable','off')
  else
    set(handles.pushbuttonUp,'Enable','on')
  end
  if (framenumber == nz)
    set(handles.pushbuttonDown,'Enable','off')
  else
    set(handles.pushbuttonDown,'Enable','on')
  end


function thispoly = getpoly(polyc,framenumber)
  thispoly = [];
  % bracket the framenumber in terms of "known" polygons
  emptyFlag = cellfun(@isempty,polyc);
  definedIndex = find(~emptyFlag);
  index = find(definedIndex < framenumber,1,'last');
  if isempty(index)
    % There are none equal or bigger; see if there is one smaller
    if isempty(definedIndex)
      return
    end
    index = 1;
  end
  if (index < length(definedIndex))
    % we have a full bracket, use linear interpolation
    index = index + [0 1];
    index = definedIndex(index);
    f = (framenumber - index(1))/diff(index);
    thispoly = (1-f)*polyc{index(1)} + f*polyc{index(2)};
  else
    % we do not have a full bracket, just copy the result
    thispoly = polyc{definedIndex(index)};
  end

function start_dragging_vertex(src,eventdata)
  set(gcbf,'WindowButtonMotionFcn',@drag_vertex);
  set(gcbf,'WindowButtonUpFcn',@release_vertex);
  
function handled = drag_vertex(hfig,event_args)
  cp = get(gca,'CurrentPoint');
  cp = cp(1,1:2);
  sender = findobj(hfig,'Tag','poly');
  % Find the closest vertex
  x = get(sender,'XData');
  y = get(sender,'YData');
  d = sqrt((x-cp(1)).^2 + (y-cp(2)).^2);
  [mind,index] = min(d);
  if (index == 1)
    index = [1 length(d)];
  end
  % Set the closest vertex to the current point
  x(index) = cp(1);
  y(index) = cp(2);
  set(sender,'XData',x,'YData',y);
  handled = true;
  
function handled = release_vertex(sender,event_args)
  hfig = gcbf;
  set(hfig,'WindowButtonMotionFcn',[]);
  set(hfig,'WindowButtonUpFcn',[]);
  handles = guidata(hfig);
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  sender = findobj(hfig,'Tag','poly');
  x = get(sender,'XData');
  y = get(sender,'YData');
  setpoly(handles,x,y,framenumber)
  handled = true;
  
function draw_poly(obj,eventdata)
  handles = guidata(obj);
  set(handles.axesImage,'ButtonDownFcn',[]);
  [px,py] = GetSelPolygon('go','r');
  framenumber = str2num(get(handles.editFrameNumber,'String'));
  setpoly(handles,px,py,framenumber);
  
function setpoly(handles,px,py,framenumber)
  px = px(1:end-1);
  py = py(1:end-1);
  polyc = getappdata(handles.figureMain,'polyc');
  polyc{framenumber} = [px(:) py(:)];
  setappdata(handles.figureMain,'polyc',polyc);
  update_display(handles)
  
 function icons = swicons
  iconsize = [30 30];
  bkgrnd = get(0,'DefaultUicontrolBackgroundColor');
  blank = zeros(iconsize) + bkgrnd(1);
  
  right0 = blank;
  right1 = blank;
  index = arrowmask(iconsize,0.8,1);
  right1(index) = 1;
  right0(index) = 0;
  % Long green arrow pointing right
  right(:,:,1) = right0;
  right(:,:,2) = right1;
  right(:,:,3) = right0;

  rightsm = blank;
  index = arrowmask(iconsize,0.5,1);
  rightsm(index) = 0;
  rightsm = repmat(rightsm,[1 1 3]);  % Short black arrow pointing right

  left = right(:,end:-1:1,:);
  leftsm = rightsm(:,end:-1:1,:);
  
  icons.right = right;
  icons.rightsm = rightsm;
  icons.left = left;
  icons.leftsm = leftsm;
  icons.up = permute(left,[2 1 3]);
  icons.upsm = permute(leftsm,[2 1 3]);
  icons.down = permute(right,[2 1 3]);
  icons.downsm = permute(rightsm,[2 1 3]);
  icons.stop = repmat(0*blank,[1 1 3]);

function index = arrowmask(iconsize,width,direction)
  [row,col] = find(ones(iconsize));
  % Make them go from 0 to 1
  row = (row-1)/(iconsize(1)-1);
  col = (col-1)/(iconsize(2)-1);
  if (direction > 0)
    index = find(2*abs(row-0.5) <= (width-col)/width);
  else
    index = find(2*abs(row-0.5) + (1 - width)/width <= col/width);
  end
  return;