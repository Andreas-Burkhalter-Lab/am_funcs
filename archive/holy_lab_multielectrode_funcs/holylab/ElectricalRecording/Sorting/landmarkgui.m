function varargout = landmarkgui(varargin)
% LANDMARKGUI M-file for landmarkgui.fig
%      LANDMARKGUI, by itself, creates a new LANDMARKGUI or raises the existing
%      singleton*.
%
%      H = LANDMARKGUI returns the handle to a new LANDMARKGUI or the handle to
%      the existing singleton*.
%
%      LANDMARKGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LANDMARKGUI.M with the given input arguments.
%
%      LANDMARKGUI('Property','Value',...) creates a new LANDMARKGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before landmarkgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to landmarkgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help landmarkgui

% Last Modified by GUIDE v2.5 21-Feb-2006 00:34:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @landmarkgui_OpeningFcn, ...
                   'gui_OutputFcn',  @landmarkgui_OutputFcn, ...
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


% --- Executes just before landmarkgui is made visible.
function landmarkgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to landmarkgui (see VARARGIN)

% Choose default command line output for landmarkgui
handles.output = hObject;

landmarkPosition = varargin{1};
snipall = varargin{2};
landmarkIndex = varargin{3};
t2V = varargin{4};

clabel = agglabel(landmarkIndex);
for i = 1:length(clabel)
  col = unique_color(i,length(clabel));
  % Draw the snippets
  line(snipall(1,clabel{i}),snipall(2,clabel{i}),...
       'Parent',handles.axesVisualizeShape,...
       'LineStyle','none',...
       'Marker','.',...
       'MarkerSize',1,...
       'Color',col);
  % Draw the landmark
  line(landmarkPosition(1,i),landmarkPosition(2,i),...
       'Parent',handles.axesVisualizeShape,...
       'LineStyle','none',...
       'Marker','o',...
       'MarkerSize',12,...
       'MarkerFaceColor',col,...
       'MarkerEdgeColor',col);
end

set(handles.edit_t2V,'String',num2str(t2V));


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes landmarkgui wait for user response (see UIRESUME)
uiwait(handles.figLM);


% --- Outputs from this function are returned to the command line.
function varargout = landmarkgui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
varargout{1} = 1-getappdata(handles.figLM,'approve');
varargout{2} = str2num(get(handles.edit_t2V,'String'));
close(handles.figLM)

function edit_t2V_Callback(hObject, eventdata, handles)
% hObject    handle to edit_t2V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_t2V as text
%        str2double(get(hObject,'String')) returns contents of edit_t2V as a double


% --- Executes during object creation, after setting all properties.
function edit_t2V_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_t2V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnRetry.
function btnRetry_Callback(hObject, eventdata, handles)
% hObject    handle to btnRetry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figLM,'approve',0);
uiresume(handles.figLM)

% --- Executes on button press in btnApprove.
function btnApprove_Callback(hObject, eventdata, handles)
% hObject    handle to btnApprove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figLM,'approve',1);
uiresume(handles.figLM)


