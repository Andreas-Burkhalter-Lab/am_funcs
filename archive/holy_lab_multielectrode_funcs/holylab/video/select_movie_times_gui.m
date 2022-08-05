function varargout = select_movie_times_gui(varargin)
% SELECT_MOVIE_TIMES_GUI: let user mark particular times in a movie
%
% This shows a movie using mplay; a user navigates to a desired frame
% and then clicks "Set."  Multiple marking points can be requested as a
% group.
%
% Syntax:
%   frametimes = select_movie_times_gui(mov,n_times)
% where
%   mov is a movie object (see mmreader or mmreader_ffmpeg)
%   n_times is the number of times you want marked in the movie
% and
%   frametimes is a vector containing the time (in seconds) for each "mark"
%     set by the user.
%
% Alternatively, you can do this:
%   frametimes = select_movie_times_gui(mov,2,'start','stop')
% in which case the different times are "named" for the user.
%
% See also: mplay, mmreader, mmreader_ffmpeg.
  
% Copyright 2010-2011 by Timothy E. Holy
    
  
%      SELECT_MOVIE_TIMES_GUI, by itself, creates a new SELECT_MOVIE_TIMES_GUI or raises the existing
%      singleton*.
%
%      H = SELECT_MOVIE_TIMES_GUI returns the handle to a new SELECT_MOVIE_TIMES_GUI or the handle to
%      the existing singleton*.
%
%      SELECT_MOVIE_TIMES_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_MOVIE_TIMES_GUI.M with the given input arguments.
%
%      SELECT_MOVIE_TIMES_GUI('Property','Value',...) creates a new SELECT_MOVIE_TIMES_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_movie_times_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_movie_times_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_movie_times_gui

% Last Modified by GUIDE v2.5 20-Aug-2010 15:20:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_movie_times_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @select_movie_times_gui_OutputFcn, ...
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


% --- Executes just before select_movie_times_gui is made visible.
function select_movie_times_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_movie_times_gui (see VARARGIN)

% Choose default command line output for select_movie_times_gui
handles.output = [];

mov = varargin{1};
n_times = Inf;
if (length(varargin) > 1)
  n_times = varargin{2};
  if length(varargin) > 2
    user_msg = varargin(3:end);
    setappdata(hObject,'user_msg',user_msg);
    set(handles.textMessage,'String',['Select ', user_msg{1}, ' time and click "Set"']);
  end
end
setappdata(hObject,'n_times',n_times);

% Set up the movie player
controls = mplay(mov);
setappdata(hObject,'controls',controls);
setappdata(hObject,'movie',mov);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_movie_times_gui wait for user response (see UIRESUME)
uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = select_movie_times_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
controls = getappdata(hObject,'controls');
s = struct(controls);
close(s.hfig);
close(handles.figMain);


% --- Executes on button press in btnSet.
function btnSet_Callback(hObject, eventdata, handles)
% hObject    handle to btnSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hfig = handles.figMain;
mov = getappdata(hfig,'movie');
t = mov.CurrentTime;
handles.output = [handles.output t];
guidata(handles.figMain,handles)
n_next = length(handles.output)+1;
n_times = getappdata(hfig,'n_times');
if (n_next > n_times)
  uiresume
end
if isappdata(hfig,'user_msg')
  user_msg = getappdata(hfig,'user_msg');
  if (length(user_msg) >= n_next)
    set(handles.textMessage,'String',['Select ', user_msg{n_next}, ' time and click "Set"']);
  end
end

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output = [];
guidata(handles.figMain,handles)
uiresume

% --- Executes on button press in btnDone.
function btnDone_Callback(hObject, eventdata, handles)
% hObject    handle to btnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume
