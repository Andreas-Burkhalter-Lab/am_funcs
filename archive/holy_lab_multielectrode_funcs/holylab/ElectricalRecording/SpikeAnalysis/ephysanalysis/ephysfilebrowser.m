function varargout = ephysfilebrowser(varargin)
% EPHYSFILEBROWSER: GUI for browsing files and channels
% Syntax:
%   ephysfilebrowser(ephys)
%   ephysfilebrowser(ephys,options)
% where
%   ephys is a structure array, typically containing ephys data for
%     different files;
%   options (optional) is passed to EPHYSCHANNELBROWSER to affect its
%     display.
%
% Note: it's assumed that the same channels are recorded in each file.
%
% See also: EPHYSCHANNELBROWSER.

% Copyright 2005 by Timothy E. Holy

% EPHYSFILEBROWSER M-file for ephysfilebrowser.fig
%      EPHYSFILEBROWSER, by itself, creates a new EPHYSFILEBROWSER or raises the existing
%      singleton*.
%
%      H = EPHYSFILEBROWSER returns the handle to a new EPHYSFILEBROWSER or the handle to
%      the existing singleton*.
%
%      EPHYSFILEBROWSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EPHYSFILEBROWSER.M with the given input arguments.
%
%      EPHYSFILEBROWSER('Property','Value',...) creates a new EPHYSFILEBROWSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ephysfilebrowser_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ephysfilebrowser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help ephysfilebrowser

% Last Modified by GUIDE v2.5 01-May-2005 14:45:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ephysfilebrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @ephysfilebrowser_OutputFcn, ...
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


% --- Executes just before ephysfilebrowser is made visible.
function ephysfilebrowser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ephysfilebrowser (see VARARGIN)

% Choose default command line output for ephysfilebrowser
handles.output = hObject;
handles.ephys = varargin{1};
handles.chanfigure = [];

set(handles.listbox1,'String',{handles.ephys.basefilename});
set(handles.listbox2,'String',num2str(handles.ephys(1).channels'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ephysfilebrowser wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ephysfilebrowser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        listbox1
  handles = epfb_displaychannel(handles);
  guidata(hObject, handles);

  
% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set(hObject,'String',{handles.ephys.basefilename});

% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2
  handles = epfb_displaychannel(handles);
  guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%set(hObject,'String',num2str(handles.ephys(1).channels));


function handles = epfb_displaychannel(handles)
  selected_file = get(handles.listbox1,'Value');
  selected_channel = get(handles.listbox2,'Value');
  close(handles.chanfigure(ishandle(handles.chanfigure)));
  etmp = ephyssubchan(handles.ephys(selected_file),...
                      handles.ephys(selected_file).channels(selected_channel));
  hfig = ephyschannelbrowser(etmp);
  handles.chanfigure = hfig;