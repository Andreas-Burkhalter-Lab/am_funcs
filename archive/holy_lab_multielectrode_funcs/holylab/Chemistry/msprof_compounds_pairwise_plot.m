function varargout = msprof_compounds_pairwise_plot(varargin)
% MSPROF_COMPOUNDS_PAIRWISE_PLOT: compare concentrations of compounds between samples
% Syntax:
%   msprof_compounds_pairwise_plot(s)
% where
%   s is a structure of the form produced by msprof_concentrations.
%
% See also: msprof_concentrations.

% Copyright 2011 by Timothy E. Holy

% MSPROF_COMPOUNDS_PAIRWISE_PLOT M-file for msprof_compounds_pairwise_plot.fig
%      MSPROF_COMPOUNDS_PAIRWISE_PLOT, by itself, creates a new MSPROF_COMPOUNDS_PAIRWISE_PLOT or raises the existing
%      singleton*.
%
%      H = MSPROF_COMPOUNDS_PAIRWISE_PLOT returns the handle to a new MSPROF_COMPOUNDS_PAIRWISE_PLOT or the handle to
%      the existing singleton*.
%
%      MSPROF_COMPOUNDS_PAIRWISE_PLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSPROF_COMPOUNDS_PAIRWISE_PLOT.M with the given input arguments.
%
%      MSPROF_COMPOUNDS_PAIRWISE_PLOT('Property','Value',...) creates a new MSPROF_COMPOUNDS_PAIRWISE_PLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msprof_compounds_pairwise_plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msprof_compounds_pairwise_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msprof_compounds_pairwise_plot

% Last Modified by GUIDE v2.5 26-Aug-2010 15:47:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msprof_compounds_pairwise_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @msprof_compounds_pairwise_plot_OutputFcn, ...
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


% --- Executes just before msprof_compounds_pairwise_plot is made visible.
function msprof_compounds_pairwise_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to msprof_compounds_pairwise_plot (see VARARGIN)

% Choose default command line output for msprof_compounds_pairwise_plot
handles.output = hObject;

s = varargin{1};
% Massage the m/z, trange fields into an array suitable for mdexplore
n_compounds = length(s.mz);
mzt = repmat(struct('mz',[],'span',[]),1,n_compounds);
for i = 1:n_compounds
  mzt(i).mz = s.mz(i);
  mzt(i).span = s.xrange(:,i);
end
s.mzt = mzt;
setappdata(handles.figMain,'compounds',s);
set(handles.popX,'String',s.stimtag,'Value',1);
set(handles.popY,'String',s.stimtag,'Value',2);

% Update handles structure
guidata(hObject, handles);

% Draw
mspcpp_plot(handles)

% UIWAIT makes msprof_compounds_pairwise_plot wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = msprof_compounds_pairwise_plot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popY.
function popY_Callback(hObject, eventdata, handles)
% hObject    handle to popY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popY contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popY
mspcpp_plot(handles)

% --- Executes during object creation, after setting all properties.
function popY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popX.
function popX_Callback(hObject, eventdata, handles)
% hObject    handle to popX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popX
mspcpp_plot(handles)

% --- Executes during object creation, after setting all properties.
function popX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnExport.
function btnExport_Callback(hObject, eventdata, handles)
% hObject    handle to btnExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkboxLog.
function checkboxLog_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLog
if get(hObject,'Value')
  set(handles.axesPair,'XScale','log','YScale','log')
else  
  set(handles.axesPair,'XScale','linear','YScale','linear')
end

function mspcpp_plot(handles)
  logFlag = get(handles.checkboxLog,'Value');
  s = getappdata(handles.figMain,'compounds');
  xIndex = get(handles.popX,'Value');
  yIndex = get(handles.popY,'Value');
  ops = struct('hax',handles.axesPair,...
    'plotfunc',@(pk) msprof_choose_peak_timecourses_gui(pk.mz,pk.span));
  cla(handles.axesPair)
  mdexplore(s.I([xIndex yIndex],:),s.mzt,ops);
  if logFlag
    set(handles.axesPair,'XScale','log','YScale','log')
  end
  
    