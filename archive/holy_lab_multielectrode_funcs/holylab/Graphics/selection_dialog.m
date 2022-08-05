function varargout = selection_dialog(varargin)
% SELECTION_DIALOG M-file for selection_dialog.fig
%      SELECTION_DIALOG, by itself, creates a new SELECTION_DIALOG or raises the existing
%      singleton*.
%
%      H = SELECTION_DIALOG returns the handle to a new SELECTION_DIALOG or the handle to
%      the existing singleton*.
%
%      SELECTION_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTION_DIALOG.M with the given input arguments.
%
%      SELECTION_DIALOG('Property','Value',...) creates a new SELECTION_DIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selection_dialog_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selection_dialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selection_dialog

% Last Modified by GUIDE v2.5 29-Jun-2006 12:05:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selection_dialog_OpeningFcn, ...
                   'gui_OutputFcn',  @selection_dialog_OutputFcn, ...
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


% --- Executes just before selection_dialog is made visible.
function selection_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selection_dialog (see VARARGIN)

% Choose default command line output for selection_dialog
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if(length(varargin)==0)
   error('usage: indices=selection_dialog(strings, colors)');
end
if(length(varargin)==1)
   colors={};
else
   colors=varargin{2};
end
strings=varargin{1};
% create checkboxes
fig=hObject;
cbGeometry=get(handles.cbFirst, 'position');
height=cbGeometry(4);
spacing=0.1;
geometry=get(fig, 'position');
set(fig, 'position', [geometry(1) geometry(2) geometry(3) geometry(4)+height*(1+spacing)*(length(strings)-1)]);
hCheckboxes=handles.cbFirst;
for idx=2:length(strings)
   curCheckbox=copyobj(handles.cbFirst, fig);
   set(curCheckbox, 'position', [cbGeometry(1) cbGeometry(2)+(idx-1)*height*(1+spacing) cbGeometry(3:4)]);
   hCheckboxes(end+1)=curCheckbox;
end
for idx=1:length(strings)
   curCheckbox=hCheckboxes(end-idx+1);
   set(curCheckbox, 'string', strings{idx});
   set(curCheckbox, 'userdata', idx);
   if(~isempty(colors))
      set(curCheckbox, 'foregroundColor', colors{idx});
   end
end

% UIWAIT makes selection_dialog wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = selection_dialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
hCheckboxes=findobj(hObject, 'style', 'checkbox', 'value', 1);
selectedIndices=get(hCheckboxes, 'userdata');
if(iscell(selectedIndices))
   selectedIndices=[selectedIndices{:}];
end
varargout{1} = selectedIndices;
close(hObject);

% --- Executes on button press in btnOk.
function btnOk_Callback(hObject, eventdata, handles)
% hObject    handle to btnOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(gcf);

% --- Executes on button press in cbFirst.
function cbFirst_Callback(hObject, eventdata, handles)
% hObject    handle to cbFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbFirst


