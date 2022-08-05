function varargout = imselrectgui(varargin)
% IMSELRECTGUI M-file for imselrectgui.fig
% Syntax:
%   [rangeout,csout] = imselrectgui(imagedata,rangein,csin)
%      IMSELRECTGUI, by itself, creates a new IMSELRECTGUI or raises the existing
%      singleton*.
%
%      H = IMSELRECTGUI returns the handle to a new IMSELRECTGUI or the handle to
%      the existing singleton*.
%
%      IMSELRECTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMSELRECTGUI.M with the given input arguments.
%
%      IMSELRECTGUI('Property','Value',...) creates a new IMSELRECTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imselrectgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imselrectgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imselrectgui

% Last Modified by GUIDE v2.5 10-Feb-2004 21:26:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imselrectgui_OpeningFcn, ...
                   'gui_OutputFcn',  @imselrectgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imselrectgui is made visible.
function imselrectgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
handles.output = hObject;
handles.imagedata = varargin{1};

axes(handles.imaxes);
handles.himage = image(handles.imagedata,'CDataMapping','scaled');
set(gca,'Visible','off');
colormap(gray)
xlim=get(gca, 'xlim');
ylim = get(gca,'YLim');
xdist=diff(xlim); ydist=diff(ylim);
handles.line_x1=line(xlim(1)+[xdist/4 xdist/4], ylim);
handles.line_x2=line(xlim(1)+[xdist*3/4 xdist*3/4], ylim);
handles.line_y1=line(xlim, ylim(1)+[ydist/4 ydist/4]);
handles.line_y2=line(xlim, ylim(1)+[ydist*3/4 ydist*3/4]);
set([handles.line_x1 handles.line_x2 handles.line_y1 handles.line_y2], 'LineWidth', 3);

% Update handles structure
guidata(hObject, handles);

drag_line(handles.line_x1);
drag_line(handles.line_x2);
drag_line(handles.line_y1);
drag_line(handles.line_y2);




% --- Outputs from this function are returned to the command line.
function varargout = imselrectgui_OutputFcn(hObject, eventdata, handles)
% UIWAIT makes imselrectgui wait for user response (see UIRESUME)
uiwait(handles.figure1);
if ~ishandle(hObject)
  % If the user closed the figure, it's like hitting cancel
  varargout{1} = [];
else
  handles = guidata(hObject);
  varargout{1} = handles.output;
end



% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% handles.cs = [];
% guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in OKbutton.
function OKbutton_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);



