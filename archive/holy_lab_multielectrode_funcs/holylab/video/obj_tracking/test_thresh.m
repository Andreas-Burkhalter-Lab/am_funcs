function varargout = test_thresh(varargin)
% TEST_THRESH GUI to set the contrast for an image
% Syntax:
%   [rangeout,csout] = test_thresh(imagedata,rangein,csin)
% where
%   imagedata is the image to be displayed/histogrammed
%   rangein is a 2-vector containing the current setting of the
%     colorlimits, i.e., the min/max values which go to black/white (or
%     bottom/top of the colormap)
%   csin is the color saturated pixels option - set to 0 for no, 1 for yes
% and
%   rangeout is a 2-vector giving the user-chosen colorlimits.
%   csout?

%      TEST_THRESH, by itself, creates a new TEST_THRESH or raises the existing
%      singleton*.
%
%      H = TEST_THRESH returns the handle to a new TEST_THRESH or the handle to
%      the existing singleton*.
%
%      TEST_THRESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_THRESH.M with the given input arguments.
%
%      TEST_THRESH('Property','Value',...) creates a new TEST_THRESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imrangegui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_thresh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_thresh

% Last Modified by GUIDE v2.5 02-Mar-2009 10:18:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_thresh_OpeningFcn, ...
                   'gui_OutputFcn',  @test_thresh_OutputFcn, ...
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


% --- Executes just before test_thresh is made visible.
function test_thresh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
handles.output = hObject;
handles.imagedata = double(varargin{1});
handles.range = double(varargin{2});
set(handles.colorsaturated,'Value',varargin{3});
% Set up the thumbnail image
axes(handles.imaxes);
handles.himage = image(handles.imagedata,'CDataMapping','scaled');
set(gca,'CLim',handles.range,'Visible','off');
colormap(gray)
axis image;
% Compute the histogram
[handles.nperbin,handles.xbin] = hist(handles.imagedata(:),100);
% Create the plot objects
axes(handles.histaxes);
hhist = bar(handles.xbin,handles.nperbin);
set(hhist,'HitTest','off'); % Keep histogram from intercepting mouse clicks
ylim = get(gca,'YLim');
ylim(1) = max(1,ylim(1));
handles.rangeline(1) = line(handles.range([1 1]),ylim); % Lower threshold line
% handles.rangeline(2) = line(handles.range([2 2]),ylim); % Upper threshold line
options.move = @imrangegui_moveline; % Set up callbacks for moving threshold lines
options.done = @imrangegui_doneline;
set(handles.rangeline,'EraseMode','xor','ButtonDownFcn',{@moveline,options})
%line([0 0],ylim,'Color',[0 0 1]); % Blue line at minimum
%line([4095 4095],ylim,'Color',[1 0 0]); % Red line at maximum
% Calculate the x-axes and set them
set(handles.histaxes,'XLim',imrangegui_xlim(handles));
set(handles.histaxes,'YScale','log','YLim',ylim);
% Fill in the edit text
set(handles.editmin,'String',num2str(handles.range(1)));
set(handles.editmax,'String',num2str(handles.range(2)));
% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = test_thresh_OutputFcn(hObject, eventdata, handles)
varargout={[], []};
return;
% UIWAIT makes test_thresh wait for user response (see UIRESUME)
uiwait(handles.figure1);
if ~ishandle(hObject)
  % If the user closed the figure, it's like hitting cancel
  varargout{1} = [];
  varargout{2} = [];
else
  handles = guidata(hObject);
  varargout{1} = handles.range;
  varargout{2} = get(handles.colorsaturated,'Value');
  close(handles.figure1);
end

% --- Executes on button press in refreshxlim.
function refreshxlim_Callback(hObject, eventdata, handles)
set(handles.histaxes,'XLim',imrangegui_xlim(handles));

% editmin
% --- Executes during object creation, after setting all properties.
function editmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editmin_Callback(hObject, eventdata, handles)
val = round(str2num(get(hObject,'String')));
if ~isempty(val)
  handles.range(1) = val;
  handles.range(2)=val+0.1;
  % Figure out which of the range lines is actually the minimum
  xd = get(handles.rangeline,'XData');
  xdata(1) = xd(1); 
  xdata(2) = xdata(1)+0.1;
  % [m,minindx] = min(xdata);
  % Set the position
  set(handles.rangeline,'XData',[val val]);
end
guidata(hObject,handles);
set(handles.imaxes,'CLim',handles.range);

% editmax
% --- Executes during object creation, after setting all properties.
function editmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function editmax_Callback(hObject, eventdata, handles)
val = str2num(get(hObject,'String'));
if ~isempty(val)
  handles.range(2) = val;
  % Figure out which of the range lines is actually the maximum
  xd = get(handles.rangeline,'XData');
  xdata(1) = xd{1}(1); xdata(2) = xd{2}(1);
  [m,maxindx] = max(xdata);
  % Set the position
  set(handles.rangeline(maxindx),'XData',[val val]);
end
guidata(hObject,handles);
set(handles.imaxes,'CLim',handles.range);

% --- Executes on button press in colorsaturated.
function colorsaturated_Callback(hObject, eventdata, handles)
% No response needed


% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
handles.range = [];
handles.cs = [];
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in OKbutton.
function OKbutton_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Callbacks for dragging the range lines
function imrangegui_moveline(hObject, eventdata)
handles = guidata(hObject);
% While dragging, just update the edit text
xd = get(handles.rangeline,'XData');
range(1) = xd(1); 
range(2) = range(1)+0.1;
% range = sort(range);  % In case the user crossed the two lines
handles.range = range;
set(handles.editmin,'String',num2str(range(1)));
set(handles.editmax,'String',num2str(range(2)));
guidata(hObject,handles);

function imrangegui_doneline(hObject,eventdata)
% When the user releases the mouse button, update the edit text, but also
% change the contrast on the image thumbnail
imrangegui_moveline(hObject,eventdata);
handles = guidata(hObject);
set(handles.imaxes,'CLim',handles.range);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xlim = imrangegui_xlim(handles)
xlim(1) = min([handles.xbin(1) handles.range(1)]);
xlim(2) = max([handles.xbin(end) handles.range(2)]);
xlim = xlim + [-1 1]*0.1*diff(xlim); % Expand by 10% in each direction
