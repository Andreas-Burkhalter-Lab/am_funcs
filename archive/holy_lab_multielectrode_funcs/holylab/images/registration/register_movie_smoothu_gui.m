function varargout = register_movie_smoothu_gui(varargin)
% Syntax:
%   [udatanew,errnew] = register_movie_smoothu_gui(udata,err,errfunc,pcoptions,pixel_spacing)
  
% REGISTER_MOVIE_SMOOTHU_GUI M-file for register_movie_smoothu_gui.fig
%      REGISTER_MOVIE_SMOOTHU_GUI, by itself, creates a new REGISTER_MOVIE_SMOOTHU_GUI or raises the existing
%      singleton*.
%
%      H = REGISTER_MOVIE_SMOOTHU_GUI returns the handle to a new REGISTER_MOVIE_SMOOTHU_GUI or the handle to
%      the existing singleton*.
%
%      REGISTER_MOVIE_SMOOTHU_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_MOVIE_SMOOTHU_GUI.M with the given input arguments.
%
%      REGISTER_MOVIE_SMOOTHU_GUI('Property','Value',...) creates a new REGISTER_MOVIE_SMOOTHU_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before register_movie_smoothu_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_movie_smoothu_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_movie_smoothu_gui

% Last Modified by GUIDE v2.5 23-Jan-2011 13:33:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @register_movie_smoothu_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @register_movie_smoothu_gui_OutputFcn, ...
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


% --- Executes just before register_movie_smoothu_gui is made visible.
function register_movie_smoothu_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_movie_smoothu_gui (see VARARGIN)

  udata = varargin{1};
  err = varargin{2};
  errfunc = varargin{3};
  pcoptions = varargin{4};
  pixel_spacing = varargin{5};
  
  funch = register_gui_utilities;

  % Find stacks that have u defined
  haveu = cellfun(@(x) ~isempty(x),udata);
  stacknums = find(haveu);
  % Determine the largest spatial grid size
  uszmax = [0 0 0];
  for i = 1:length(stacknums)
    thissz = size(udata{stacknums(i)}{1});
    thissz(end+1:3) = 1;
    uszmax = max(uszmax,thissz);
  end
  % Interpolate all up to the maximum size
  sz = cat(1,pcoptions.pyramid.sz);
  desired_level = find(all(sz == repmat(uszmax,size(sz,1),1),2));
  for i = 1:length(stacknums)
    thisstacknum = stacknums(i);
    udata{thisstacknum} = funch.match_u(udata{thisstacknum},desired_level,pcoptions); %#ok<FNDSB>
  end
  
  setappdata(handles.figMain,'udata0',udata);
  setappdata(handles.figMain,'err0',err);
  setappdata(handles.figMain,'errfunc',errfunc);
  setappdata(handles.figMain,'stacknums',stacknums);
  setappdata(handles.figMain,'pcoptions',pcoptions);
  setappdata(handles.figMain,'pixel_spacing',pixel_spacing);
  setappdata(handles.figMain,'funch',funch);
  setappdata(handles.figMain,'usz',uszmax);
  
  % Prepare graphs
  set([handles.axesUDiffOriginal handles.axesUDiffNew handles.axesErr],'XLim',[0 length(udata)+1],'TickDir','out');
  set([handles.axesUDiffOriginal handles.axesUDiffNew],'XTick',[])
  funch.plot_udiff(handles.axesUDiffOriginal,udata,pcoptions,pixel_spacing)
  line(1:length(err),err,'LineStyle','none','Marker','x','Parent',handles.axesErr);
  
  % Choose default command line output for register_movie_smoothu_gui
  handles.output = hObject;
  
  % Update handles structure
  guidata(hObject, handles);
  
  % UIWAIT makes register_movie_smoothu_gui wait for user response (see UIRESUME)
  uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = register_movie_smoothu_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
  varargout = handles.output;
  delete(handles.figMain);
else
  varargout = {[],[]};
end


% --- Executes on button press in btnFilter.
function btnFilter_Callback(hObject, eventdata, handles)
% hObject    handle to btnFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  filter_size = str2double(get(handles.editFilterSize,'String'));
  udata = getappdata(handles.figMain,'udata0');
  n_stacks = length(udata);
  usz = getappdata(handles.figMain,'usz');
  haveu = cellfun(@(x) ~isempty(x),udata);
  udata = udata(haveu);
  for dimIndex = 1:3
    ucoord = cellfun(@(c) c{dimIndex}(:),udata,'UniformOutput',false);
    u = cat(2,ucoord{:});
    uf = u;
    for i = 1:size(u,1)
      uf(i,:) = medfilt(u(i,:),filter_size);
    end
    for stackIndex = 1:length(udata)
      udata{stackIndex}{dimIndex} = reshape(uf(:,stackIndex),usz);
    end
  end
  udatanew = cell(1,n_stacks);
  udatanew(haveu) = udata;
  setappdata(handles.figMain,'udata',udatanew);
  cla(handles.axesUDiffNew);
  pcoptions = getappdata(handles.figMain,'pcoptions');
  pixel_spacing = getappdata(handles.figMain,'pixel_spacing');
  funch = getappdata(handles.figMain,'funch');
  funch.plot_udiff(handles.axesUDiffNew,udatanew,pcoptions,pixel_spacing)


function editFilterSize_Callback(hObject, eventdata, handles)
% hObject    handle to editFilterSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFilterSize as text
%        str2double(get(hObject,'String')) returns contents of editFilterSize as a double


% --- Executes during object creation, after setting all properties.
function editFilterSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFilterSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnUpdateError.
function btnUpdateError_Callback(hObject, eventdata, handles)
% hObject    handle to btnUpdateError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  err0 = getappdata(handles.figMain,'err0');
  err = err0;
  udata = getappdata(handles.figMain,'udata');
  haveu = cellfun(@(x) ~isempty(x),udata);
  stacknums = find(haveu);
  errfunc = getappdata(handles.figMain,'errfunc');
  pgoptions = struct('max',length(stacknums),'progress',0);
  pgoptions = progress_bar(pgoptions);
  tic
  for thisstack = stacknums
    u = udata{thisstack};
    err(thisstack) = errfunc(thisstack,u);
    if (toc > 3 || thisstack == stacknums(end))
      cla(handles.axesErr);
      line(1:length(err0),[err0(:) err(:)],'LineStyle','none','Marker','x','Parent',handles.axesErr);
      pgoptions.progress = find(thisstack == stacknums);
      pgoptions = progress_bar(pgoptions);
      tic
    end
  end
  setappdata(handles.figMain,'err',err);

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles.output = {[],[]};
  guidata(handles.figMain,handles);
  uiresume(handles.figMain);

% --- Executes on button press in btnDone.
function btnDone_Callback(hObject, eventdata, handles)
% hObject    handle to btnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  udata = getappdata(handles.figMain,'udata');
  err = getappdata(handles.figMain,'err');
  handles.output = {udata,err};
  guidata(handles.figMain,handles);
  uiresume(handles.figMain);
  
