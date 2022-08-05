function varargout = pd_gaussian_gui(varargin)
% PD_GAUSSIAN_GUI: fitting PSFs with Gaussian phase aberrations
% 
% This GUI gives the user a chance to tweak the initial guess for the
% voltage-dependent aberration using a single Gaussian.  It then can 
% Syntax:
%   [parameters,rawphase] = pd_gaussian_gui(imseries,v,pupildata)
%   [parameters,rawphase] = pd_gaussian_gui(imseries,v,pupildata,options)
% where
%   imseries is a m-by-n-by-K set of images obtained by manipulating a
%     single actuator;
%   v is a vector of length K listing the voltage associated with each
%     image;
%   pupildata is a structure with fields H0, rho, and theta;
%   options may have the following fields:
%     paramsin: a 4-vector containing the initial guess for the parameters
%       of the Gaussian (amplitude, mean, and Omega = 1/sigma^2)
%     interactive (default true): if false, is non-blocking and the GUI can
%       be controlled remotely/automatically (see generic GUI help below)
% and
%   parameters is a structure containing the amplitude, mean, and Omega
%     of the voltage-dependent Gaussian
%   rawphase is the voltage-dependent Gaussian itself (in the pupil)
%
% Usage: adjust the amplitude and width of the gaussian using the edit
% boxes. Click inside the phase plot to set the position of the center of
% the gaussian. Test your changes by pressing "Play". When you feel like
% you've adjusted the parameters appropriately, click "Done".

% Copyright 2009 by Timothy E. Holy

% PD_GAUSSIAN_GUI M-file for pd_gaussian_gui.fig
%      PD_GAUSSIAN_GUI, by itself, creates a new PD_GAUSSIAN_GUI or raises the existing
%      singleton*.
%
%      H = PD_GAUSSIAN_GUI returns the handle to a new PD_GAUSSIAN_GUI or the handle to
%      the existing singleton*.
%
%      PD_GAUSSIAN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PD_GAUSSIAN_GUI.M with the given input arguments.
%
%      PD_GAUSSIAN_GUI('Property','Value',...) creates a new PD_GAUSSIAN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pd_gaussian_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pd_gaussian_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pd_gaussian_gui

% Last Modified by GUIDE v2.5 11-Jul-2009 05:53:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pd_gaussian_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pd_gaussian_gui_OutputFcn, ...
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


% --- Executes just before pd_gaussian_gui is made visible.
function pd_gaussian_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pd_gaussian_gui (see VARARGIN)

% Choose default command line output for pd_gaussian_gui
handles.output = {hObject};  % this will allow the figure handle to be returned for "remote control"

% Update handles structure
guidata(hObject, handles);

% Parse command-line arguments
im = varargin{1};
v = varargin{2};
pupildata = varargin{3};
if (length(varargin) > 3)
  options = varargin{4};
else
  options = struct;
end
options = default(options,'interactive',true);

K = length(v);
imsz = size3(im);
if (imsz(3) ~= K)
  error('Mismatch between size of im and the number of distinct voltages');
end

% Create the params2phi function
config = struct('mode','vgaussian singlechannel','rho',pupildata.rho,'theta',pupildata.theta,'v',v);
phimodel = phi_parametrizations(config);

% Store the data
setappdata(handles.figPDG,'phimodel',phimodel);
setappdata(handles.figPDG,'im',im);
setappdata(handles.figPDG,'v',v);
setappdata(handles.figPDG,'pupildata',pupildata);
setappdata(handles.figPDG,'options',options);
% Set up parameters
if isfield(options,'paramsin')
  gaussianP = options.paramsin;
  store_gaussianP(gaussianP,handles);
else
  mu = [0 0];
  setappdata(handles.figPDG,'mu',mu);
  update_gaussianP(handles);
end

% Show the phase
show_phase(handles)
set(handles.axesPhase,'NextPlot','replacechildren','ButtonDownFcn',@phase_click);

% Show the "movie"
play_images(handles)

% UIWAIT makes pd_gaussian_gui wait for user response (see UIRESUME)
if options.interactive
  uiwait(handles.figPDG);
end


% --- Outputs from this function are returned to the command line.
function varargout = pd_gaussian_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
  varargout = handles.output;
  options = getappdata(handles.figPDG,'options');
  if options.interactive
    close(handles.figPDG);
  end
else
  varargout = {[],[]};
end


function editAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to editAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAmplitude as text
%        str2double(get(hObject,'String')) returns contents of editAmplitude as a double


% --- Executes during object creation, after setting all properties.
function editAmplitude_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editWidth_Callback(hObject, eventdata, handles)
% hObject    handle to editWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWidth as text
%        str2double(get(hObject,'String')) returns contents of editWidth as a double


% --- Executes during object creation, after setting all properties.
function editWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnDone.
function [gaussianP,phiV0] = btnDone_Callback(hObject, eventdata, handles)
% hObject    handle to btnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  gaussianP = getappdata(handles.figPDG,'gaussianP');
  phimodel = getappdata(handles.figPDG,'phimodel');
  options = getappdata(handles.figPDG,'options');
  phiV0 = phimodel.phicoef(gaussianP);
  handles.output = {gaussianP,phiV0};
  guidata(hObject,handles);
  if options.interactive
    uiresume(handles.figPDG);
  end

% --- Executes on button press in btnPlay.
function btnPlay_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_gaussianP(handles)
show_phase(handles)
play_images(handles)

% --- Executes on button press in btnOptimize.
function btnOptimize_Callback(hObject, eventdata, handles)
% hObject    handle to btnOptimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  gaussianP = getappdata(handles.figPDG,'gaussianP');
  phimodel = getappdata(handles.figPDG,'phimodel');
  im = getappdata(handles.figPDG,'im');
  pupildata = getappdata(handles.figPDG,'pupildata');
  gaussianP = calcphi2d(im,pupildata.H0,gaussianP,phimodel);
  store_gaussianP(gaussianP,handles);
  show_phase(handles)
  play_images(handles)
  
function phase_click(hObject, eventdata)
  handles = guidata(hObject);
  cp = get(hObject,'CurrentPoint');
  setappdata(handles.figPDG,'mu',cp(1,1:2));
  update_gaussianP(handles)
  show_phase(handles)
  play_images(handles)

%% User-interaction functions
function play_images(handles)
  options = getappdata(handles.figPDG,'options');
  if ~options.interactive
    return
  end
  im = getappdata(handles.figPDG,'im');
  pupildata = getappdata(handles.figPDG,'pupildata');
  gaussianP = getappdata(handles.figPDG,'gaussianP');
  phimodel = getappdata(handles.figPDG,'phimodel');
  
  
  imsz = size3(im);
  K = imsz(3);
  %obj = getappdata(handles.figPDG,'object');
  obj = zeros(imsz(1:2));
  obj(round(imsz(1)/2),round(imsz(2)/2)) = 1;
  
  phik = phimodel.param2phi(gaussianP);
  %[val,grad,obj] = pdpenalty(phik,im,pupildata.H0);
  [Hk,sk,imc] = pd_forward_model_2d(phik,pupildata.H0,obj);
  
  for i = 1:K
    axes(handles.axesData)
    imshowsc(im(:,:,i));
    axes(handles.axesCalculated)
    imshowsc(imc(:,:,i));
    drawnow
  end
  
function show_phase(handles)
  options = getappdata(handles.figPDG,'options');
  if ~options.interactive
    return
  end
  phimodel = getappdata(handles.figPDG,'phimodel');
  gaussianP = getappdata(handles.figPDG,'gaussianP');
  pupildata = getappdata(handles.figPDG,'pupildata');

  phi = fftshift(phimodel.phicoef(gaussianP));
  mxphi = max(abs(phi(:)));
  if (mxphi == 0)
    mxphi = 1;
  end
  phinorm = phi/mxphi;
  cm = jet;
  M = size(cm,1);
  phi_int = round(phinorm * (M-1)/2 + (M+1)/2);
  phi_color = zeros([size(phi) 3]);
  for i = 1:3
    phi_color(:,:,i) = reshape(cm(phi_int(:),i),size(phi));
  end

  xmax = max(pupildata.rho(1,:));
  ymax = max(pupildata.rho(:,1));
    
  axes(handles.axesPhase);
  him = image([-1 1]*xmax,[-1 1]*ymax,phi_color);
  set(him,'HitTest','off');
  
  
%% Utility functions
function update_gaussianP(handles)
  mu = getappdata(handles.figPDG,'mu');
  gaussianP = [getValue(handles.editAmplitude) mu ...
	       1/getValue(handles.editWidth)^2];
  setappdata(handles.figPDG,'gaussianP',gaussianP);

function store_gaussianP(gaussianP,handles)
  setappdata(handles.figPDG,'mu',gaussianP(2:3));
  set(handles.editAmplitude,'String',num2str(gaussianP(1)))
  set(handles.editWidth,'String',num2str(1/sqrt(gaussianP(4))));
  setappdata(handles.figPDG,'gaussianP',gaussianP);

function sz = size3(A)
  sz = size(A);
  sz(end+1:3) = 1;

function val = getValue(h)
  val = str2double(get(h,'String'));
  



