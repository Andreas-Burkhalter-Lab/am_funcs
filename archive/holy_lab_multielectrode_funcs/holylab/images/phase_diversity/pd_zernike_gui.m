function varargout = pd_zernike_gui(varargin)
% PD_ZERNIKE_GUI: fitting PSFs with gaussian phase aberrations
% Syntax:
%   Zvalue = pd_zernike_gui(imbase,imseries,pupildata,Zindex,Zvalue0,options)
% where
%   imbase is a m-by-n "reference image" of the object (should be
%     approximately unaberrated)
%   imseries is a m-by-n-by-K set of images (e.g., those obtained by
%     manipulating a single actuator);
%   pupildata is a structure with fields H0, rho, and theta;
%   Zindex is a vector of P Zernike indices (single-index numbering scheme)
%     specifying which basis functions to use
%   Zvalue0 is a K-by-P matrix giving the starting guess for the P Zernike
%     coefficients for each of the K images compared to imbase.
%   options is a structure which may have the following fields:
%     Zxreg (default [0 0]) is a 1-by-2 vector giving the first-order Zernike coefficients
%       needed to align imbase to imseries
%     interactive (default true): if false, is non-blocking and the GUI can
%       be controlled remotely/automatically (see generic GUI help below)
% and
%   Zvalue is the optimized set of Zernike coefficients for each image, a
%     K-by-P matrix.
%
% Usage notes: if there is some drift in the unaberrated image, you can
% just subtract off the corresponding first-order Zernikes from the result.
  
% Copyright 2009 by Timothy E. Holy

% PD_ZERNIKE_GUI M-file for pd_zernike_gui.fig
%      PD_ZERNIKE_GUI, by itself, creates a new PD_ZERNIKE_GUI or raises the existing
%      singleton*.
%
%      H = PD_ZERNIKE_GUI returns the handle to a new PD_ZERNIKE_GUI or the handle to
%      the existing singleton*.
%
%      PD_ZERNIKE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PD_ZERNIKE_GUI.M with the given input arguments.
%
%      PD_ZERNIKE_GUI('Property','Value',...) creates a new PD_ZERNIKE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pd_zernike_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pd_zernike_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pd_zernike_gui

% Last Modified by GUIDE v2.5 11-Jul-2009 06:33:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pd_zernike_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pd_zernike_gui_OutputFcn, ...
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


% --- Executes just before pd_zernike_gui is made visible.
function pd_zernike_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pd_zernike_gui (see VARARGIN)

% Choose default command line output for pd_zernike_gui
handles.output = {hObject};

% Update handles structure
guidata(hObject, handles);

% Parse command-line arguments
imbase = varargin{1};
imseries = varargin{2};
pupildata = varargin{3};
Zindex = varargin{4};
Zvalue0 = varargin{5};
if nargin > 5
  options = varargin{6};
else
  options = struct;
end
options = default(options,'Zxreg',[0 0],'interactive',true);
Zxreg = options.Zxreg;

zernv = zernike_values(pupildata.rho,pupildata.theta,Zindex);

setappdata(handles.figPDZ,'imbase',imbase);
setappdata(handles.figPDZ,'imseries',imseries)
setappdata(handles.figPDZ,'pupildata',pupildata);
setappdata(handles.figPDZ,'zernv',zernv);
setappdata(handles.figPDZ,'Zvalue',Zvalue0);
setappdata(handles.figPDZ,'Zxreg',Zxreg);
setappdata(handles.figPDZ,'Zindex',Zindex);
setappdata(handles.figPDZ,'recalc',true);
setappdata(handles.figPDZ,'options',options);

% Display
show_Zernikes(handles);

% UIWAIT makes pd_zernike_gui wait for user response (see UIRESUME)
if options.interactive
  uiwait(handles.figPDZ);
end


% --- Outputs from this function are returned to the command line.
function varargout = pd_zernike_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
  varargout = handles.output;
  options = getappdata(handles.figPDZ,'options');
  if options.interactive
    close(handles.figPDZ);
  end
else
  varargout = {[]};
end

% --- Executes on button press in btnDone.
function Zvalue = btnDone_Callback(hObject, eventdata, handles)
% hObject    handle to btnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Zvalue = getappdata(handles.figPDZ,'Zvalue');
options = getappdata(handles.figPDZ,'options');
handles.output = {Zvalue};
guidata(hObject,handles);
if options.interactive
  uiresume(handles.figPDZ);
end


% --- Executes on button press in btnMedian.
function btnMedian_Callback(hObject, eventdata, handles)
% hObject    handle to btnMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  Zvalue = getappdata(handles.figPDZ,'Zvalue');
  medSize = str2double(get(handles.editMedianSize,'String'));
  medSize = round(medSize);
  if (medSize < 3)
    medSize = 3;
    set(handles.editMedianSize,'String',num2str(medSize));
  end
  Zvalue = medfilt1improved(Zvalue,medSize);
  setappdata(handles.figPDZ,'Zvalue',Zvalue);
  setappdata(handles.figPDZ,'recalc',true);
  show_Zernikes(handles);

% --- Executes on button press in btnOptimize.
function fval = btnOptimize_Callback(hObject, eventdata, handles)
% hObject    handle to btnOptimize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  imbase = getappdata(handles.figPDZ,'imbase');
  imseries = getappdata(handles.figPDZ,'imseries');
  pupildata = getappdata(handles.figPDZ,'pupildata');
  Zvalue = getappdata(handles.figPDZ,'Zvalue');
  Zxreg = getappdata(handles.figPDZ,'Zxreg');
  options = getappdata(handles.figPDZ,'options');

  % Set up a Zernike model of the phases
  params.mode = 'Zn pair';
  params.Zindex = getappdata(handles.figPDZ,'Zindex');
  params.rho = pupildata.rho;
  params.theta = pupildata.theta;
  zf = phi_parametrizations(params);
  zf.normalize_grad = false;
  zf.tol = 1e-8;
  zf.iter_max = 1000;

  if options.interactive
    % Turn the current Zvalue plot into dashed lines
    axes(handles.axesZernike)
    hZold = findobj('type','line');
    set(hZold,'LineStyle','--');
    drawnow
  end
  
  % Optimize the Zernikes
  imtot = imbase;
  imtot(:,:,2:size(imseries,3)+1) = imseries;
  Zxregrep = repmat(Zxreg,size(Zvalue,1),1);
  Zvalue(:,[1 2]) = Zvalue(:,[1 2]) + Zxregrep;  % compensate for shift
  [Zvalue,fval] = calcphi2d(imtot,pupildata.H0,Zvalue,zf);
  Zvalue(:,[1 2]) = Zvalue(:,[1 2]) - Zxregrep;
  setappdata(handles.figPDZ,'Zvalue',Zvalue);
  
  if options.interactive
    % Plot the new values and, after a delay, delete the old lines
    hold on
    plot(Zvalue)
    hold off
    axis tight
    pause(2)
    if all(ishandle(hZold))
      delete(hZold)
    end
  end

% --- Executes on button press in btnPlay.
function btnPlay_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  recalc = getappdata(handles.figPDZ,'recalc');
  if recalc
    calculate_images(handles)
  end
  play_images(handles)

  

function editMedianSize_Callback(hObject, eventdata, handles)
% hObject    handle to editMedianSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMedianSize as text
%        str2double(get(hObject,'String')) returns contents of editMedianSize as a double


% --- Executes during object creation, after setting all properties.
function editMedianSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMedianSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnLinRegress.
function btnLinRegress_Callback(hObject, eventdata, handles)
% hObject    handle to btnLinRegress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  Zvalue = getappdata(handles.figPDZ,'Zvalue');
  n = size(Zvalue,1);
  [Zv,Zs] = linregress((1:n)',Zvalue);
  Zvalue = (1:n)' * Zv + repmat(Zs,n,1);
  setappdata(handles.figPDZ,'Zvalue',Zvalue);
  show_Zernikes(handles)
  
%% User-interaction functions
function play_images(handles)
  imseries = getappdata(handles.figPDZ,'imseries');
  imscalc = getappdata(handles.figPDZ,'imscalc');
  
  imsz = size3(imseries);
  K = imsz(3);
  
  for i = 1:K
    axes(handles.axesData)
    imshowsc(imseries(:,:,i));
    axes(handles.axesCalculated)
    imshowsc(imscalc(:,:,i));
    drawnow
  end
  
function show_Zernikes(handles)
  Zvalue = getappdata(handles.figPDZ,'Zvalue');
  axes(handles.axesZernike);
  plot(Zvalue);
  axis tight
  
%% Utility functions
function calculate_images(handles)
  imbase = getappdata(handles.figPDZ,'imbase');
  imseries = getappdata(handles.figPDZ,'imseries');
  imtot = imbase;
  imtot(:,:,2:size(imseries,3)+1) = imseries;
  zernv = getappdata(handles.figPDZ,'zernv');
  Zvalue = getappdata(handles.figPDZ,'Zvalue');
  Zxreg = getappdata(handles.figPDZ,'Zxreg');
  Zxregrep = repmat(Zxreg,size(Zvalue,1),1);
  Zvalue(:,[1 2]) = Zvalue(:,[1 2]) + Zxregrep;
  phik = zeros(size(imtot));
  for frameIndex = 1:size(imseries,3)
    for Zindex = 1:size(zernv,3)
      phik(:,:,frameIndex+1) = phik(:,:,frameIndex+1) + zernv(:,:,Zindex) * ...
        Zvalue(frameIndex,Zindex);
    end
  end
  pupildata = getappdata(handles.figPDZ,'pupildata');
  [val,grad,obj] = pdpenalty(phik,imtot,pupildata.H0);
  [Hk,sk,imc] = pd_forward_model_2d(phik,pupildata.H0,obj);
  setappdata(handles.figPDZ,'imscalc',imc(:,:,2:end));
  setappdata(handles.figPDZ,'recalc',false);

function sz = size3(A)
  sz = size(A);
  sz(end+1:3) = 1;

function val = getValue(h)
  val = str2double(get(h,'String'));
  


