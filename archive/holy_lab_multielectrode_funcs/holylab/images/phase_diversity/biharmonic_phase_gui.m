function varargout = biharmonic_phase_gui(varargin)
  % biharmonic_phase_gui(pupilinfo,imsz)
  
% BIHARMONIC_PHASE_GUI M-file for biharmonic_phase_gui.fig
%      BIHARMONIC_PHASE_GUI, by itself, creates a new BIHARMONIC_PHASE_GUI or raises the existing
%      singleton*.
%
%      H = BIHARMONIC_PHASE_GUI returns the handle to a new BIHARMONIC_PHASE_GUI or the handle to
%      the existing singleton*.
%
%      BIHARMONIC_PHASE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BIHARMONIC_PHASE_GUI.M with the given input arguments.
%
%      BIHARMONIC_PHASE_GUI('Property','Value',...) creates a new BIHARMONIC_PHASE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before biharmonic_phase_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to biharmonic_phase_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help biharmonic_phase_gui

% Last Modified by GUIDE v2.5 19-Nov-2009 08:02:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @biharmonic_phase_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @biharmonic_phase_gui_OutputFcn, ...
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


% --- Executes just before biharmonic_phase_gui is made visible.
function biharmonic_phase_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to biharmonic_phase_gui (see VARARGIN)

% Choose default command line output for biharmonic_phase_gui
handles.output = hObject;

pupilinfo = varargin{1};
imsz = varargin{2};
[H,rho,theta] = pupil_initialize(pupilinfo,imsz);
[act_grid,act_xy] = mirao52_layout;
act_xy = act_xy - 4.5;  % to put the center in the center of the grid

inpupil = H > 0;
th = theta(inpupil(:));
u = repmat(rho(inpupil(:)),[1 2]) .* [cos(th),sin(th)];
u = u';

setappdata(handles.figMain,'act_xy',act_xy);
setappdata(handles.figMain,'u',u);
NA_image = pupilinfo.NA/pupilinfo.M;
res_factor = pupilinfo.pixel_spacing * NA_image / pupilinfo.lambda;
th = linspace(0,2*pi,101);
ucontour = res_factor(1) * [cos(th); sin(th)];
contourIndex = convhull(u(1,:),u(2,:));
ucontour = u(:,contourIndex);
setappdata(handles.figMain,'pupilContour',ucontour);
setappdata(handles.figMain,'inpupil',inpupil);

% Lay out the axis grids
% The plot of the whole membrane
figure
hax = SplitGrid((1:7)/8,(1:7)/8,ones(1,8),ones(1,8));
hax = hax';
killFlag = isnan(act_grid);
delete(hax(killFlag))
handles.axesComplete = hax(~killFlag);
% The grid with contour plot
figure
for i = 1:52
  text(act_xy(i,2),act_xy(i,1),sprintf('%d',i))
end
hl = line(ucontour(1,:),ucontour(2,:));
handles.pupilContourLine = hl;
set(gca,'XLim',[-5 5],'YLim',[-5 5],'YDir','reverse');
% The membrane through the pupil
figure
hax = SplitGrid((1:7)/8,(1:7)/8,ones(1,8),ones(1,8));
hax = hax';
killFlag = isnan(act_grid);
delete(hax(killFlag))
handles.axesThroughPupil = hax(~killFlag);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes biharmonic_phase_gui wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = biharmonic_phase_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editR_Callback(hObject, eventdata, handles)
% hObject    handle to editR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editR as text
%        str2double(get(hObject,'String')) returns contents of editR as a double
  bhpg_update(handles)

% --- Executes during object creation, after setting all properties.
function editR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editA11_Callback(hObject, eventdata, handles)
% hObject    handle to editA11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editA11 as text
%        str2double(get(hObject,'String')) returns contents of editA11 as a double
  bhpg_update(handles)


% --- Executes during object creation, after setting all properties.
function editA11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editA21_Callback(hObject, eventdata, handles)
% hObject    handle to editA21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editA21 as text
%        str2double(get(hObject,'String')) returns contents of editA21 as a double
  bhpg_update(handles)


% --- Executes during object creation, after setting all properties.
function editA21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editA12_Callback(hObject, eventdata, handles)
% hObject    handle to editA12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editA12 as text
%        str2double(get(hObject,'String')) returns contents of editA12 as a double
  bhpg_update(handles)


% --- Executes during object creation, after setting all properties.
function editA12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editA22_Callback(hObject, eventdata, handles)
% hObject    handle to editA22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editA22 as text
%        str2double(get(hObject,'String')) returns contents of editA22 as a double
  bhpg_update(handles)


% --- Executes during object creation, after setting all properties.
function editA22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editOffset1_Callback(hObject, eventdata, handles)
% hObject    handle to editOffset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOffset1 as text
%        str2double(get(hObject,'String')) returns contents of editOffset1 as a double
  bhpg_update(handles)


% --- Executes during object creation, after setting all properties.
function editOffset1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOffset1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editOffset2_Callback(hObject, eventdata, handles)
% hObject    handle to editOffset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editOffset2 as text
%        str2double(get(hObject,'String')) returns contents of editOffset2 as a double
  bhpg_update(handles)


% --- Executes during object creation, after setting all properties.
function editOffset2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editOffset2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function bhpg_update(handles)
  R = str2double(get(handles.editR,'String'));
  A = zeros(2,2);
  A(1,1) = str2double(get(handles.editA11,'String'));
  A(2,1) = str2double(get(handles.editA21,'String'));
  A(1,2) = str2double(get(handles.editA12,'String'));
  A(2,2) = str2double(get(handles.editA22,'String'));
  xi0 = zeros(1,2);
  xi0(1) = str2double(get(handles.editOffset1,'String'));
  xi0(2) = str2double(get(handles.editOffset2,'String'));

  % Display the membrane shape
  act_xi = getappdata(handles.figMain,'act_xy');
  for i = 1:size(act_xi,1)
    b = biharmonic_pointsource(100,act_xi(i,:)/R);
    imagesc(b,'Parent',handles.axesComplete(i));
  end
  set(handles.axesComplete,'XTick',[],'YTick',[]);
  
  % Draw the pupil contour
  pupilContour = getappdata(handles.figMain,'pupilContour');
  xi = A*pupilContour + repmat(xi0(:),1,size(pupilContour,2));
  set(handles.pupilContourLine,'XData',xi(1,:),'YData',xi(2,:));
  
  % Display the phase
  u = getappdata(handles.figMain,'u');
  inpupil = getappdata(handles.figMain,'inpupil');
  xi = A*u + repmat(xi0(:),1,size(u,2));
  xi = xi';
  for i = 1:size(act_xi,1)
    ph = biharmonic_pointsource(xi,act_xi(i,:),R);
    ph = ph - mean(ph(:));  % Subtract the piston shift
    phim = zeros(size(inpupil));
    phim(inpupil) = ph;
    imagesc(fftshift(phim),'Parent',handles.axesThroughPupil(i));
  end
  set(handles.axesThroughPupil,'XTick',[],'YTick',[]);
  