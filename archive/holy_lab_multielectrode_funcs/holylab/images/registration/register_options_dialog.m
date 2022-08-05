function varargout = register_options_dialog(varargin)
% REGISTER_OPTIONS_DIALOG M-file for register_options_dialog.fig
%      REGISTER_OPTIONS_DIALOG, by itself, creates a new REGISTER_OPTIONS_DIALOG or raises the existing
%      singleton*.
%
%      H = REGISTER_OPTIONS_DIALOG returns the handle to a new REGISTER_OPTIONS_DIALOG or the handle to
%      the existing singleton*.
%
%      REGISTER_OPTIONS_DIALOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_OPTIONS_DIALOG.M with the given input arguments.
%
%      REGISTER_OPTIONS_DIALOG('Property','Value',...) creates a new REGISTER_OPTIONS_DIALOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before register_options_dialog_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_options_dialog_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_options_dialog

% Last Modified by GUIDE v2.5 08-Jan-2011 06:21:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @register_options_dialog_OpeningFcn, ...
                   'gui_OutputFcn',  @register_options_dialog_OutputFcn, ...
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


% --- Executes just before register_options_dialog is made visible.
function register_options_dialog_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_options_dialog (see VARARGIN)

defoptions = {'presmooth',[0 0 0],...
  'refit',0,...
  'allowed_initialguess',{'Unwarped','Previous solution','Temporal interpolation'},...
  'allowed_method',{'Phase correlation','Multigrid least-squares'},...
  'initialguess','Unwarped',...
  'method','Phase correlation',...
  'lambda',0,...
  'wcycle',true,...
  'gap_data',0};
options = varargin{1};
options = default(options,defoptions{:});

rmod_set_presmooth(handles,options.presmooth);
issmoothing = any(options.presmooth > 0);
offon = {'off','on'};
set(handles.checkboxRefitUnsmoothed,'Value',options.refit & issmoothing,'Enable',offon{issmoothing+1});
igvalue = strmatch(options.initialguess,options.allowed_initialguess,'exact');
if isempty(igvalue)
  error(['Initial guess "' options.initialguess '" is not one of the allowed choices']);
end
set(handles.popupmenuInitialGuess,'String',options.allowed_initialguess,...
  'Value',igvalue);
mvalue = strmatch(options.method,options.allowed_method,'exact');
if isempty(mvalue)
  error(['Method "' options.initialguess '" is not one of the allowed choices']);
end
set(handles.popupmenuMethod,'String',options.allowed_method,...
  'Value',mvalue);
set(handles.editLambda,'String',num2str(options.lambda));
set(handles.checkboxWCycle,'Value',options.wcycle);
set(handles.editGapData,'String',num2str(options.gap_data));
if ~isnan(options.mspixval)
  set(handles.textMeanSquareImageValue,'String',sprintf('%0.2e',options.mspixval));
else
  set([handles.textMsqPV handles.textMeanSquareImageValue],'Visible','off');
end
setappdata(handles.figMain,'options',options);
rmod_update_gridstats(handles)
rod_setstate(handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes register_options_dialog wait for user response (see UIRESUME)
uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = register_options_dialog_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
  if isempty(handles)
    % User closed the dialog box
    varargout = {0};
    return
  end
  canceled = getappdata(handles.figMain,'canceled');
  if canceled
    varargout = {0};
  else
    allowed_initialguess = get(handles.popupmenuInitialGuess,'String');
    allowed_method = get(handles.popupmenuMethod,'String');
    output = struct('initialguess',allowed_initialguess{get(handles.popupmenuInitialGuess,'Value')},...
      'method',allowed_method{get(handles.popupmenuMethod,'Value')},...
      'presmooth',rmod_get_presmooth(handles),...
      'refit',get(handles.checkboxRefitUnsmoothed,'Value'),...
      'lambda',str2double(get(handles.editLambda,'String')),...
      'gap_data',str2double(get(handles.editGapData,'String')),...
      'wcycle',get(handles.checkboxWCycle,'Value'),...
      'createmask',get(handles.checkboxCreateMask,'Value'));
    varargout = {output};
  end
  close(handles.figMain)


function rmod_set_presmooth(handles,val)
  set(handles.editPreSmooth,'String',sprintf('[%g %g %g]',val))
  
function val = rmod_get_presmooth(handles)
  str = get(handles.editPreSmooth,'String');
  str = str(2:end-1);
  val = sscanf(str,'%g')';

function rmod_update_gridstats(handles)
  options = getappdata(handles.figMain,'options');
  gridsize = options.gridsize;
  image_gridkeep = all(gridsize >= options.min_pixels,2);
  u_gridkeep = all(gridsize >= options.min_g_pixels,2);
  image_gridsize = gridsize(image_gridkeep,:);
  u_gridsize = gridsize(u_gridkeep,:);
  gap_data = str2double(get(handles.editGapData,'String'));
  u_gridsize = u_gridsize(1+gap_data:end,:);
  max_grids = min(size(image_gridsize,1),size(u_gridsize,1));
  % Only show the full image grid if the method is multigrid
  methods = get(handles.popupmenuMethod,'String');
  method = methods{get(handles.popupmenuMethod,'Value')};
  if strcmp(method,'Multigrid least-squares')
    set(handles.textImageGrid,'String',num2str(image_gridsize(1:max_grids,:)))
    set(handles.textUGrid,'String',num2str(u_gridsize(1:max_grids,:)))
  else
    set(handles.textImageGrid,'String',num2str(image_gridsize(1,:)))
    set(handles.textUGrid,'String',num2str(u_gridsize))
  end

function rod_setstate(handles)
  methods = get(handles.popupmenuMethod,'String');
  method = methods{get(handles.popupmenuMethod,'Value')};
  offon = {'off','on'};
  h = [handles.checkboxWCycle handles.checkboxCreateMask handles.textPreSmooth ...
    handles.editPreSmooth handles.checkboxRefitUnsmoothed];
  set(h,'Visible',offon{strcmp(method,'Multigrid least-squares')+1});
  rmod_update_gridstats(handles)
  
  
function editPreSmooth_Callback(hObject, eventdata, handles)
% hObject    handle to editPreSmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPreSmooth as text
%        str2double(get(hObject,'String')) returns contents of editPreSmooth as a double
  val = rmod_get_presmooth(handles);
  if any(val > 0)
    set(handles.checkboxRefitUnsmoothed,'Enable','on')
  else
    set(handles.checkboxRefitUnsmoothed,'Enable','off','Value',0)
  end

% --- Executes during object creation, after setting all properties.
function editPreSmooth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPreSmooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxRefitUnsmoothed.
function checkboxRefitUnsmoothed_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxRefitUnsmoothed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxRefitUnsmoothed



function editGapData_Callback(hObject, eventdata, handles)
% hObject    handle to editGapData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editGapData as text
%        str2double(get(hObject,'String')) returns contents of editGapData as a double
  rmod_update_gridstats(handles)

% --- Executes during object creation, after setting all properties.
function editGapData_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editGapData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuInitialGuess.
function popupmenuInitialGuess_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInitialGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInitialGuess contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInitialGuess


% --- Executes during object creation, after setting all properties.
function popupmenuInitialGuess_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInitialGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editLambda_Callback(hObject, eventdata, handles)
% hObject    handle to editLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLambda as text
%        str2double(get(hObject,'String')) returns contents of editLambda as a double


% --- Executes during object creation, after setting all properties.
function editLambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxWCycle.
function checkboxWCycle_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxWCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxWCycle


% --- Executes on button press in btnOK.
function btnOK_Callback(hObject, eventdata, handles)
% hObject    handle to btnOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figMain,'canceled',false);
uiresume(handles.figMain)

% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.figMain,'canceled',true);
uiresume(handles.figMain)


% --- Executes on selection change in popupmenuMethod.
function popupmenuMethod_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMethod
  rod_setstate(handles);


% --- Executes during object creation, after setting all properties.
function popupmenuMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxCreateMask.
function checkboxCreateMask_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxCreateMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
