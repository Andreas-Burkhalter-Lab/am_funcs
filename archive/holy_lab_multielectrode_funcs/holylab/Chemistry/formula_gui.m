function varargout = formula_gui(varargin)
% FORMULA_GUI: inspect mass spec data related to chemical formula
% Syntax:
%    formula_gui(mass,isotopeinfo,nI,nIcov,options)

% FORMULA_GUI M-file for formula_gui.fig
%      FORMULA_GUI, by itself, creates a new FORMULA_GUI or raises the existing
%      singleton*.
%
%      H = FORMULA_GUI returns the handle to a new FORMULA_GUI or the handle to
%      the existing singleton*.
%
%      FORMULA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FORMULA_GUI.M with the given input arguments.
%
%      FORMULA_GUI('Property','Value',...) creates a new FORMULA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before formula_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to formula_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help formula_gui

% Last Modified by GUIDE v2.5 18-Sep-2009 06:35:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @formula_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @formula_gui_OutputFcn, ...
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


% --- Executes just before formula_gui is made visible.
function formula_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to formula_gui (see VARARGIN)

% Choose default command line output for formula_gui
handles.output = hObject;

mass = varargin{1};
isotopeinfo = varargin{2};
nI = varargin{3};
nIcov = varargin{4};
options = struct;
if (length(varargin) > 4)
  options = varargin{5};
end
options = default(options,'mass_fractional_error',1e-5);

setappdata(handles.figMain,'mass',mass);
setappdata(handles.figMain,'isotopeinfo',isotopeinfo);
setappdata(handles.figMain,'nI',nI);
setappdata(handles.figMain,'nIcov',nIcov);
set(handles.editMassAccuracy,'String',num2str(options.mass_fractional_error));

% Update handles structure
guidata(hObject, handles);
update_display(handles);

% UIWAIT makes formula_gui wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = formula_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editMassAccuracy_Callback(hObject, eventdata, handles)
% hObject    handle to editMassAccuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMassAccuracy as text
%        str2double(get(hObject,'String')) returns contents of editMassAccuracy as a double
update_display(handles);

% --- Executes during object creation, after setting all properties.
function editMassAccuracy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMassAccuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function update_display(handles)
  hax = findobj(handles.figMain,'type','axes');
  hax_extra = setdiff(hax,handles.axesPlot);
  delete(hax_extra);
  ops = struct('mass_fractional_error',str2double(get(handles.editMassAccuracy,'String')));
  mass = getappdata(handles.figMain,'mass');
  isotopeinfo = getappdata(handles.figMain,'isotopeinfo');
  nI = getappdata(handles.figMain,'nI');
  nIcov = getappdata(handles.figMain,'nIcov');
  [f,me] = mass2formula(mass,isotopeinfo,ops);
  [sme,sortOrder] = sort(abs(me));
  me = me(sortOrder);
  f = f(sortOrder,:);
  chi2 = formula_chisq(f,nI,nIcov);
  axes(handles.axesPlot)
  hl = plot(me/mass,'b');
  ylabel(handles.axesPlot,'Fractional mass error','Color','b');
  xlabel(handles.axesPlot,'Formula candidate','Color','k');
  ax2 = axes('Position',get(handles.axesPlot,'Position'),...
    'XAxisLocation','bottom',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k');
  hl(2) = line(1:length(chi2),chi2,'Parent',ax2,'Color','r');
  ylabel(ax2,'Isotopologue fitting error','Color','r');
  setappdata(handles.figMain,'formula',f);
  set(hl,'HitTest','off')
  set([handles.axesPlot ax2],'ButtonDownFcn',@update_formula_string,'Box','off','TickDir','out');

function update_formula_string(src,event)
  handles = guidata(src);
  cp = get(src,'CurrentPoint');
  fIndex = round(cp(1));
  formula = getappdata(handles.figMain,'formula');
  formula = formula(fIndex,:);
  isotopeinfo = getappdata(handles.figMain,'isotopeinfo');
  fstr = '';
  for i = 1:length(formula)
    if (formula(i) ~= 0)
      fstr = [fstr isotopeinfo(i).symbol num2str(formula(i)) ' '];
    end
  end
  set(handles.textFormulaString,'String',fstr);

