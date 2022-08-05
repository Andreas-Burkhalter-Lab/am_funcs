function varargout = triage_templates(varargin)
% TRIAGE_TEMPLATES Eliminate redundant templates
% Syntax:
%    triage_templates(filein)
%    triage_templates(filein,options)
% where
%   filein is the name of a fine_cluster file. By default, at the end this
%     fine_cluster file will be moved to fine_cluster0, and the new file
%     will take its place
%   options is a structure which may have the following fields:
%      thresh (default 0.8): threshold for the dot product between
%        normalized, smoothed template pairs needed to trigger examination
%      lowpass (default 0.2): cutoff (in Nyquist units) of lowpassing used
%        to smooth the templates.  This helps insure that noise does not
%        dominate the comparison between templates.


% TRIAGE_TEMPLATES M-file for triage_templates.fig
%      TRIAGE_TEMPLATES, by itself, creates a new TRIAGE_TEMPLATES or raises the existing
%      singleton*.
%
%      H = TRIAGE_TEMPLATES returns the handle to a new TRIAGE_TEMPLATES or the handle to
%      the existing singleton*.
%
%      TRIAGE_TEMPLATES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRIAGE_TEMPLATES.M with the given input arguments.
%
%      TRIAGE_TEMPLATES('Property','Value',...) creates a new TRIAGE_TEMPLATES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before triage_templates_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to triage_templates_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help triage_templates

% Last Modified by GUIDE v2.5 19-Jul-2007 22:53:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @triage_templates_OpeningFcn, ...
                   'gui_OutputFcn',  @triage_templates_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
 if (nargin > 2) && ischar(varargin{1})
     gui_State.gui_Callback = str2func(varargin{1});
 end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before triage_templates is made visible.
function triage_templates_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to triage_templates (see VARARGIN)
  if ~isunix
    error('This works only on UNIX platforms');
  end
  filename = varargin{1};
  if (length(varargin) < 2)
    options = struct;
  else
    options = varargin{1};
  end
  options = default(options,'thresh',0.8);
  options = default(options,'lowpass',0.2);
  
  load('-mat',filename);
  n_templates = length(templates);
  
  % Create smoothed, normalized templates (so that noise doesn't dominate
  % estimates of correlation)
  [b,a] = butter(2,options.lowpass,'low');
  tpl_sm = cell(1,n_templates);
  for i = 1:n_templates
    tpl_sm{i} = filtfilt(b,a,templates{i});  % smoothed template
    tpl_sm{i} = tpl_sm{i} / sqrt(sum(tpl_sm{i}.^2));  % normalized, smoothed
  end
  
  % Compute the dot product between all pairs
  cth = zeros(n_templates,n_templates);  % "cos(theta)"
  for i = 1:n_templates
    for j = i+1:n_templates
      cth(i,j) = sum(tpl_sm{i}.*tpl_sm{j});
    end
  end
  
  % Sort them into decreasing order, in a way that we still know which
  % templates are in each pair
  cthposIndex = find(cth > options.thresh);
  [cthposI,cthposJ] = find(cth > options.thresh);  
  [sort_cth,sortIndex] = sort(cth(cthposIndex),'descend');
  cthpos = [cthposI(:) cthposJ(:)];
  cthpos = cthpos(sortIndex,:);
  
  isalive = true(1,n_templates);
  currentIndex = 1;
  
  setappdata(handles.figTriageTemplates,'templates',templates);
  setappdata(handles.figTriageTemplates,'templ_pair',cthpos);
  setappdata(handles.figTriageTemplates,'isalive',isalive);
  setappdata(handles.figTriageTemplates,'currentIndex',currentIndex);
  setappdata(handles.figTriageTemplates,'filename',filename);
  setappdata(handles.figTriageTemplates,'template_length',length(templates{1}));
  
  set(handles.axesTemplates,'ButtonDownFcn',{@trtmpl_zoom,handles},...
    'NextPlot','replacechildren');
  
  trtempl_draw(handles);
  
% Choose default command line output for triage_templates
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes triage_templates wait for user response (see UIRESUME)
% uiwait(handles.figTriageTemplates);


% --- Outputs from this function are returned to the command line.
function varargout = triage_templates_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Kill1.
function Kill1_Callback(hObject, eventdata, handles)
% hObject    handle to Kill1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  trtempl_kill(handles,1);

% --- Executes on button press in NextBtn.
function NextBtn_Callback(hObject, eventdata, handles)
% hObject    handle to NextBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  trtempl_inc(handles);
  trtempl_draw(handles);


% --- Executes on button press in Kill2.
function Kill2_Callback(hObject, eventdata, handles)
% hObject    handle to Kill2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  trtempl_kill(handles,2);
  

function trtempl_kill(handles,killindex)
  templ_pair = getappdata(handles.figTriageTemplates,'templ_pair');
  isalive = getappdata(handles.figTriageTemplates,'isalive');
  currentIndex = getappdata(handles.figTriageTemplates,'currentIndex');
  
  isalive(templ_pair(currentIndex,killindex)) = false;
  setappdata(handles.figTriageTemplates,'isalive',isalive);

  trtempl_inc(handles);
  trtempl_draw(handles);
  
function trtempl_inc(handles)
  % Increment the currentIndex, skipping pairs that have at least one
  % non-alive member
  templ_pair = getappdata(handles.figTriageTemplates,'templ_pair');
  isalive = getappdata(handles.figTriageTemplates,'isalive');
  currentIndex = getappdata(handles.figTriageTemplates,'currentIndex')+1;
  n_pairs = size(templ_pair,1);
  
  if (currentIndex > n_pairs)
    trtempl_save(handles);
  end
  c_pair = templ_pair(currentIndex,:);
  while (~isalive(c_pair(1)) || ~isalive(c_pair(2))) && currentIndex < n_pairs
    currentIndex = currentIndex+1;
    c_pair = templ_pair(currentIndex,:);
  end
  if (currentIndex > n_pairs)
    trtempl_save(handles);
  end
  setappdata(handles.figTriageTemplates,'currentIndex',currentIndex);
  currentIndex
  
function trtempl_draw(handles)
  templ_pair = getappdata(handles.figTriageTemplates,'templ_pair');
  isalive = getappdata(handles.figTriageTemplates,'isalive');
  currentIndex = getappdata(handles.figTriageTemplates,'currentIndex');
  templates = getappdata(handles.figTriageTemplates,'templates');
  n_pairs = size(templ_pair,1);
  c_pair = templ_pair(currentIndex,:);
  
  set(handles.textProgress,'String',sprintf('%d%%',round(100*currentIndex/n_pairs)));
  t1 = templates{c_pair(1)};
  t1 = t1 / sqrt(sum(t1.^2));
  t2 = templates{c_pair(2)};
  t2 = t2 / sqrt(sum(t2.^2));
  hline = plot(handles.axesTemplates,[t1 t2]);
  set(hline,'HitTest','off')
  axis tight
  
function trtempl_save(handles)
  isalive = getappdata(handles.figTriageTemplates,'isalive');
  % Get the original filename
  filename = getappdata(handles.figTriageTemplates,'filename');
  % Find all the files with name = 'filename%d', where %d is an integer.
  % We're going to move the input file to an unclaimed file of that format
  oldfilenames = dirbyname([filename '*']);
  currentBackup = 0;
  for i = 1:length(oldfilenames)
    thisBackup = sscanf(oldfilenames{i},[filename '%d']);
    if ~isempty(thisBackup)
      currentBackup = max(currentBackup,thisBackup+1);
    end
  end
  load('-mat',filename);
  [status,result] = system(['mv ' filename ' ' filename num2str(currentBackup)]);
  if (status == 0)
    fineClusters = fineClusters(isalive);
    templates = templates(isalive);
    save(filename,'channels','fineClusterOptions','fineClusters',...
      'medv','rawClusterFilename','rawClusterStartIndex','templates',...
      'thresh','-mat');
  else
    errdlg('Had trouble moving the old file out of the way. Is this is permissions problem? Fix the problem, and then click "Next" again',...
      'File I/O trouble');
  end
  
  
function trtmpl_zoom(hObject,eventdata,handles)
  cp = get(handles.axesTemplates,'CurrentPoint');
  seltype = get(handles.figTriageTemplates,'SelectionType');
  cp = cp(1,1);
  xlim = get(handles.axesTemplates,'XLim');
  dx = diff(xlim);
  template_length = getappdata(handles.figTriageTemplates,'template_length');
  
  switch(seltype)
    case 'normal'
      % zoom in
      xlim = cp + dx/4*[-1 1];
    case 'alt'
      xlim = cp + dx*[-1 1];
  end
  xlim(1) = max(1,xlim(1));
  xlim(2) = min(template_length,xlim(2));
  set(handles.axesTemplates,'XLim',xlim);
  