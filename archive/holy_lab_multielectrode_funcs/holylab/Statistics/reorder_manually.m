function varargout = reorder_manually(varargin)
% Click on columns and drag them to a new position. Click "Done" when
% satisfied.
%
%  columnOrder = reorder_manually(x)
%  columnOrder = reorder_manually(x,xerr)
%  columnOrder = reorder_manually(...,options)
%  where options may have the following fields:
%    columnOrder: a starting ordering for the columns
%    clim: color limits
%    colormap: the color map to use
%    yticklabels: a cell array of strings, used to label the y ticks
%    ytick: a vector of #s
% If you supply xerr, then x and its error is plotted in a separate figure
% when you click on it.

% REORDER_MANUALLY MATLAB code for reorder_manually.fig
%      REORDER_MANUALLY, by itself, creates a new REORDER_MANUALLY or raises the existing
%      singleton*.
%
%      H = REORDER_MANUALLY returns the handle to a new REORDER_MANUALLY or the handle to
%      the existing singleton*.
%
%      REORDER_MANUALLY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REORDER_MANUALLY.M with the given input arguments.
%
%      REORDER_MANUALLY('Property','Value',...) creates a new REORDER_MANUALLY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reorder_manually_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reorder_manually_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reorder_manually

% Last Modified by GUIDE v2.5 30-Nov-2011 12:22:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reorder_manually_OpeningFcn, ...
                   'gui_OutputFcn',  @reorder_manually_OutputFcn, ...
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


% --- Executes just before reorder_manually is made visible.
function reorder_manually_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reorder_manually (see VARARGIN)
  % Parse arguments
  x = varargin{1};
  handles.have_xerr = false;
  options = struct;
  if (length(varargin) > 1 && isnumeric(varargin{2}))
    xerr = varargin{2};
    setappdata(handles.figMain,'xerr',xerr)
    handles.have_xerr = true;
    figure
    handles.hax_xerr = gca;
  end
  if (length(varargin) > 1 + handles.have_xerr)
    options = varargin{2+handles.have_xerr};
  end
  options = default(options,...
    'columnOrder',1:size(x,2),...
    'clim',[min(x(:)) max(x(:))],...
    'colormap',jet);
  
  % Store data in the figure
  setappdata(handles.figMain,'x',x)
  setappdata(handles.figMain,'columnOrder',options.columnOrder)
  setappdata(handles.figMain,'sourceColumn',[])

  % Set up the main figure
  handles.img = imagesc(x(:,options.columnOrder),'Parent',handles.axesMatrix,...
    'ButtonDownFcn',@rm_buttondown,...
    'CDataMapping','scaled');
  set(handles.axesMatrix,'CLim',options.clim,'TickDir','out');
  colormap(handles.axesMatrix,options.colormap);
  if isfield(options,'YTick')
    set(handles.axesMatrix,'YTick',options.YTick);
  end
  if isfield(options,'YTickLabel')
    set(handles.axesMatrix,'YTickLabel',options.YTickLabel);
  end
  set(handles.figMain,'WindowButtonUpFcn',@rm_buttonup)
  
  % Choose default command line output for reorder_manually
  handles.output = hObject;
  
  % Update handles structure
  guidata(hObject, handles);
  
  % UIWAIT makes reorder_manually wait for user response (see UIRESUME)
  uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = reorder_manually_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
  if ~isempty(handles)
    varargout{1} = handles.output;
    if handles.have_xerr && ishandle(handles.hax_xerr)
      close(get(handles.hax_xerr,'Parent'))
    end
    close(handles.figMain)
  else
    varargout{1} = [];
  end


function rm_buttondown(src,~)
  hax = get(src,'Parent');
  cp = get(hax,'CurrentPoint');
  hfig = get(hax,'Parent');
  handles = guidata(hfig);
  columnIndex = round(cp(1,1));
  columnIndex0 = columnIndex;
  setappdata(handles.figMain,'sourceColumn',columnIndex);
  if handles.have_xerr
    x = getappdata(handles.figMain,'x');
    xerr = getappdata(handles.figMain,'xerr');
    columnOrder = getappdata(handles.figMain,'columnOrder');
    columnIndex = columnOrder(columnIndex);
    errorbar(handles.hax_xerr,x(:,columnIndex),xerr(:,columnIndex))
    title(handles.hax_xerr,num2str(columnIndex0))
  end
  
function rm_buttonup(src,~)
  handles = guidata(src);
  cp = get(handles.axesMatrix,'CurrentPoint');
  yl = get(handles.axesMatrix,'YLim');
  if (cp(1,2) < yl(1) || cp(1,2) > yl(2))
    % If clicked outside axes, don't do anything
    setappdata(handles.figMain,'sourceColumn',[]);
    return
  end
  destColumn = round(cp(1,1));
  sourceColumn = getappdata(handles.figMain,'sourceColumn');
  if (sourceColumn ~= destColumn)
    columnOrder = getappdata(handles.figMain,'columnOrder');
    if (sourceColumn < destColumn)
      columnOrder = columnOrder([1:sourceColumn-1 sourceColumn+1:destColumn sourceColumn destColumn+1:end]);
    else
      columnOrder = columnOrder([1:destColumn sourceColumn destColumn+1:sourceColumn-1 sourceColumn+1:end]);
    end
    setappdata(handles.figMain,'columnOrder',columnOrder)
    x = getappdata(handles.figMain,'x');
    set(handles.img,'CData',x(:,columnOrder));
  end


% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(~, ~, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles.output = getappdata(handles.figMain,'columnOrder');
  guidata(handles.figMain,handles)
  uiresume(handles.figMain)
