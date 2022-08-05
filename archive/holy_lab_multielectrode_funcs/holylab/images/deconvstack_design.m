function varargout = deconvstack_design(varargin)
% DECONVSTACK_DESIGN: viewing tool for 4D image stacks
% Syntax:
%   deconvstack_design
% The basic call
%   deconvstack_design(basename)
% Start with a particular file open
%   deconvstack_design(basename,'checkboxUpdateT',0)
% Start with a file, and don't ever draw the time response window (it's
% the drawing of this window that, over the network especially, takes the
% most time)

% DECONVSTACK_DESIGN M-file for deconvstack_design.fig
%      DECONVSTACK_DESIGN, by itself, creates a new DECONVSTACK_DESIGN or raises the existing
%      singleton*.
%
%      H = DECONVSTACK_DESIGN returns the handle to a new DECONVSTACK_DESIGN or the handle to
%      the existing singleton*.
%
%      DECONVSTACK_DESIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DECONVSTACK_DESIGN.M with the given input arguments.
%
%      DECONVSTACK_DESIGN('Property','Value',...) creates a new DECONVSTACK_DESIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before deconvstack_design_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to deconvstack_design_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help deconvstack_design

% Last Modified by GUIDE v2.5 29-Mar-2006 06:20:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deconvstack_design_OpeningFcn, ...
                   'gui_OutputFcn',  @deconvstack_design_OutputFcn, ...
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


% --- Executes just before deconvstack_design is made visible.
function deconvstack_design_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to deconvstack_design (see VARARGIN)

% Choose default command line output for deconvstack_design
handles.output = hObject;

% Draw the play/stepping icons
icons = swicons;
buttonTag = {'btnPlayUp','btnStepUp','btnStepDown','btnPlayDown',...
  'btnStop'};
iconField = {'up','upsm','downsm','down',...
  'stop'};
for i = 1:length(buttonTag)
  set(handles.(buttonTag{i}),'CData',icons.(iconField{i}));
end

% Set up a stack movie, if a basename was supplied; otherwise, set up a
% blank movie (user will have to click "File->Open")
if ~isempty(varargin)
  smm = stackmm(varargin{1});
else
  smm = [];
end
setappdata(handles.figDSD,'stack',smm);

% If there are additional arguments, parse them
% Currently this only supports setting checkboxUpdateT and
% checkboxUpdateZ values.
if (length(varargin) > 1)
  carg = 2;
  while (carg < length(varargin))
    hname = varargin{carg};
    hvalue = varargin{carg+1};
    set(handles.(hname),'Value',hvalue);
    carg = carg+2;
  end
end

% Set up blank images
handles.imageXY = image('parent',handles.axesXY,...
                        'ButtonDownFcn',@mousedownXY);
handles.imageXZ = image('parent',handles.axesXZ);
set([handles.axesXY handles.axesXZ],...
  'CLim',[500 3000]);
set([handles.imageXY],...
  'CDataMapping','scaled','EraseMode','none');
set([handles.imageXZ],...
  'CDataMapping','scaled','EraseMode','normal');

% Install mouse handler & draglines
%install_mouse_event_handler(handles.axesXY,'down',@mousedownXY);
%set(handles.axesXY,'ButtonDownFcn',@mousedownXY);
%install_mouse_event_handler(handles.axesXY,'move',@mousemoveXY);
%install_mouse_event_handler(handles.axesXY,'up',@mouseupXY);

% Place holders for draglines
handles.dragZ1 = [];
handles.dragZ2 = [];
handles.selbox = [];

% Update handles structure
guidata(hObject, handles);

% Initialize default locations
dsd_initialize_positions(handles);
colormap(gray(256));

% Display images data
dsd_display(handles,{'XY','XZ'});

% UIWAIT makes deconvstack_design wait for user response (see UIRESUME)
% uiwait(handles.figDSD);


% --- Outputs from this function are returned to the command line.
function varargout = deconvstack_design_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function handled = mousedownXY(sender,event)
  hax = get(sender,'Parent');
  hfig = get(hax,'Parent');
  handles = guidata(hfig);
  if ishandle(handles.selbox)
    delete(handles.selbox);
  end
  cpstart = get(hax,'CurrentPoint');
  cpstart = round(cpstart(1,1:2));
  rb = rbbox;
  cpend = get(hax,'CurrentPoint');
  cpend = round(cpend(1,1:2));
  p1 = min(cpstart,cpend);
  offset = abs(cpstart-cpend);
  setappdata(handles.figDSD,'selectXY',[p1 offset]);
  set(handles.btnCalc,'Enable','on')
  x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
  y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
  handles.selbox = line(x,y,'Color','r','LineWidth',3);
  cp = getappdata(handles.figDSD,'XYZT');
  cp = [round(p1+offset/2) 1 1];
  setappdata(handles.figDSD,'XYZT',cp);
  guidata(hfig,handles);
  dsd_display(handles,{'XZ'});
  dsd_updatenpix(handles);
  handled = 1;

function handled = mouseupXY(sender,event)
  hax = sender;
  cpfinal = get(hax,'CurrentPoint');
  cpfinal = round(cpfinal(1,1:2));
  hfig = get(hax,'Parent');
  cpstart = getappdata(hfig,'rbstart');
  handled = 1;

% Callback when user clicks on an image, to choose a new set of coordinates
function dsd_choose(src,evt,hfig,coord)
  hax = get(src,'Parent');
  handles = guidata(hfig);
  cp = get(hax,'CurrentPoint');
  cp = round(cp(1,1:2));
  cp_all = getappdata(handles.figDSD,'XYZT');
  switch coord
   case 'XY'
    setappdata(handles.figDSD,'XYZT',[cp cp_all(3:4)]);
    %set(handles.textXY,'String',sprintf('(%d,%d)',cp));
    set(handles.dragX,'XData',[1 1]*cp(1));
    set(handles.dragY,'YData',[1 1]*cp(2));
    dsd_display(handles,{'XZ','TY'});
   case 'Z'
    setappdata(handles.figDSD,'XYZT',[cp_all(1:2) cp(2) cp_all(4)]);
    set(handles.textFrame,'String',sprintf('Frame %d',cp(2)));
    set(handles.dragZ,'YData',[1 1]*cp(2));
    dsd_display(handles,{'XY','TY'});
   case 'T'
    setappdata(handles.figDSD,'XYZT',[cp_all(1:3) cp(1)]);
    set(handles.textStack,'String',sprintf('Stack %d',cp(1)));
    set(handles.dragT,'XData',[1 1]*cp(1));
    dsd_display(handles,{'XY','XZ'});
  end

function dsd_updatenpix(handles,evt)
  if ishandle(handles)
    hfig = get_parent_fig(handles);
    handles = guidata(hfig);
  end
  selectXY = getappdata(handles.figDSD,'selectXY');
  yd = get([handles.dragZ1,handles.dragZ2],'YData');
  zrange = round([yd{1}(1) yd{2}(1)]);
  nvox = selectXY(3)*selectXY(4)*(diff(zrange)+1);
  set(handles.textNPix,'String',sprintf('%d voxels selected',nvox));
  
% Give starting positions for X,Y,Z,T based on the stack parameters
function dsd_initialize_positions(handles)
% Clear out any old objects
if ishandle(handles.selbox)
  delete(handles.selbox);
end
if ishandle(handles.dragZ1)
  delete(handles.dragZ1);
end
if ishandle(handles.dragZ2)
  delete(handles.dragZ2);
end
handles.selbox = [];
handles.dragZ1 = [];
handles.dragZ2 = [];
% Set data for new stack
stack = getappdata(handles.figDSD,'stack');
if isempty(stack)
  cp = [1 1 1 1];
else
  sz = stack.size;
  % Start in the middle of the very first frame acquired
  cp = [round(sz(1)/2) round(sz(2)/2) 1 1];
end
setappdata(handles.figDSD,'XYZT',cp);
set(handles.textFrame,'String',sprintf('Frame %d',cp(3)));
if ~isempty(stack)
  % Set axis limits
  set(handles.axesXY,'XLim',[1 sz(1)],'YLim',[1 sz(2)]);
  set(handles.axesXZ,'XLim',[1 sz(1)],'YLim',[1 sz(3)]);
  % Install draglines
  handles.dragZ1 = line([1 sz(2)],[2 2],'Color','r','Parent',handles.axesXZ);
  handles.dragZ2 = line([1 sz(2)],[1 1]*(sz(3)-1),'Color','r','Parent',handles.axesXZ);
  set([handles.dragZ1 handles.dragZ2],...
      'EraseMode','xor','LineWidth',3);
  %set([handles.dragZ handles.dragT],...
  %    'EraseMode','xor','HitTest','off');
  drag_line(handles.dragZ1,struct('onDragDone',@dsd_updatenpix));
  drag_line(handles.dragZ2,struct('onDragDone',@dsd_updatenpix));
  %drag_line(handles.dragT,struct('onDragDone',{@dsd_display,handles,{'XY','XZ'}}));
  %drag_line(handles.dragZ,struct('onDragDone',{@dsd_display,handles,{'XY','TY'}}));
  %drag_line(handles.dragT,struct('onDragDone',{@dsd_display,handles,{'XY','XZ'}}));
  guidata(handles.figDSD, handles);
end

function dsd_display(handles,panels)
stack = getappdata(handles.figDSD,'stack');
if isempty(stack)
  return
end
cp = getappdata(handles.figDSD,'XYZT');
for i = 1:length(panels)
  switch panels{i}
   case 'XY'
    im = stack(:,:,cp(3),cp(4));
    set(handles.imageXY,'CData',squeeze(im)');
   case 'XZ'
    %if get(handles.checkboxUpdateZ,'Value')
      im = stack(:,cp(2),:,cp(4));
      set(handles.imageXZ,'CData',squeeze(im)');
    %end
   otherwise
    error(['Panel ' panels{i} ' not recognized.']);
  end
  %hdrag = [handles.dragX handles.dragY handles.dragZ handles.dragT];
  %vstate = get(hdrag,'Visible');
  %set(hdrag,'Visible','off');
  %drawnow
  %set(hdrag,{'Visible'},vstate);
  %drawnow
end



% --- Executes on button press in btnPlayUp.
function btnPlayUp_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  cp = getappdata(handles.figDSD,'XYZT');
  stack = getappdata(handles.figDSD,'stack');
  sz = stack.size;
  setappdata(handles.figDSD,'isPlaying',1);
  while (cp(3) < sz(3) && getappdata(handles.figDSD,'isPlaying'))
    cp(3) = cp(3)+1;
    im = stack(:,:,cp(3),cp(4));
    set(handles.textFrame,'String',sprintf('Frame %d',cp(3)));
    %set(handles.dragZ,'YData',[1 1]*cp(3));
    set(handles.imageXY,'CData',squeeze(im)');
    drawnow
  end
  setappdata(handles.figDSD,'XYZT',cp);


% --- Executes on button press in btnStepUp.
function btnStepUp_Callback(hObject, eventdata, handles)
% hObject    handle to btnStepUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp = getappdata(handles.figDSD,'XYZT');
stack = getappdata(handles.figDSD,'stack');
cp(3) = cp(3)+1;
sz = stack.size;
if (cp(3) <= sz(3))
  setappdata(handles.figDSD,'XYZT',cp);
  set(handles.textFrame,'String',sprintf('Frame %d',cp(3)));
  %set(handles.dragZ,'YData',[1 1]*cp(3));
  dsd_display(handles,{'XY'});
end


% --- Executes on button press in btnStepDown.
function btnStepDown_Callback(hObject, eventdata, handles)
% hObject    handle to btnStepDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp = getappdata(handles.figDSD,'XYZT');
stack = getappdata(handles.figDSD,'stack');
cp(3) = cp(3)-1;
if (cp(3) > 0)
  setappdata(handles.figDSD,'XYZT',cp);
  set(handles.textFrame,'String',sprintf('Frame %d',cp(3)));
  %set(handles.dragZ,'YData',[1 1]*cp(3));
  dsd_display(handles,{'XY'});
end




% --- Executes on button press in btnPlayDown.
function btnPlayDown_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  cp = getappdata(handles.figDSD,'XYZT');
  stack = getappdata(handles.figDSD,'stack');
  setappdata(handles.figDSD,'isPlaying',1);
  while (cp(3) > 1 && getappdata(handles.figDSD,'isPlaying'))
    cp(3) = cp(3)-1;
    im = stack(:,:,cp(3),cp(4));
    set(handles.textFrame,'String',sprintf('Frame %d',cp(3)));
    %set(handles.dragZ,'YData',[1 1]*cp(3));
    set(handles.imageXY,'CData',squeeze(im)');
    drawnow
  end
  setappdata(handles.figDSD,'XYZT',cp);


% --- Executes on button press in btnPlayBack.
function btnPlayBack_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  cp = getappdata(handles.figDSD,'XYZT');
  stack = getappdata(handles.figDSD,'stack');
  setappdata(handles.figDSD,'isPlaying',1);
  while (cp(4) > 1 && getappdata(handles.figDSD,'isPlaying'))
    cp(4) = cp(4)-1;
    im = stack(:,:,cp(3),cp(4));
    set(handles.textStack,'String',sprintf('Stack %d',cp(4)));
    set(handles.dragT,'XData',[1 1]*cp(4));
    set(handles.imageXY,'CData',squeeze(im)');
    drawnow
  end
  setappdata(handles.figDSD,'XYZT',cp);
  dsd_display(handles,{'XZ'});


% --- Executes on button press in btnStepBack.
function btnStepBack_Callback(hObject, eventdata, handles)
% hObject    handle to btnStepBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp = getappdata(handles.figDSD,'XYZT');
stack = getappdata(handles.figDSD,'stack');
cp(4) = cp(4)-1;
if (cp(4) > 0)
  setappdata(handles.figDSD,'XYZT',cp);
  set(handles.textStack,'String',sprintf('Stack %d',cp(4)));
  set(handles.dragT,'XData',[1 1]*cp(4));
  dsd_display(handles,{'XY','XZ'});
end




% --- Executes on button press in btnStepForward.
function btnStepForward_Callback(hObject, eventdata, handles)
% hObject    handle to btnStepForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cp = getappdata(handles.figDSD,'XYZT');
stack = getappdata(handles.figDSD,'stack');
cp(4) = cp(4)+1;
sz = stack.size;
if (cp(4) <= sz(4))
  setappdata(handles.figDSD,'XYZT',cp);
  set(handles.textStack,'String',sprintf('Stack %d',cp(4)));
  set(handles.dragT,'XData',[1 1]*cp(4));
  dsd_display(handles,{'XY','XZ'});
end



% --- Executes on button press in btnPlayForward.
function btnPlayForward_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayForward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  cp = getappdata(handles.figDSD,'XYZT');
  stack = getappdata(handles.figDSD,'stack');
  sz = stack.size;
  setappdata(handles.figDSD,'isPlaying',1);
  while (cp(4) < sz(4) && getappdata(handles.figDSD,'isPlaying'))
    cp(4) = cp(4)+1;
    im = stack(:,:,cp(3),cp(4));
    im = squeeze(im)';
    set(handles.textStack,'String',sprintf('Stack %d',cp(4)));
    set(handles.dragT,'XData',[1 1]*cp(4));
    set(handles.imageXY,'CData',im);
    drawnow
  end
  setappdata(handles.figDSD,'XYZT',cp);
  dsd_display(handles,{'XZ'});
  

% --- Executes on button press in btnStop.
function btnStop_Callback(hObject, eventdata, handles)
% hObject    handle to btnStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  setappdata(handles.figDSD,'isPlaying',0);


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenu_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[basename,pathname] = uigetfile('*','Select header file');
smm = stackmm([pathname filesep basename]);
setappdata(handles.figDSD,'stack',smm);
setappdata(handles.figDSD,'basename',basename);
setappdata(handles.figDSD,'pathname',pathname);
% Initialize default locations
dsd_initialize_positions(handles);
% Display images data
dsd_display(handles,{'XY','XZ'});


% --------------------------------------------------------------------
function QuitMenu_Callback(hObject, eventdata, handles)
% hObject    handle to QuitMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figDSD);


% A helper function to draw the paging icons
function icons = swicons
  iconsize = [17 17];
  bkgrnd = get(0,'DefaultUicontrolBackgroundColor');
  blank = zeros(iconsize) + bkgrnd(1);
  
  right0 = blank;
  right1 = blank;
  index = mask(iconsize,0.8,1);
  right1(index) = 1;
  right0(index) = 0;
  % Long green arrow pointing right
  right(:,:,1) = right0;
  right(:,:,2) = right1;
  right(:,:,3) = right0;

  rightsm = blank;
  index = mask(iconsize,0.5,1);
  rightsm(index) = 0;
  rightsm = repmat(rightsm,[1 1 3]);  % Short black arrow pointing right

  left = right(:,end:-1:1,:);
  leftsm = rightsm(:,end:-1:1,:);
  
  icons.right = right;
  icons.rightsm = rightsm;
  icons.left = left;
  icons.leftsm = leftsm;
  icons.up = permute(left,[2 1 3]);
  icons.upsm = permute(leftsm,[2 1 3]);
  icons.down = permute(right,[2 1 3]);
  icons.downsm = permute(rightsm,[2 1 3]);
  icons.stop = repmat(0*blank,[1 1 3]);
  
% This function takes parameters from swicons
% to do the actual icon "drawing"
function index = mask(iconsize,width,direction)
  [row,col] = find(ones(iconsize));
  % Make them go from 0 to 1
  row = (row-1)/(iconsize(1)-1);
  col = (col-1)/(iconsize(2)-1);
  if (direction > 0)
    index = find(2*abs(row-0.5) <= (width-col)/width);
  else
    index = find(2*abs(row-0.5) + (1 - width)/width <= col/width);
  end
  return;


% --- Executes on button press in checkboxUpdateT.
function checkboxUpdateT_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUpdateT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUpdateT
if get(hObject,'Value')
  dsd_display(handles,{'TY'});
end

% --- Executes on button press in checkboxUpdateZ.
function checkboxUpdateZ_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxUpdateZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxUpdateZ
if get(hObject,'Value')
  dsd_display(handles,{'XZ'});
end


% --- Executes on button press in btnImageContrast.
function btnImageContrast_Callback(hObject, eventdata, handles)
% hObject    handle to btnImageContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  newclim = imrangegui(get(handles.imageXY,'CData'), ...
                       get(handles.axesXY,'CLim'), 0);
  set([handles.axesXY handles.axesXZ handles.axesTY],...
      'CLim',newclim);
  


% --- Executes on button press in checkboxShowCrossHairs.
function checkboxShowCrossHairs_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowCrossHairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% checkboxShowCrossHairs
  visible_value = {'off','on'};
  set([handles.dragX handles.dragY],'Visible',...
                    visible_value{get(hObject,'Value')+1});
  


% --- Executes on button press in btnCalc.
function btnCalc_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Pull out the selected region
selectXY = getappdata(handles.figDSD,'selectXY');
yd = get([handles.dragZ1 handles.dragZ2],'YData');
zrange = [yd{1}(1) yd{2}(1)];
stack = getappdata(handles.figDSD,'stack');
m = stack(selectXY(1):selectXY(1)+selectXY(3),...
          selectXY(2):selectXY(2)+selectXY(4),...
          zrange(1):zrange(2),1);
mz = stack2zcol(m);  % Convert to a zcolumn matrix
% If the user wants only a subset of the data, choose that
keep_frac = str2num(get(handles.editUsePercent,'String'))/100;
if (keep_frac <= 0 || keep_frac > 1)
  errordlg('Use percentage must be between 0 and 100');
end
stride = round(1/keep_frac);
mzs = mz(:,1:stride:end);
% Create the starting filter
filtsize = str2num(get(handles.editFiltSize,'String'));
halfsize = round((filtsize-1)/2);
x = -halfsize:halfsize;
f = exp(-x.^2/(2*(halfsize/2)^2));   % Initialize with a gaussian
f = f'/sum(f);
% Do the blind deconvolution
mzst = edgetaper(mzs,f);
nIterations = str2num(get(handles.editIterations,'String'));
psf = {f};
img = {mzst};
progops = progress_bar(struct('max',nIterations,...
  'progress',0,...
  'what','Calculating PSF...'));
for i = 1:nIterations
  progops.progress = i;
  progress_bar(progops);
  [img,psf] = deconvblind(img,psf,1,0.1);
end
figure; colormap(gray(256));
clim = get(handles.axesXY,'CLim');
subplot(1,2,1); imagesc(mzs,clim)
subplot(1,2,2); imagesc(img{2},[clim(1) clim(2)*1.5])
psf
figure
plot([f psf{2}])
title('PSF')
xlabel('z')
setappdata(handles.figDSD,'PSF',psf);
set(handles.btnSave,'Enable','on')


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
psf = getappdata(handles.figDSD,'PSF');
basename = getappdata(handles.figDSD,'basename');
pathname = getappdata(handles.figDSD,'pathname');
[filename,pathname] = uiputfile('*.psf','Choose an output file',...
                                [pathname filesep basename]);
nIterations = str2num(get(handles.editIterations,'String'));
save([pathname filesep filename],'psf','nIterations');
                    


function editUsePercent_Callback(hObject, eventdata, handles)
% hObject    handle to editUsePercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editUsePercent as text
%        str2double(get(hObject,'String')) returns contents of editUsePercent as a double


% --- Executes during object creation, after setting all properties.
function editUsePercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editUsePercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFiltSize_Callback(hObject, eventdata, handles)
% hObject    handle to editFiltSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFiltSize as text
%        str2double(get(hObject,'String')) returns contents of editFiltSize as a double


% --- Executes during object creation, after setting all properties.
function editFiltSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFiltSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editIterations_Callback(hObject, eventdata, handles)
% hObject    handle to editIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editIterations as text
%        str2double(get(hObject,'String')) returns contents of editIterations as a double


% --- Executes during object creation, after setting all properties.
function editIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


