function varargout = register_movie_gui(varargin)
% REGISTER_MOVIE_GUI: register grayscale volumes over time
%
% This GUI allows you to interactively register "4D" recordings (volumes +
% time). It works on grayscale images.
%
% Syntax:
%   register_movie_gui
% This will prompt you to pick a file to analyze
%
%   register_movie_gui(filename)
%   register_movie_gui(smm)
% With this syntax you specify the file, or stackmm object, directly.
%
%   register_movie_gui(...,options)
% allows you to control the behavior more explicitly, this may be the 
% output from find_bad_frames:
% isbad: an n_frames-by-n_stacks logical array true whenever a given frame
%   is bad. If you supply this field, then the default "stacknums" will
%   exclude any stack that has a bad frame. See find_bad_frames (or
%   register_dfof_validate_stack) for a function that can calculate "isbad"
%   for you.
% shift: supplies the rigid shift if previously determined.
% stacknums: the initial set of stacks that are selected for possible
%   registration. The default value is the set of all stacks immediately
%   preceding a stimulus. The GUI allows you to add, delete, or shift
%   (using the left/right arrow keys) the selected stacks.
% base_stacknum: the reference stack used as the "fixed" image in all
%   registrations. By default, this is the middle item in "stacknums."
%
% See also: register_movie_dfof_gui, find_bad_frames.

% Copyright 2010-2011 by Timothy E. Holy

% REGISTER_MOVIE_GUI M-file for register_movie_gui.fig
%      REGISTER_MOVIE_GUI, by itself, creates a new REGISTER_MOVIE_GUI or raises the existing
%      singleton*.
%
%      H = REGISTER_MOVIE_GUI returns the handle to a new REGISTER_MOVIE_GUI or the handle to
%      the existing singleton*.
%
%      REGISTER_MOVIE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_MOVIE_GUI.M with the given input arguments.
%
%      REGISTER_MOVIE_GUI('Property','Value',...) creates a new REGISTER_MOVIE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before register_movie_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_movie_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_movie_gui

% Last Modified by GUIDE v2.5 24-Jan-2011 06:03:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @register_movie_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @register_movie_gui_OutputFcn, ...
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


% --- Executes just before register_movie_gui is made visible.
function register_movie_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_movie_gui (see VARARGIN)

options = struct;
if isempty(varargin)
  filename = uigetfile('*','Pick file to register');
  if isequal(filename,0)
    delete(handles.figMain)
    return
  end
  smm = stackmm(filename);
else
  if ischar(varargin{1}) || iscell(varargin{1})
    filename = varargin{1};
    smm = stackmm(filename);
  elseif isa(varargin{1},'stackmm')
    smm = varargin{1};
  end
end
if (length(varargin) > 1)
  options = varargin{2};
end

sz = smm.size;
header = smm.header;
if isfield(header,'stim_lookup') && ~isempty(header.stim_lookup)
  if any(header.stim_lookup == 0)
    prestim = find(header.stim_lookup(1:end-1) == 0 & header.stim_lookup(2:end) > 0);
  else
    prestim = find(diff(header.stim_lookup));
  end
else
  if ~isfield(options,'stack_skip')
    error('When there is no stimulus information, you must supply a "stack_skip" field in options');
  end
  prestim = 1:options.stack_skip:sz(4);
end
prestim = prestim(:)';
if isfield(options,'isbad')
  badstacks = find(any(options.isbad,1));
  prestim = setdiff(prestim,badstacks);
end
options = default(options,'stacknums',prestim);
options = default(options,'base_stacknum',options.stacknums(round(length(options.stacknums)/2)));
% Remove the base_stacknum from stacknums
indx = find(options.stacknums == options.base_stacknum);
if ~isempty(indx)
  options.stacknums(indx) = [];
end

setappdata(handles.figMain,'base_stacknum',options.base_stacknum)
setappdata(handles.figMain,'stacknums',options.stacknums)
selected = false(1,length(options.stacknums));
setappdata(handles.figMain,'selected',selected);
setappdata(handles.figMain,'smm',smm);
setappdata(handles.figMain,'sz',sz);
setappdata(handles.figMain,'options',options);
funch = register_gui_utilities;
setappdata(handles.figMain,'funch',funch);

% Pre-allocate results
if isfield(options,'err')
  err = nanmean(options.err,1)';
else
  err = NaN(sz(4),1);
end
erru = NaN(sz(4),1);
if isfield(options,'shift')
  shift = options.shift';
else
  shift = NaN(sz(4),3);
end
udata = cell(sz(4),1);
setappdata(handles.figMain,'err',err)
setappdata(handles.figMain,'erru',erru)
setappdata(handles.figMain,'shift',shift)
setappdata(handles.figMain,'udata',udata)

% Prepare graphs
set([handles.axesSel handles.axesErr handles.axesShift handles.axesUDiff],'XLim',[-1 sz(end)+1],'TickDir','out')
set([handles.axesErr handles.axesShift handles.axesUDiff],'XTick',[])
set(handles.axesSel,'YTick',[])
ylabel(handles.axesErr,'Error')
ylabel(handles.axesShift,'Shift')
ylabel(handles.axesUDiff,'u-difference')
% Plot bars corresponding to stimulation
stim = header.stim_lookup;
col = repmat(0.8,max(stim),3);  % use gray for all bars
h = plot_stimulus_bars(handles.axesSel,stim,col,[0 1]);
set(h,'HitTest','off')
% Make start and stop lines
h = line([1 1],[0 1],'Parent',handles.axesSel,'Color','k','LineWidth',2);
drag_line(h)
handles.dragStart = h;
h = line(sz(end)+[0.5 0.5],[0 1],'Parent',handles.axesSel,'Color','k','LineWidth',2);
drag_line(h)
handles.dragStop = h;
% Set up buttondownfcn for de-selecting
set(handles.axesSel,'ButtonDownFcn',@rmg_select);
% Set a default frame
set(handles.editFrame,'String',num2str(round(sz(3)/2)));

% Set up defaults for parameters
maskops = struct('margin_frac',0.1,'margin_pixel',2);
setappdata(handles.figMain,'maskops',maskops);

% Set parameters for warping
stk = single(smm(:,:,:,options.base_stacknum));
% mask = true(size(stk));  % a "fake" mask
pyramid = array_restrict_schedule(size(stk),header);
% rmg_params = register_multigrid_options(stk,mask,struct('pixel_spacing',header.pixel_spacing));
warpops.gridsize = cat(1,pyramid.sz);
warpops.mspixval = mean(stk(:).^2);
warpops.min_pixels = 7;
warpops.min_g_pixels = 2;
warpops.gap_data = 3;
% warpops.min_pixels = rmg_params.min_pixels;
% warpops.min_g_pixels = rmg_params.min_g_pixels;
% warpops.gap_data = rmg_params.gap_data;
setappdata(handles.figMain,'warpops',warpops);
pcoptions = register_phasecorr_initialize(stk,struct('pyramid',pyramid));
h = fspecial('gaussian',5,1);
pcoptions.smooth = @(im)imfilter(im,h);
setappdata(handles.figMain,'pcoptions',pcoptions);

% Choose default command line output for register_movie_gui
handles.output = hObject;

% Indicate the "parent" GUI so register_gui_utilities knows what to do
handles.caller = 'register_movie_gui';

% Update handles structure
guidata(hObject, handles);

% Do initial plotting
rmg_plot_stacknums(handles)
rmg_plot_results(handles)

% UIWAIT makes register_movie_gui wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = register_movie_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
  varargout{1} = handles.output;
end



function rmg_plot_stacknums(handles)
  stacknums = getappdata(handles.figMain,'stacknums');
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  h = findobj(handles.axesSel,'type','line');
  h = setdiff(h,[handles.dragStart handles.dragStop]);
  if ~isempty(h)
    delete(h)
  end
  hline = line([stacknums; stacknums],repmat([0;1],1,length(stacknums)),'Color','b','Parent',handles.axesSel,'ButtonDownFcn',@rmg_select);
  selected = getappdata(handles.figMain,'selected');
  offon = {'off';'on'};
  set(hline,{'Selected'},offon(1+selected(:)));
  line([base_stacknum base_stacknum],[0 1],'Color','r','Parent',handles.axesSel);
  if (sum(selected) == 1)
    set(handles.btnPlayStack,'Enable','on');
  else
    set(handles.btnPlayStack,'Enable','off');
  end
  btns = [handles.btnPlayMovie handles.btnRegisterRigid handles.btnRegisterWarp];
  if (sum(selected) > 0)
    set(btns,'Enable','on');
  else
    set(btns,'Enable','off');
  end
     
function rmg_plot_results(handles)
  err = getappdata(handles.figMain,'err');
  erru = getappdata(handles.figMain,'erru');
  shift = getappdata(handles.figMain,'shift');
  smm = getappdata(handles.figMain,'smm');
  header = smm.header;
  cla(handles.axesErr)
  line(1:length(err),[err erru],'Parent',handles.axesErr,'Marker','x')
  cla(handles.axesShift)
  line(1:length(err),bsxfun(@times,shift,header.pixel_spacing),'Parent',handles.axesShift,'Marker','x')
  funch = getappdata(handles.figMain,'funch');
  funch.plot_udiff_gui(handles);
%   btns = handles.btnMask;
%   if any(~isnan(shift(:)))
%     set(btns,'Enable','on')
%   else
%     set(btns,'Enable','off')
%   end

function rmg_select(src,evt)
  hfig = get_parent_fig(src);
  seltype = get(hfig,'SelectionType');
  handles = guidata(hfig);
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  isline = false;
  if strcmp(get(src,'type'),'line')
    isline = true;
  end
  if isline
    x = get(src,'XData');
    x = x(1);
    indx = find(stacknums == x);
  end
  switch seltype
    case 'normal'
      selected = false(size(selected));
      if isline
        selected(indx) = true;
      end
    case {'extend','alt'}
      if isline
        selected(indx) = ~selected(indx);
      end
  end
  setappdata(handles.figMain,'selected',selected);
  rmg_plot_stacknums(handles)
  
function u = rmg_match_u(u,pcoptions,desired_level)
  if isempty(u) || isempty(u{1})
    u = cell(1,3);
    for dimIndex = 1:3
      u{dimIndex} = zeros(pcoptions.pyramid(desired_level).sz,pcoptions.class);
    end
  else
    % Match to a pyramid level
    sz = cat(1,pcoptions.pyramid.sz);
    usz = size(u{1});
    if all(usz == 1)
      for dimIndex = 1:3
        u{dimIndex} = repmat(u{dimIndex},pcoptions.pyramid(desired_level).sz);
      end
      return
    else
      eq = sz == repmat(usz,size(sz,1),1);
      current_level = find(all(eq,2));
      if isempty(current_level)
        error('Didn''t match the level');
      end
    end
    while (current_level > desired_level)
      % Prolong u
      current_level = current_level-1;
      for dimIndex = 1:3
        u{dimIndex} = array_prolong(u{dimIndex},pcoptions.pyramid(current_level).sz);
      end
    end
    while (current_level < desired_level)
      % Restrict u
      current_level = current_level+1;
      for dimIndex = 1:3
        u{dimIndex} = array_restrict(u{dimIndex},pcoptions.pyramid(current_level).restrict);
      end
    end
  end

function umg = rmg_upc2mg(u,imsz)
  % Convert phasecorr u representation to multigrid
  umg = cat(4,u{:});
  szu = size(umg);
  umg = bsxfun(@rdivide,umg,reshape(imsz./szu(1:3),[1 1 1 3]));

function upc = rmg_umg2pc(umg,imsz)
  % Convert multigrid u representation to phasecorr
  szu = size(umg);
  umg = bsxfun(@times,umg,reshape(imsz./szu(1:3),[1 1 1 3]));
  upc = mat2cell(umg,szu(1),szu(2),szu(3),[1 1 1]);
  upc = upc(:)';
  
function shift = rmg_getshift(handles,stacknum)
  shift = getappdata(handles.figMain,'shift');
  shift = shift(stacknum,:);
  if any(isnan(shift))
    warndlg('Shift not defined for this stack')
  end
  shift(isnan(shift)) = 0;
  
  
% function rmg_create_mask(handles)
%   % Collect information about maximum shift across selected time range, so
%   % we can make a mask that excludes the pixels that will go over the edge
%   sz = getappdata(handles.figMain,'sz');
%   shift = getappdata(handles.figMain,'shift');
%   nanflag = isnan(shift(:,1));
%   keepIndex = find(~nanflag);
%   trange = rmg_get_trange(handles);
%   keepFlag = keepIndex >= trange(1) & keepIndex <= trange(2);
%   keepIndex = keepIndex(keepFlag);
%   shift = shift(keepIndex,:);
%   % Get mask parameters from user
%   maskops = getappdata(handles.figMain,'maskops');
%   fn = fieldnames(maskops);
%   defans = struct2cell(maskops);
%   for i = 1:length(defans)
%     defans{i} = num2str(defans{i});
%   end
%   answer = inputdlg(fieldnames(maskops),'Enter mask-creation options',1,defans);
%   if isempty(answer)
%     return
%   end
%   for i = 1:length(answer)
%     answer{i} = str2double(answer{i});
%   end
%   maskops = cell2struct(answer,fn,1);
%   setappdata(handles.figMain,'maskops',maskops);
%   % Create the mask
%   mask = register_shift2mask(sz,shift,maskops);
%   msgbox(sprintf('Masking out %d%% of the pixels',round(100*sum(mask(:) == 0)/numel(mask))));
%   % Create rmg_params
%   smm = getappdata(handles.figMain,'smm');
%   header = smm.header;
%   base_stacknum = getappdata(handles.figMain,'base_stacknum');
%   base_stack = smm(:,:,:,base_stacknum);
%   rmg_params = register_multigrid_options(double(base_stack),mask,struct('pixel_spacing',header.pixel_spacing));
%   setappdata(handles.figMain,'rmg_params',rmg_params);

% --- Executes on button press in btnRegisterRigid.
function btnRegisterRigid_Callback(hObject, eventdata, handles)
% hObject    handle to btnRegisterRigid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Get the selected stack numbers
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  % Get the base stack
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  smm = getappdata(handles.figMain,'smm');
  base_stack = single(smm(:,:,:,base_stacknum));
  % Load information about the error & shift
  err = getappdata(handles.figMain,'err');
  shift = getappdata(handles.figMain,'shift');
  % Perform the registration
  pgoptions = struct('max',length(stacknums));
  tic
  for i = 1:length(stacknums)
    thisstack = stacknums(i);
    stk = single(smm(:,:,:,thisstack));
    [stkr,thisshift] = register_rigid(base_stack,stk);
    shift(thisstack,:) = thisshift;
    imdiff = base_stack - stkr;
%     err(thisstack) = nanmean(abs(imdiff(:)))/nanmean(stkr(:));
    err(thisstack) = nanmean(abs(imdiff(:).^2));
    if (toc > 3 || i == length(stacknums))
      pgoptions.progress = i;
      pgoptions = progress_bar(pgoptions);
      tic
    end
  end
  % Save the results
  setappdata(handles.figMain,'err',err)
  setappdata(handles.figMain,'shift',shift)
  rmg_plot_results(handles)


% --- Executes on button press in btnToggleSelect.
function btnToggleSelect_Callback(hObject, eventdata, handles)
% hObject    handle to btnToggleSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  selected = getappdata(handles.figMain,'selected');
  stacknums = getappdata(handles.figMain,'stacknums');
  funch = getappdata(handles.figMain,'funch');
  trange = funch.get_trange(handles);
  inrange = stacknums >= trange(1) & stacknums <= trange(2);
  selected(inrange) = ~selected(inrange);
  setappdata(handles.figMain,'selected',selected);
  rmg_plot_stacknums(handles)
    

% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Get the data we want to save
  err = getappdata(handles.figMain,'err');
  erru = getappdata(handles.figMain,'erru');
  shift = getappdata(handles.figMain,'shift');
  udata = getappdata(handles.figMain,'udata');
  stacknums = getappdata(handles.figMain,'stacknums');
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  options = getappdata(handles.figMain,'options');
  smm = getappdata(handles.figMain,'smm');
  filename = smm.filename;
  header = smm.header;
  pixel_spacing = header.pixel_spacing;
  funch = getappdata(handles.figMain,'funch');
  trange = funch.get_trange(handles);
  [savefilename,pathname] = uiputfile;
  if ~isequal(filename,0)
    save([pathname savefilename],'filename','pixel_spacing','stacknums','base_stacknum','err','erru','shift','udata','trange','options', '-v7.3');
  end
    
% --- Executes on button press in btnQuit.
function btnQuit_Callback(hObject, eventdata, handles)
% hObject    handle to btnQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnPlayStack.
function btnPlayStack_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  funch = getappdata(handles.figMain,'funch');
  funch.play_stack(handles)
  return
  
  smm = getappdata(handles.figMain,'smm');
  % Get the base stack
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  base_stack = single(smm(:,:,:,base_stacknum));
  % Get the selected stack
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  stk = single(smm(:,:,:,stacknums));
    
  % Register the selected stack
  playType = get(handles.popupmenuPlayType,'Value');
  switch playType
    case 1
      stkr = stk;
    case 2
      shift = rmg_getshift(handles,stacknums);
      stkr = image_shift(stk,shift);
    case 3
      shift = rmg_getshift(handles,stacknums);
      udata = getappdata(handles.figMain,'udata');
      u = udata{stacknums};
      pcoptions = getappdata(handles.figMain,'pcoptions');
      pcoptions.shift = shift;
      stkr = register_phasecorr_warp(u,stk,pcoptions);
  end
  % Generate the movie
  sz = size(stkr);
  rgb = zeros([sz(1:2) 3 sz(3)]);
  rgb(:,:,1,:) = base_stack;
  rgb(:,:,2,:) = stkr;
  rgb = double(rgb) / double(max(rgb(:)));
  mplay(rgb)

% --- Executes on button press in btnPlayMovie.
function btnPlayMovie_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  funch = getappdata(handles.figMain,'funch');
  funch.play_movie(handles)
  return
  
  frame = str2double(get(handles.editFrame,'String'));
  smm = getappdata(handles.figMain,'smm');
%   % Get the base frame
%   base_stacknum = getappdata(handles.figMain,'base_stacknum');
%   base_frame = single(smm(:,:,frame,base_stacknum));
  % Get the selected stack number(s)
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  if isempty(stacknums)
    return
  end
  % Get the shifts or deformation
  playType = get(handles.popupmenuPlayType,'Value');
  if (playType > 1)
    shift = getappdata(handles.figMain,'shift');
    shift = shift(stacknums,:);
    if any(isnan(shift(:)))
      errordlg('The shift (rigid registration) is not defined for all selected stacks');
      return
    end
    if (playType == 3)
      udata = getappdata(handles.figMain,'udata');
      udata = udata(stacknums);
      if any(cellfun(@isempty,udata))
        errordlg('The deformation (warp registration) is not defined for all selected stacks');
        return
      end
      pcoptions = getappdata(handles.figMain,'pcoptions');
    end
  end
  % Iteratively get stacks, register them, and select the frame
  sz = smm.size;
  mov = zeros([sz(1:2) length(stacknums)],smm.type);
  pgoptions = struct('max',length(stacknums));  % progress bar
  tic
  for i = 1:length(stacknums)
    thisstk = stacknums(i);
    if playType == 2
      stk = smm(:,:,:,thisstk);
      % Register the selected stack
      stkr = image_shift(stk,shift(i,:));
      mov(:,:,i) = stkr(:,:,frame);
      if (toc > 3 || i == length(stacknums))
        pgoptions.progress = i;
        pgoptions = progress_bar(pgoptions);
        tic
      end
    elseif playType == 3
      stk = smm(:,:,:,thisstk);
      % Register the selected stack
      u = udata{i};
      pcoptions.shift = shift(i,:);
      stkr = register_phasecorr_warp(u,single(stk),pcoptions);
      mov(:,:,i) = stkr(:,:,frame);
      if (toc > 3 || i == length(stacknums))
        pgoptions.progress = i;
        pgoptions = progress_bar(pgoptions);
        tic
      end
    else
      mov(:,:,i) = smm(:,:,frame,thisstk);
    end
  end
  % Play the movie
  movmax = max(mov(:));
  mov = uint8(mov*(255/double(movmax)));
  mplay(mov)

% --- Executes on button press in btnAdd.
function btnAdd_Callback(hObject, eventdata, handles)
% hObject    handle to btnAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  set([handles.axesErr handles.axesShift],'Color',[0.5 0.5 0.5]);
  btns = [handles.btnRegisterRigid handles.btnRegisterWarp handles.btnToggleSelect handles.btnPlayStack ...
    handles.btnPlayMovie handles.btnSave handles.btnLoad];
  state = get(btns,'Enable');
  set(btns,'Enable','off');
  hc = findobj(handles.axesSel,'type','line');
  set(hc,'HitTest','off','Selected','off');
  w = waitforbuttonpress;
  if (w == 0 && gca == handles.axesSel)
    cp = get(handles.axesSel,'CurrentPoint');
    stacknums = getappdata(handles.figMain,'stacknums');
    selected = getappdata(handles.figMain,'selected');
    [stacknums,sortOrder] = sort([stacknums round(cp(1))]);
    selected = [selected false];
    selected = selected(sortOrder);
    setappdata(handles.figMain,'stacknums',stacknums);
    setappdata(handles.figMain,'selected',selected);
  end
  set([handles.axesErr handles.axesShift],'Color',[1 1 1]);
  set(btns,{'Enable'},state);
  rmg_plot_stacknums(handles)
  

% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [filename,pathname] = uigetfile;
  if isequal(filename,0)
    return
  end
  if ~isequal(filename,0)
    s = load([pathname filename]);
  end
  setappdata(handles.figMain,'err',s.err);
  setappdata(handles.figMain,'shift',s.shift);
  setappdata(handles.figMain,'stacknums',s.stacknums);
  setappdata(handles.figMain,'base_stacknum',s.base_stacknum);
  if isfield(s,'trange')
    set(handles.dragStart,'XData',s.trange(1)*[1 1]);
    set(handles.dragStop,'XData',s.trange(2)*[1 1]);
  end
  if isfield(s,'erru')
    setappdata(handles.figMain,'erru',s.erru);
  end
  if isfield(s,'udata')
    setappdata(handles.figMain,'udata',s.udata);
  end
  selected = false(size(s.stacknums));
  setappdata(handles.figMain,'selected',selected);
  rmg_plot_stacknums(handles)
  rmg_plot_results(handles)

function editFrame_Callback(hObject, eventdata, handles)
% hObject    handle to editFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrame as text
%        str2double(get(hObject,'String')) returns contents of editFrame as a double


% --- Executes during object creation, after setting all properties.
function editFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figMain and none of its controls.
function figMain_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

  switch eventdata.Key
    case {'backspace','delete','leftarrow','rightarrow'}
      % Note: these next lines assume you're going to modify some aspect of the
      % stack numbers.
      stacknums = getappdata(handles.figMain,'stacknums');
      selected = getappdata(handles.figMain,'selected');
      err = getappdata(handles.figMain,'err');
      erru = getappdata(handles.figMain,'erru');
      shift = getappdata(handles.figMain,'shift');
      udata = getappdata(handles.figMain,'udata');
      selsn = stacknums(selected);
      err(selsn) = nan;
      erru(selsn) = nan;
      shift(selsn,:) = nan;
      udata(selsn) = {[]};
      switch eventdata.Key
        case {'backspace','delete'}
          stacknums(selected) = [];
          selected = false(size(stacknums));
        case 'leftarrow'
          stacknums(selected) = selsn - 1;
          stacknums(stacknums < 1) = 1;
        case 'rightarrow'
          sz = getappdata(handles.figMain,'sz');
          stacknums(selected) = selsn + 1;
          stacknums(stacknums > sz(4)) = sz(4);
        otherwise
          %display(eventdata.Key)
      end
      setappdata(handles.figMain,'stacknums',stacknums);
      setappdata(handles.figMain,'selected',selected);
      setappdata(handles.figMain,'err',err);
      setappdata(handles.figMain,'erru',erru);
      setappdata(handles.figMain,'shift',shift);
      setappdata(handles.figMain,'udata',udata);
      rmg_plot_stacknums(handles)
      rmg_plot_results(handles)
  end


% --- Executes on selection change in popupmenuPlayType.
function popupmenuPlayType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuPlayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuPlayType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuPlayType


% --- Executes during object creation, after setting all properties.
function popupmenuPlayType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuPlayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnRegisterWarp.
function btnRegisterWarp_Callback(hObject, eventdata, handles)
% hObject    handle to btnRegisterWarp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

  % Check to see if all selected stacks have been rigid-registered
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  shift = getappdata(handles.figMain,'shift');
  shift = shift(stacknums,:);
  if any(isnan(shift(:)))
    errordlg('Not all stacks have been rigid-registered, you need to do that first');
    return
  end
  
  funch = getappdata(handles.figMain,'funch');
  
  % Set the options via the dialog
  ok = funch.set_warpops(handles);
  if ~ok
    return
  end
  warpops = getappdata(handles.figMain,'warpops');
  
  drawnow
  
  % Choose the stack order
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  [~,sortOrder] = sort(abs(stacknums - base_stacknum)); % order by distance from base
  stacknums = stacknums(sortOrder);
  shift = shift(sortOrder,:);

  % Load other needed variables
  pcoptions = getappdata(handles.figMain,'pcoptions');
  udata = getappdata(handles.figMain,'udata');
  erru = getappdata(handles.figMain,'erru');
  
  % Prepare the fixed image data
  smm = getappdata(handles.figMain,'smm');
  if strcmp(warpops.method,'Multigrid least-squares')
    % If the user is doing multigrid, check to see if we have rmg_params
    % and create it if not
    if (warpops.createmask || ~isappdata(handles.figMain,'rmg_params'))
      funch.create_mask(handles)
    end
    rmg_params = getappdata(handles.figMain,'rmg_params');
    rmg_params.gap_data = warpops.gap_data;
    rmg_params.wcycle = warpops.wcycle;
    if any(warpops.presmooth > 0)
      rmg_params_smooth = register_multigrid_options_compress(rmg_params);
      base_stack = rmg_params.image_grid(1).imFixed;
      rmg_params_smooth = register_multigrid_options_expand(rmg_params_smooth,imfilter_gaussian_mex(base_stack,warpops.presmooth));
    end
  else
    % Phasecorr
    fixed = single(smm(:,:,:,base_stacknum));
    % phasecorr does not benefit from pre-smoothing, so don't worry about
    % supporting it
  end

  pgoptions = struct('max',length(stacknums),'progress',0);
  pgoptions = progress_bar(pgoptions);
  tic
  % Loop over selected stacks
  % Phasecorr should utilize provided badstack info
   hfig = get(hObject,'Parent');
   hfig_options = getappdata(hfig,'options');
   if isfield(hfig_options, 'isbad');
       isbad = hfig_options.isbad;
   else isbad = [];
   end
   if ~isempty(isbad) && gt(max(stacknums),size(isbad,2))
       warning('Error in size of isbad');
   end
   
  for i = 1:length(stacknums)
    thisstacknum = stacknums(i);
    % Determine the initial guess
    u = funch.u_initial_guess(udata,thisstacknum,base_stacknum,pcoptions,warpops);
    % Do the registration
    stk = smm(:,:,:,thisstacknum);
    switch warpops.method
      case 'Multigrid least-squares'
        rmg_params.shift = shift(i,:);
        stk = double(stk);
        umg = double(rmg_upc2mg(u,pcoptions.sz_spatial));
        if any(warpops.presmooth > 0)
          stksm = imfilter_gaussian_mex(stk,warpops.presmooth);
          [umg,thiserr] = register_multigrid_vcycle(umg,stksm,warpops.lambda,rmg_params_smooth);
          if warpops.refit
            [umg,thiserr] = register_multigrid_vcycle(umg,stk,warpops.lambda,rmg_params);
          end
        else
          [umg,thiserr] = register_multigrid_vcycle(umg,stk,warpops.lambda,rmg_params);
        end
        u = rmg_umg2pc(single(umg),pcoptions.sz_spatial);
      case 'Phase correlation'
        pcoptions.shift = shift(i,:);
        stk = single(stk);
        error_flag = false; % Mark for rewarping with nan values
        if ~isempty(isbad)  
            if any(isbad(:,thisstacknum))
                stk(:,:,isbad(:,thisstacknum)) = 0;
                error_flag = true;
            end
        end
        stk(isnan(stk)) = 0;
        pcoptions.lambda = warpops.lambda;
        if isempty(u) || isempty(u{1})
          multires_params = struct('desired_ulevel',warpops.gap_data+1,...
            'display',true);
          [u,~,thiserr] = register_phasecorr_multires(fixed,stk,pcoptions,multires_params);
          thiserr = thiserr(end);
        else
          u = register_phasecorr_improve(u,fixed,stk,pcoptions);
          stkw = register_phasecorr_warp(u,stk,pcoptions);
          thiserr = nanmean((fixed(:)-stkw(:)).^2);
        end
        if error_flag
            stk(:,:,isbad(:,thisstacknum)) = nan;
            stkw = register_phasecorr_warp(u,stk,pcoptions);
            thiserr = nanmean((fixed(:)-stkw(:)).^2);
        end
    end
    if any(isnan(u{1}(:)))
      error('NaN u');
    end
    udata{thisstacknum} = u;
    erru(thisstacknum) = thiserr;
    if (toc > 3 || thisstacknum == stacknums(end))
      setappdata(handles.figMain,'erru',erru);
      setappdata(handles.figMain,'udata',udata);
      rmg_plot_results(handles)
      pgoptions.progress = find(thisstacknum == stacknums);
      pgoptions = progress_bar(pgoptions);
      tic
    end
  end
  setappdata(handles.figMain,'udata',udata);
  setappdata(handles.figMain,'erru',erru);
%   set(btns,{'Enable'},state);
  rmg_plot_results(handles)
          
        
  

% --- Executes on button press in btnMask.
function btnMask_Callback(hObject, eventdata, handles)
% hObject    handle to btnMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  rmg_create_mask(handles)
  


% --- Executes on button press in btnClearU.
function btnClearU_Callback(hObject, eventdata, handles)
% hObject    handle to btnClearU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  stacknums = getappdata(handles.figMain,'stacknums');
  selected = getappdata(handles.figMain,'selected');
  stacknums = stacknums(selected);
  udata = getappdata(handles.figMain,'udata');
  erru = getappdata(handles.figMain,'erru');
  for i = stacknums
    udata{i} = [];
  end
  erru(stacknums) = NaN;
  setappdata(handles.figMain,'udata',udata);
  setappdata(handles.figMain,'erru',erru);
  rmg_plot_results(handles);


% --- Executes on button press in btnSmoothU.
function btnSmoothU_Callback(hObject, eventdata, handles)
% hObject    handle to btnSmoothU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  udata = getappdata(handles.figMain,'udata');
  err = getappdata(handles.figMain,'erru');
  pcoptions = getappdata(handles.figMain,'pcoptions');
  smm = getappdata(handles.figMain,'smm');
  header = smm.header;
  pixel_spacing = header.pixel_spacing;
  funch = getappdata(handles.figMain,'funch');
  errfunc = funch.calculate_mismatch;
  thisfunc = @(stacknum,u) errfunc(handles,stacknum,'Registered, warped',u);
  [udatanew,errnew] = register_movie_smoothu_gui(udata,err,thisfunc,pcoptions,pixel_spacing);
  if ~isempty(udatanew)
    setappdata(handles.figMain,'udata',udatanew);
    setappdata(handles.figMain,'erru',errnew);
    rmg_plot_results(handles)
  end
  