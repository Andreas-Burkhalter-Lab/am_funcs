function varargout = register_movie_dfof_gui(varargin)
% register_movie_dfof_gui: register a movie using stimulus responses
% Syntax:
%   register_movie_dfof_gui(smm,params)
% where
%   smm is a stackmm object
%   params is a structure that is the output of
%     register_dfof_validate_stacks. 
% If the DeltaF/F images look noisy to you, consider adding a field
% "dfof_smooth" to params. The syntax of this field is either a vector,
% such as [1 1 0] specifying the width along each axis for gaussian
% smoothing, or a function handle such as
%   h = fspecial('gaussian',5,1);
%   params.dfof_smooth = @(im) imfilter(im,h);
% For small filter sizes, the latter will be faster.
%
%   register_movie_dfof_gui(smm,filename)
% launches the GUI using a file that has been saved using the "Save"
% button.
%
% See also: register_dfof_validate_stacks, register_movie_gui.

% Copyright 2011 by Timothy E. Holy

% REGISTER_MOVIE_DFOF_GUI M-file for register_movie_dfof_gui.fig
%      REGISTER_MOVIE_DFOF_GUI, by itself, creates a new REGISTER_MOVIE_DFOF_GUI or raises the existing
%      singleton*.
%
%      H = REGISTER_MOVIE_DFOF_GUI returns the handle to a new REGISTER_MOVIE_DFOF_GUI or the handle to
%      the existing singleton*.
%
%      REGISTER_MOVIE_DFOF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_MOVIE_DFOF_GUI.M with the given input arguments.
%
%      REGISTER_MOVIE_DFOF_GUI('Property','Value',...) creates a new REGISTER_MOVIE_DFOF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before register_movie_dfof_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_movie_dfof_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_movie_dfof_gui

% Last Modified by GUIDE v2.5 23-Jan-2011 10:55:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @register_movie_dfof_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @register_movie_dfof_gui_OutputFcn, ...
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


% --- Executes just before register_movie_dfof_gui is made visible.
function register_movie_dfof_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_movie_dfof_gui (see VARARGIN)

smm = varargin{1};
options = varargin{2};
loadsaved = false;
if ischar(options)
  loadsaved = true;
  savefile = options;
  stmp = load(savefile,'options');
  options = stmp.options;
end
[~,basename] = fileparts(smm.filename);
options = default(options,'dfof_dir',['/tmp/register-' getenv('USER') '/' basename]);
if ~exist(options.dfof_dir,'dir')
  mkdir(options.dfof_dir)
else
  answer = questdlg(['Directory ' options.dfof_dir ' already exists. Re-use any existing files?']);
  switch answer
    case {'','Cancel'}
      delete(handles.figMain)
      return
    case 'No'
      rmdir(options.dfof_dir,'s');
      mkdir(options.dfof_dir);
  end
end

sz = smm.size;
header = smm.header;
funch = register_gui_utilities;
setappdata(handles.figMain,'smm',smm);
setappdata(handles.figMain,'sz',sz);
setappdata(handles.figMain,'options',options);
setappdata(handles.figMain,'base_stacknum',options.base_stacknum);
setappdata(handles.figMain,'funch',funch);

%% Find the trials with enough good data to proceed
[~,onset] = funch.choose_good_onsets(smm,options,options.prestim,options.peristim);
onset = onset(options.include);
stacknums = sort(cat(1,onset{:}))';
% % Remove the base_stacknum from stacknums
% indx = find(stacknums == options.base_stacknum);
% if ~isempty(indx)
%   stacknums(indx) = [];
% end
setappdata(handles.figMain,'stacknums',stacknums)
setappdata(handles.figMain,'onset',onset)


%% Pre-allocate working storage
selected = false(1,length(stacknums));
err = NaN(sz(4),1);
erru = err;
% err(stacknums,:) = nanmean(options.err(:,stacknums))';
shift = nan(size(options.shift'));
firststacknums = funch.choose_good_firststacks(options,stacknums,options.peristim);
shift(stacknums,:) = options.shift(:,firststacknums)';
udata = cell(sz(4),1);
setappdata(handles.figMain,'selected',selected);
setappdata(handles.figMain,'err',err)
setappdata(handles.figMain,'erru',erru)
setappdata(handles.figMain,'shift',shift)
setappdata(handles.figMain,'udata',udata)

%% Prepare graphs
set([handles.axesSel handles.axesErr handles.axesShift handles.axesUDiff],'XLim',[-1 sz(end)+1],'TickDir','out')
set([handles.axesErr handles.axesShift handles.axesUDiff],'XTick',[])
set(handles.axesSel,'YTick',[])
ylabel(handles.axesErr,'Error')
ylabel(handles.axesShift,'Shift')
ylabel(handles.axesUDiff,'u-difference')
% Plot bars corresponding to stimulation
stim = header.stim_lookup;
if ~isempty(stim)
  n_stim = max(stim);
else
  n_stim = 1;
end
col = repmat(0.8,n_stim,3);  % use gray for all bars
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

%% Prepare GUI elements
set(handles.editFrame,'String',num2str(round(sz(3)/2)));

% Set up defaults for parameters
maskops = struct('margin_frac',0.1,'margin_pixel',2);
setappdata(handles.figMain,'maskops',maskops);

% Set parameters for warping
stk = single(smm(:,:,:,options.base_stacknum));
pyramid = array_restrict_schedule(size(stk),header);
warpops.gridsize = cat(1,pyramid.sz);
warpops.mspixval = NaN;%mean(stk(:).^2);
warpops.min_pixels = 7;
warpops.min_g_pixels = 2;
warpops.gap_data = 3;
setappdata(handles.figMain,'warpops',warpops);
pcoptions = register_phasecorr_initialize(stk,struct('pyramid',pyramid));
h = fspecial('gaussian',5,1);
pcoptions.smooth = @(im)imfilter(im,h);
setappdata(handles.figMain,'pcoptions',pcoptions);

% Choose default command line output for register_movie_dfof_gui
handles.output = hObject;

% Indicate the "parent" GUI so register_gui_utilities knows what to do
handles.caller = 'register_movie_dfof_gui';

% Update handles structure
guidata(hObject, handles);

% If user specified a save file at launch, load the data
if loadsaved
  rmg_loadsaved(handles,savefile);
end

% Do initial plotting
rmg_plot_stacknums(handles)
rmg_plot_results(handles)

% UIWAIT makes register_movie_dfof_gui wait for user response (see UIRESUME)
% uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = register_movie_dfof_gui_OutputFcn(hObject, eventdata, handles) 
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
  options = getappdata(handles.figMain,'options');
  h = findobj(handles.axesSel,'type','line');
  h = setdiff(h,[handles.dragStart handles.dragStop]);
  if ~isempty(h)
    delete(h)
  end
  hline = line([stacknums; stacknums],repmat([0;1],1,length(stacknums)),'Color','b','Parent',handles.axesSel,'ButtonDownFcn',@rmg_select);
  % Set color on prototypes to magenta
  prototypeIndex = findainb(options.prototype(options.include),stacknums);
  set(hline(prototypeIndex),'Color','m')
  % Set color on base stack to red
  if ~isempty(intersect(options.base_stacknum,stacknums))
    baseIndex = findainb(options.base_stacknum,stacknums);
    set(hline(baseIndex),'Color','m')
  else
    line([0 0]+options.base_stacknum,[0 1],'Color','r','HitTest','off','Parent',handles.axesSel); % not selectable if not in stacknums
  end
  % Set selection state
  selected = getappdata(handles.figMain,'selected');
  offon = {'off';'on'};
  set(hline,{'Selected'},offon(1+selected(:)));
  if (sum(selected) == 1)
    set(handles.btnPlayStack,'Enable','on');
  else
    set(handles.btnPlayStack,'Enable','off');
  end
  btns = [handles.btnPlayMovie handles.btnRegisterWarp];
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
  
function trange = rmg_get_trange(handles)
  % Get start/stop values
  tmp = get(handles.dragStart,'XData');
  trange = tmp(1);
  tmp = get(handles.dragStop,'XData');
  trange(2) = tmp(1);


function rmg_loadsaved(handles,filename)
  s = load(filename);
  setappdata(handles.figMain,'err',s.err);
  setappdata(handles.figMain,'options',s.options);
  setappdata(handles.figMain,'stacknums',s.stacknums);
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
   
  


% --- Executes on button press in btnRegisterRigid.
function btnRegisterRigid_Callback(hObject, eventdata, handles)
% hObject    handle to btnRegisterRigid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Get the selected stack numbers
  funch = getappdata(handles.figMain,'funch');
  funch.register_rigid(handles);
  rmg_plot_results(handles)


% --- Executes on button press in btnToggleSelect.
function btnToggleSelect_Callback(hObject, eventdata, handles)
% hObject    handle to btnToggleSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  selected = getappdata(handles.figMain,'selected');
  stacknums = getappdata(handles.figMain,'stacknums');
  trange = rmg_get_trange(handles);
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
  options = getappdata(handles.figMain,'options');
  udata = getappdata(handles.figMain,'udata');
  stacknums = getappdata(handles.figMain,'stacknums');
  smm = getappdata(handles.figMain,'smm');
  filename = smm.filename;
  header = smm.header;
  pixel_spacing = header.pixel_spacing;
  trange = rmg_get_trange(handles);
  [savefilename,pathname] = uiputfile;
  if ~isequal(savefilename,0)
    save([pathname savefilename],'filename','pixel_spacing','stacknums','err','erru','options','udata','trange');
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

% --- Executes on button press in btnPlayMovie.
function btnPlayMovie_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  funch = getappdata(handles.figMain,'funch');
  funch.play_movie(handles)

  

% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  [filename,pathname] = uigetfile;
  if isequal(filename,0)
    return
  end
  if isequal(filename,0)
    return
  end
  rmg_loadsaved(handles,[pathname filename]);
  selected = getappdata(handles.figMain,'selected');
  selected(:) = false;
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
    case {'backspace','delete'}
      % Note: these next lines assume you're going to modify some aspect of the
      % stack numbers.
      stacknums = getappdata(handles.figMain,'stacknums');
      selected = getappdata(handles.figMain,'selected');
      err = getappdata(handles.figMain,'err');
      shift = getappdata(handles.figMain,'shift');
      selsn = stacknums(selected);
      err(selsn) = nan;
      shift(selsn,:) = nan;
      switch eventdata.Key
        case {'backspace','delete'}
          stacknums(selected) = [];
          selected = false(size(stacknums));
        otherwise
          %display(eventdata.Key)
      end
      setappdata(handles.figMain,'stacknums',stacknums);
      setappdata(handles.figMain,'selected',selected);
      setappdata(handles.figMain,'err',err);
      setappdata(handles.figMain,'shift',shift);
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
  % Eliminate the prototypes
  options = getappdata(handles.figMain,'options');
  stacknums = setdiff(stacknums,options.prototype);
  shift_all = getappdata(handles.figMain,'shift');
  shift = shift_all(stacknums,:);
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
  err = getappdata(handles.figMain,'err');
  
  % If the user is doing multigrid, check to see if we have rmg_params
  % and create it if not.
  % Note that because we use multiple bases, the image data will change,
  % but the mask will not.
  warpops = getappdata(handles.figMain,'warpops');
  if strcmp(warpops.method,'Multigrid least-squares')
    if (warpops.createmask || ~isappdata(handles.figMain,'rmg_params'))
      funch.create_mask(handles)
    end
    rmg_params = getappdata(handles.figMain,'rmg_params');
    rmg_params.gap_data = warpops.gap_data;
    rmg_params.wcycle = warpops.wcycle;
  end
  
  smm = getappdata(handles.figMain,'smm');
  last_base = [];

  pgoptions = struct('max',length(stacknums),'progress',0);
  pgoptions = progress_bar(pgoptions);
  tic
  % Loop over selected stacks
  for i = 1:length(stacknums)
    thisstacknum = stacknums(i);
    % Determine the initial guess
    u = funch.u_initial_guess(udata,thisstacknum,base_stacknum,pcoptions,warpops);
    % Get the stack data
    stk = funch.fetch_dfof(smm,options,thisstacknum);
    % Get the base/prototype data, if necessary
    prototype_stacknum = funch.get_prototype_stacknum(handles,thisstacknum);
    if ~isequal(prototype_stacknum,last_base)
      fixed = funch.fetch_dfof(smm,options,prototype_stacknum);
      if strcmp(warpops.method,'Multigrid least-squares')
        % Replace the fixed image data
        rmg_params = register_multigrid_options_expand(rmg_params,double(fixed));
        if any(warpops.presmooth > 0)
          rmg_params_smooth = register_multigrid_options_compress(rmg_params);
          rmg_params_smooth = register_multigrid_options_expand(rmg_params_smooth,imfilter_gaussian_mex(double(fixed),warpops.presmooth));
        end
      end
    end
    % Compute the relative shift, and the rigid-registration error
    thisshift = shift(i,:)-shift_all(prototype_stacknum,:);
    if isnan(err(thisstacknum))
      stkrr = image_shift(stk,thisshift);
      err(thisstacknum) = nanmean((stkrr(:)-fixed(:)).^2);
    end
    
    % Do the registration
    pcoptions.shift = thisshift;
    switch warpops.method
      case 'Multigrid least-squares'
        rmg_params.shift = thisshift;
        stk = double(stk);
        umg = double(funch.upc2mg(u,pcoptions.sz_spatial));
        if any(warpops.presmooth > 0)
          stksm = imfilter_gaussian_mex(stk,warpops.presmooth);
          [umg,thiserr] = register_multigrid_vcycle(umg,stksm,warpops.lambda,rmg_params_smooth);
          if warpops.refit
            [umg,thiserr] = register_multigrid_vcycle(umg,stk,warpops.lambda,rmg_params);
          end
        else
          [umg,thiserr] = register_multigrid_vcycle(umg,stk,warpops.lambda,rmg_params);
        end
        u = funch.umg2pc(single(umg),pcoptions.sz_spatial);
        % Recalculate the error with only the data portion; this is for
        % consistency with other methods.
        stkw = register_phasecorr_warp(u,single(stk),pcoptions);
        thiserr = nanmean((fixed(:)-stkw(:)).^2);
      case 'Phase correlation'
        stk = single(stk);
        stk(isnan(stk)) = 0;
        pcoptions.lambda = warpops.lambda;
        if isempty(u) || isempty(u{1})
          % Do multiresolution phasecorr
          multires_params = struct('desired_ulevel',warpops.gap_data+1,...
            'display',true);
          [u,~,thiserr] = register_phasecorr_multires(fixed,stk,pcoptions,multires_params);
          thiserr = thiserr(end);
        else
          % Do phasecorr only at highest resolution
          u = register_phasecorr_improve(u,fixed,stk,pcoptions);
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
      setappdata(handles.figMain,'err',err);
      setappdata(handles.figMain,'erru',erru);
      setappdata(handles.figMain,'udata',udata);
      rmg_plot_results(handles)
      pgoptions.progress = find(thisstacknum == stacknums);
      pgoptions = progress_bar(pgoptions);
      tic
    end
  end



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



function editCoRegisterTShift_Callback(hObject, eventdata, handles)
% hObject    handle to editCoRegisterTShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCoRegisterTShift as text
%        str2double(get(hObject,'String')) returns contents of editCoRegisterTShift as a double


% --- Executes during object creation, after setting all properties.
function editCoRegisterTShift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCoRegisterTShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxShowDfof.
function checkboxShowDfof_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowDfof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowDfof


% --- Executes on button press in btnSelect.
function btnSelect_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%   smm = getappdata(handles.figMain,'smm');
  options = getappdata(handles.figMain,'options');
  choices = options.stimuli(options.include);
  [selection,ok] = listdlg('ListString',choices);
  if ok
    selected = getappdata(handles.figMain,'selected');
    selected(:) = false;
    stacknums = getappdata(handles.figMain,'stacknums');
    onset = getappdata(handles.figMain,'onset');
    selstk = cat(1,onset{selection});
    selected(findainb(unique(selstk),stacknums)) = true;
    setappdata(handles.figMain,'selected',selected);
    rmg_plot_stacknums(handles);
  end



function editDfofMax_Callback(hObject, eventdata, handles)
% hObject    handle to editDfofMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDfofMax as text
%        str2double(get(hObject,'String')) returns contents of editDfofMax as a double


% --- Executes during object creation, after setting all properties.
function editDfofMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDfofMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
  