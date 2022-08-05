function varargout = pd_grid_gui(varargin)
% PD_GRID_GUI: tune all the actuators based on imaging data
% Syntax:
%   pd_grid_gui(specfilename (or structure),camera_DM_filename,camera_filename,snipindex,imagingConfiguration,options)
% OR
%   pd_grid_gui(prevfitfilename)
%
% Options:
%   Zindex (default 1:20): if empty, skip Zernike fitting and just do Gaussians
%   layout_fcn (default @mirao52_layout): DM layout function
%   old fits? (act_fits)
% To avoid blocking, no outputs: need "save to file"
  
% Zvalue for each actuator is stored in a cell array; that way it can be
% empty, a vector (meaning we have made a guess from the grid), or a
% matrix (meaning we have fit it to the data).
  
% PD_GRID_GUI M-file for pd_grid_gui.fig
%      PD_GRID_GUI, by itself, creates a new PD_GRID_GUI or raises the existing
%      singleton*.
%
%      H = PD_GRID_GUI returns the handle to a new PD_GRID_GUI or the handle to
%      the existing singleton*.
%
%      PD_GRID_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PD_GRID_GUI.M with the given input arguments.
%
%      PD_GRID_GUI('Property','Value',...) creates a new PD_GRID_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pd_grid_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pd_grid_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pd_grid_gui

% Last Modified by GUIDE v2.5 17-Jul-2009 11:30:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pd_grid_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pd_grid_gui_OutputFcn, ...
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


% --- Executes just before pd_grid_gui is made visible.
function pd_grid_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pd_grid_gui (see VARARGIN)

% Choose default command line output for pd_grid_gui
handles.output = hObject;

if (length(varargin) == 1)
  % We're loading an old saved file
  s = load(varargin{1});
  fn = fieldnames(s);
  for indx = 1:length(fn)
    setappdata(handles.figGrid,fn{indx},s.(fn{indx}));
  end
  load(varargin{1}); % reload to get into workspace
  [ny,nx] = size(act_grid);
else
  % We're starting from scratch, parse the arguments
  imagingConfiguration = varargin{5};
  snipindex = varargin{4};
  if (length(varargin) < 6)
    options = struct;
  else
    options = varargin{6};
  end
  options = default(options,'Zindex',1:20,'layout_fcn',@mirao52_layout);
  Zindex = options.Zindex;

  % load the specfile, containing information about disk-storage format
  data_specs = varargin{1};
  if ~isstruct(data_specs)
    data_specs = load(data_specs);
  end
  
  setappdata(handles.figGrid,'data_specs',data_specs);
  setappdata(handles.figGrid,'snipindex',snipindex);
  setappdata(handles.figGrid,'DMfilename',varargin{2});
  setappdata(handles.figGrid,'basefilename',varargin{3});
  setappdata(handles.figGrid,'Zindex',Zindex);

  % Set up the pupil
  imsz = [length(snipindex{1}),length(snipindex{2})];
  [pupildata.H0,pupildata.rho,pupildata.theta] = pupil_initialize(imagingConfiguration,imsz);
  setappdata(handles.figGrid,'pupildata',pupildata);
  
  % Create the layout of the DM
  [act_grid,act_xy] = options.layout_fcn();
  setappdata(handles.figGrid,'act_grid',act_grid);
  setappdata(handles.figGrid,'act_xy',act_xy);
  [ny,nx] = size(act_grid);

  % Create a blank space for actuator fitting results  
  if isfield(options,'act_fits')
    setappdata(handles.figGrid,'act_fits',options.act_fits);
  else
    setappdata(handles.figGrid,'act_fits',cell(ny,nx));
  end
end

% Some convenience parameters
v = data_specs.vol_tune;
setappdata(handles.figGrid,'v',v);
[tmp,baseIndex] = mindist(0,v);
setappdata(handles.figGrid,'baseIndex',baseIndex);

% Compute Zernike basis
zernv = zernike_values(pupildata.rho,pupildata.theta,Zindex);
setappdata(handles.figGrid,'zernv',zernv);

% Prepare for any Gaussian computations
config = struct('mode','vgaussian singlechannel','rho',pupildata.rho,'theta',pupildata.theta,'v',v);
gaussian_phimodel = phi_parametrizations(config);
setappdata(handles.figGrid,'gaussian_phimodel',gaussian_phimodel);

% Create grid axes
splitx = SplitAxesEvenly(nx,0);
splity = SplitAxesEvenly(ny,0);
hax = SplitGrid(splitx,splity,ones(1,nx),ones(1,ny),handles.axesGridShell);
handles.grid_axes = hax;

% Set the click properties of the grid axes
set(handles.grid_axes,'NextPlot','replacechildren',...
		  'ButtonDownFcn',@grid_click);
delete(hax(isnan(act_grid)));

% Update handles structure
guidata(hObject, handles);

replot_grid(handles)

% UIWAIT makes pd_grid_gui wait for user response (see UIRESUME)
% uiwait(handles.figGrid);


% --- Outputs from this function are returned to the command line.
function varargout = pd_grid_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnAlign.
function btnAlign_Callback(hObject, eventdata, handles)
% hObject    handle to btnAlign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  act_fits = getappdata(handles.figGrid,'act_fits');
  act_grid = getappdata(handles.figGrid,'act_grid');
  act_xy = getappdata(handles.figGrid,'act_xy');
  pupildata = getappdata(handles.figGrid,'pupildata');
  keepFlag = ~isnan(act_grid);
  hvalid = handles.grid_axes(keepFlag);
  fit_valid = act_fits(keepFlag);
  selected = get(hvalid,'Selected');
  isselected = cellfun(@(s) strcmp(s,'on'),selected) & cellfun(@(s) isfield(s,'gaussianP'),fit_valid);
  gaussianPC = cellfun(@(s) s.gaussianP,fit_valid(isselected),'UniformOutput',false);
  gaussianP = cat(1,gaussianPC{:});
  actnum = act_grid(keepFlag);
  actnum = actnum(isselected);  % as an index, this assumes actuators are numbered 1, 2, 3, 4, ...
  A = grid_cp2A(act_xy(actnum,:),gaussianP(:,[2 3]));  % fit the centers to the grid of actuators
  T = maketform('affine',A);
  % Calculate position of all actuators
  act_phase = tformfwd(T,act_xy);
  % Display
  figure
  pupilimg = fftshift(pupildata.H0);
  pupilimg = repmat(pupilimg,[1 1 3]);
  pupilimg(pupilimg == 0) = 0.8;
  xr = [-1 1] * max(pupildata.rho(:));
  image(xr,xr,pupilimg)
  for k = 1:size(act_phase,1)
    text(act_phase(k,1),act_phase(k,2),num2str(k),'HorizontalAlignment','center','VerticalAlignment','middle');
  end
  phaserange = [min(act_phase); max(act_phase)];
  set(gca,'XLim',phaserange(:,1),'YLim',phaserange(:,2),'Visible','off')
  axis equal
  % Create a guess for all actuators that haven't yet been fit
  meangaussian = mean(gaussianP,1);
  for indx = 1:numel(act_fits)
    thisActuator = act_grid(indx);
    if ~isnan(thisActuator) && ~isfield(act_fits{indx},'gaussianP')
      act_fits{indx}.gaussianP = [meangaussian(1) act_phase(thisActuator,:), meangaussian(4)];
    end
  end
  setappdata(handles.figGrid,'act_fits',act_fits);
  setappdata(handles.figGrid,'actxy2phase',T);
  

% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Collect the input information
  s.data_specs = getappdata(handles.figGrid,'data_specs');
  s.DMfilename = getappdata(handles.figGrid,'DMfilename');
  s.basefilename = getappdata(handles.figGrid,'basefilename');
  s.snipindex = getappdata(handles.figGrid,'snipindex');
  s.Zindex = getappdata(handles.figGrid,'Zindex');
  s.pupildata = getappdata(handles.figGrid,'pupildata');
  s.act_grid = getappdata(handles.figGrid,'act_grid');
  s.act_xy = getappdata(handles.figGrid,'act_xy');
  % Collect the results
  s.act_fits = getappdata(handles.figGrid,'act_fits');
  % Get the filename from the user
  [filename,pathname] = uiputfile('*.mat');
  if filename ~= 0
    save([pathname filesep filename],'-struct','s')
  end
  
% --- Executes on button press in btnAutofit.
function btnAutofit_Callback(hObject, eventdata, handles)
% hObject    handle to btnAutofit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  act_fits = getappdata(handles.figGrid,'act_fits');
  act_grid = getappdata(handles.figGrid,'act_grid');
  %gaussian_phimodel = getappdata(handles.figGrid,'gaussian_phimodel');
  pupildata = getappdata(handles.figGrid,'pupildata');
  data_specs = getappdata(handles.figGrid,'data_specs');
  Zindex = getappdata(handles.figGrid,'Zindex');
  v = getappdata(handles.figGrid,'v');
  n_axes = numel(handles.grid_axes);
  for indx = 1:n_axes
    this_hax = handles.grid_axes(indx);
    if ~ishandle(this_hax)
      continue;  % NaN-corresponding grid positions
    end
    this_fit = act_fits{indx};
    this_actuator = act_grid(indx);
    data_present = any(data_specs.act_num == this_actuator);
    if ~data_present
      continue
    end
    % Get data
    actnum = act_grid(indx);
    [imbase,imDM] = load_image_data(actnum,handles);
    this_fit.interactive = false;
    this_fit.Zindex = Zindex;
    this_fit = pd_tune_single_actuator(imbase,imDM,v,pupildata,this_fit);
    act_fits{indx} = this_fit;
    % Update each time, in case this gets interrupted
    setappdata(handles.figGrid,'act_fits',act_fits);
    replot_grid(handles)
  end
    
      

  % Stuff from gaussian_gui
%   gaussianP = getappdata(handles.figPDG,'gaussianP');
%   im = getappdata(handles.figPDG,'im');
%   pupildata = getappdata(handles.figPDG,'pupildata');
%   gaussianP = calcphi2d(im,pupildata.H0,gaussianP,phimodel);
% 

% --- Executes on button press in btnClearGrid.
function btnClearGrid_Callback(hObject, eventdata, handles)
% hObject    handle to btnClearGrid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in btnClearSelected.
function btnClearSelected_Callback(hObject, eventdata, handles)
% hObject    handle to btnClearSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  act_fits = getappdata(handles.figGrid,'act_fits');
  act_grid = getappdata(handles.figGrid,'act_grid');
  keepFlag = ~isnan(act_grid);
  hvalid = handles.grid_axes(keepFlag);
  fit_valid = act_fits(keepFlag);
  selected = get(hvalid,'Selected');
  isselected = cellfun(@(s) strcmp(s,'on'),selected) & cellfun(@(s) isfield(s,'gaussianP'),fit_valid);
  fit_valid{isselected} = struct;
  act_fits(keepFlag) = fit_valid;
  setappdata(handles.figGrid,'act_fits',act_fits);
  replot_grid(handles);

  
%% More user interaction functions
function grid_click(object,src,event)
  % Determine the click type
  handles = guidata(object);
  seltype = get(handles.figGrid,'SelectionType');
  % Get the actuator number
  gridIndex = handle2index(object,handles);
  actnum = index2actnum(gridIndex,handles);
  switch seltype
    case 'normal'
      % Collect any previously-defined solution for the actuator
      act_fits = getappdata(handles.figGrid,'act_fits');
      prev_fit = act_fits{gridIndex};
      if isempty(prev_fit)
        prev_fit = struct;
      end
      prev_fit = default(prev_fit,'Zindex',getappdata(handles.figGrid,'Zindex'));
      % Load the image data
      [imbase,imDM] = load_image_data(actnum,handles);
      % Get other things we'll need
      v = getappdata(handles.figGrid,'v');
      pupildata = getappdata(handles.figGrid,'pupildata');
      % Tune this actuator
      fit = pd_tune_single_actuator(imbase,imDM,v,pupildata,prev_fit);
%       if ~isfield(fit,'Zvalue') || isempty(fit.Zvalue)
%         % User closed figure without optimizing
%         return
%       end
      act_fits{gridIndex} = fit; %% START HERE (change act_fits to struct, elim act_gaussian, flags to say what has happened)
      setappdata(handles.figGrid,'act_fits',act_fits);
      replot_grid(handles)
    case 'alt'
      state = get(handles.grid_axes(gridIndex),'Selected');
      if (strcmp(state,'off'))
        state = 'on';
      else
        state = 'off';
      end
      set(handles.grid_axes(gridIndex),'Selected',state);
  end

%% Plotting functions
function replot_grid(handles)
  N = numel(handles.grid_axes);
  act_fits = getappdata(handles.figGrid,'act_fits');
  act_grid = getappdata(handles.figGrid,'act_grid');
  data_specs = getappdata(handles.figGrid,'data_specs');
  gaussian_phimodel = getappdata(handles.figGrid,'gaussian_phimodel');
  for indx = 1:N
    this_hax = handles.grid_axes(indx);
    if ~ishandle(this_hax)
      continue;  % NaN-corresponding grid positions
    end
    this_fit = act_fits{indx};
    this_actuator = act_grid(indx);
    data_present = any(data_specs.act_num == this_actuator);
    if data_present && isfield(this_fit,'Zvalue') && ~isempty(this_fit.Zvalue)
      % There has been a full Zernike-based fit. Extract the part that
      % depends linearly upon voltage
      Zvalue = this_fit.Zvalue;
      v = getappdata(handles.figGrid,'v');
      Zv = linregress(v,Zvalue);
      % Compute the phase corresponding to the voltage-dependent part
      nZ = length(Zv);
      zernv = getappdata(handles.figGrid,'zernv');
      szz = size(zernv);
      phiV = zeros(szz(1:2));
      for k = 1:nZ
        phiV = phiV + Zv(k) * zernv(:,:,k);
      end
    elseif isfield(this_fit,'gaussianP') && data_present
      % We have a Gaussian estimate of the aberration. Convert to a full
      % phase
      phiV = gaussian_phimodel.phicoef(this_fit.gaussianP);
    else
      % No fit yet, just display the actuator number
      set(this_hax,'XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
      text(0.5,0.5,num2str(this_actuator),'Parent',this_hax,'HorizontalAlignment','center',...
        'VerticalAlignment','middle','HitTest','off');
      if ~data_present
        set(this_hax,'Color',[0.8 0.8 0.8],'ButtonDownFcn','');
      end
      continue; % skip the rest of the plotting stuff
    end
    % Plot it
    axes(this_hax);
    himg = imagesc(fftshift(phiV));
    set(himg,'HitTest','off');  % Don't block clicks from axes
    axis tight
    set(gca,'YDir','reverse','XTick',[],'YTick',[])
    % Test to see if this is just a gaussian-estimate, or a real fit;
    % if it's just a gaussian-estimate, turn the grid on as a visual
    % indicator.
%     if (length(this_fit) == 1)
%       % Gaussian-estimate
%       grid on
%     else
%       % A real fit
%       grid off
%       set(this_hax,'XTick',[],'YTick',[]);
%     end
  end

      
  
%% Utility functions
function [imbase,imDM] = load_image_data(actuator,handles)
  data_specs = getappdata(handles.figGrid,'data_specs');
  imbase = pdocpi_fread(getappdata(handles.figGrid,'basefilename'), ...
    data_specs,actuator);
  imDM = pdocpi_fread(getappdata(handles.figGrid,'DMfilename'), ...
    data_specs,actuator);
  snipindex = getappdata(handles.figGrid,'snipindex');
  imbase = double(imbase(snipindex{:},:));
  imDM = double(imDM(snipindex{:},:));

function [iy,ix] = handle2index(hax,handles)
  if (nargout == 1)
    % Single-index syntax
    iy = find(hax == handles.grid_axes);
  else
    % Double-index syntax
    [iy,ix] = find(hax == handles.grid_axes);
  end
  
function actnum = index2actnum(varargin)
% either index2actnum(index,handles)
%     or index2actnum(iy,ix,handles)
  act_grid = getappdata(varargin{end}.figGrid,'act_grid');
  if (length(varargin) == 2)
    % Single-index syntax
    actnum = act_grid(varargin{1});
  else
    % Double-index syntax
    actnum = act_grid(varargin{1:2});
  end
  
function A = grid_cp2A(input,base)
  T = cp2tform(input,base,'similarity');
  A = T.tdata.T(:,1:2);
  % Now fit to a model in which each coordinate is first allowed to be
  % scaled before applying the rotation
  s = sqrt(sum(A(1:2,:).^2,1));
  theta = atan2(A(2,1),A(1,1));
  p0 = [s theta A(3,:)];
  p = fminunc(@(x) grid_cp_misalignment(input,base,x),p0);
  A = grid_p2A(p);
  
function A = grid_p2A(p)
  s = p(1:2);
  theta = p(3);
  t = p(4:5);
  st = sin(theta);
  ct = cos(theta);
  A = [s(1)*ct -s(2)*st; s(1)*st s(2)*ct; t];
  
function err = grid_cp_misalignment(input,base,p)
  A = grid_p2A(p);
  input_transformed = input * A(1:2,1:2) + repmat(A(3,:),size(input,1),1);
  tmp = input_transformed - base;
  err = sum(tmp(:).^2);
  


