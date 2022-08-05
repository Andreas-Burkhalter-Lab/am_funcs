function varargout = register_dfof_validate_stacks(varargin)
% register_dfof_validate_stacks: a GUI for checking for bad frames in base stacks
%
% This GUI determines (via find_bad_frames) the bad frames in a movie. It
% then gives the user the opportunity to select a "prototypical"
% response for each stimulus.
%
% Syntax:
%   s = register_dfof_validate_stacks
%   s = register_dfof_validate_stacks(smm)
%   s = register_dfof_validate_stacks(filename)
% where the argument, if supplied, determines the data file that will
% be analyzed. 
%
% The first thing you will be prompted to do is examine the base stack, to
% make sure it has no irregularities.  By default, the base stack
% will be the stack preceding a stimulus near the middle of the experiment.
% This stack will be shown in the red channel; the two stacks immediately
% before this stack will be shown in green and blue.  A movie window will
% pop up, and you should play through the movie to make sure you don't see
% big splotches of red or cyan.  If you think it's fine, answer the prompt
% with "y" (yes). Otherwise, you will be prompted to try another base
% stack.
%
% The output "s" is a structure that lists the stimuli, the base
% stack # (indicated as the valve onset), the choices of temporal
% shifts for "pre-stimulus" stacks ("background") and "peri-stimulus"
% stacks ("foreground"), and the complete record of information about bad
% frames.
%
%  s = register_dfof_validate_stacks(...,options)
% This version allows you to supply the output of find_bad_frames as
% fields of the structure "options": base_stacknum, err, shift, and isbad.
% This can save you from having to call find_bad_frames a second time
% (which is the slow step in this GUI).
%
% See also: find_bad_frames, register_movie_dfof_gui.

% Copyright 2011 by Timothy E. Holy

% REGISTER_DFOF_VALIDATE_STACKS M-file for register_dfof_validate_stacks.fig
%      REGISTER_DFOF_VALIDATE_STACKS, by itself, creates a new REGISTER_DFOF_VALIDATE_STACKS or raises the existing
%      singleton*.
%
%      H = REGISTER_DFOF_VALIDATE_STACKS returns the handle to a new REGISTER_DFOF_VALIDATE_STACKS or the handle to
%      the existing singleton*.
%
%      REGISTER_DFOF_VALIDATE_STACKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in REGISTER_DFOF_VALIDATE_STACKS.M with the given input arguments.
%
%      REGISTER_DFOF_VALIDATE_STACKS('Property','Value',...) creates a new REGISTER_DFOF_VALIDATE_STACKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before register_dfof_validate_stacks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to register_dfof_validate_stacks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help register_dfof_validate_stacks

% Last Modified by GUIDE v2.5 25-Jan-2011 08:53:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @register_dfof_validate_stacks_OpeningFcn, ...
                   'gui_OutputFcn',  @register_dfof_validate_stacks_OutputFcn, ...
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


% --- Executes just before register_dfof_validate_stacks is made visible.
function register_dfof_validate_stacks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to register_dfof_validate_stacks (see VARARGIN)

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

%% Validate base stack and do bad-frame analysis
if ~isfield(options,'err') || ~isfield(options,'isbad') || ~isfield(options,'base_stacknum')
  [isbad,base_stacknum,err,shift] = find_bad_frames(smm,options);
  if isempty(isbad)
    delete(handles.figMain);
    return;
  end
  options.base_stacknum = base_stacknum;
  options.err = err;
  options.shift = shift;
  options.isbad = isbad;
end

header = smm.header;
funch = register_gui_utilities;
% prestim = find(header.stim_lookup(1:end-1) == 0 & header.stim_lookup(2:end) > 0);
% prestim = prestim(:);

%% Find the default "stimulus-specific prototype stacks"
% These should have no bad frames, and be as close to the base stack as
% possible
% deltat = [rdvs_get_deltat(handles.editPre) rdvs_get_deltat(handles.editPeri)];
[prototype,onset,ustim] = funch.choose_good_onsets(smm,options,rdvs_get_deltat(handles.editPre),rdvs_get_deltat(handles.editPeri));
stimuli = header.stim_labels(ustim);
n_stimuli = length(onset);
include = true(1,n_stimuli);

setappdata(handles.figMain,'smm',smm);
setappdata(handles.figMain,'options',options);
setappdata(handles.figMain,'onset',onset);
setappdata(handles.figMain,'prototype',prototype);
setappdata(handles.figMain,'stimuli',stimuli);
setappdata(handles.figMain,'include',include);
setappdata(handles.figMain,'funch',funch);

%% Set up graphical elements
set(handles.popupmenuStimulus,'String',stimuli);
if (n_stimuli > 0)
  yl = [0.5 n_stimuli+0.5];
else
  yl = [0.5 1.5];
end
set(handles.axesErr,'YTick',1:n_stimuli,'YTickLabel',stimuli);
set(handles.axesShift,'YTick',[]);
set([handles.axesErr handles.axesShift],'YLim',yl,'YDir','reverse');
xlabel(handles.axesErr,'Mismatch')
xlabel(handles.axesShift,'Shift')
rdvs_set_popup_stacknums(handles,1)

rdvs_plot_results(handles)

% UIWAIT makes register_dfof_validate_stacks wait for user response (see UIRESUME)
uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = register_dfof_validate_stacks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if nargin > 2 && ~isempty(handles) && isfield(handles,'output')
  varargout = {handles.output};
  if ishandle(handles.figMain)
    delete(handles.figMain);
  end
else
  varargout = {[]};
end


function sn = rdvs_get_deltat(h)
  sn = str2num(get(h,'String')); %#ok<ST2NM>
  if isempty(sn)
    errordlg('The entry cannot be interpreted as a list of integers.')
  end

function rdvs_set_popup_stacknums(handles,stimIndex)
  prototype = getappdata(handles.figMain,'prototype');
  onset = getappdata(handles.figMain,'onset');
  trial = find(prototype(stimIndex)==onset{stimIndex});
  set(handles.popupmenuPrototype,'String',num2str(onset{stimIndex}(:)),...
    'Value',trial);

function stk = rdvs_get_current_stacks(handles,index)
  smm = getappdata(handles.figMain,'smm');
  prototype = getappdata(handles.figMain,'prototype');
  stimIndex = get(handles.popupmenuStimulus,'Value');
  curstack = prototype(stimIndex);
  stk = smm(:,:,:,curstack+index);
 
function err = rdvs_calculate_err(err_tot,prototype,deltat)
  n_stimuli = length(prototype);
  n_frames = size(err_tot,1);
  err = zeros(n_frames,n_stimuli);
  for stimIndex = 1:n_stimuli
    err(:,stimIndex) = mean(err_tot(:,prototype(stimIndex)+deltat),2);
  end
  
function rdvs_plot_results(handles)
  cla(handles.axesErr)
  cla(handles.axesShift)
  options = getappdata(handles.figMain,'options');
  prototype = getappdata(handles.figMain,'prototype');
  deltat = rdvs_get_deltat(handles.editPre);
  errPre = rdvs_calculate_err(options.err,prototype,deltat);
  deltat = rdvs_get_deltat(handles.editPeri);
  errPeri = rdvs_calculate_err(options.err,prototype,deltat);
  % Plot the error
  yval = repmat((1:length(prototype))',1,size(errPre,1))';
  line(errPre,yval,'Parent',handles.axesErr,'LineStyle','none','Marker','x','Color','g');
  line(errPeri,yval,'Parent',handles.axesErr,'LineStyle','none','Marker','x','Color','r');
  % Plot the shift, converting to physical units
  shift = options.shift(:,prototype);
  smm = getappdata(handles.figMain,'smm');
  header = smm.header;
  pixel_spacing = header.pixel_spacing;
  line(shift(1,:)*pixel_spacing(1),1:length(prototype),'Parent',handles.axesShift,'LineStyle','none','Marker','x','Color','b');
  line(shift(2,:)*pixel_spacing(2),1:length(prototype),'Parent',handles.axesShift,'LineStyle','none','Marker','x','Color','g');
  line(shift(3,:)*pixel_spacing(3),1:length(prototype),'Parent',handles.axesShift,'LineStyle','none','Marker','x','Color','r');

function rdvs_reset_stacknums_prototypes(handles)
  prototype = getappdata(handles.figMain,'prototype');
  smm = getappdata(handles.figMain,'smm');
  options = getappdata(handles.figMain,'options');
  funch = getappdata(handles.figMain,'funch');
  deltat = [rdvs_get_deltat(handles.editPre) rdvs_get_deltat(handles.editPeri)];
  [~,onset] = funch.choose_good_onsets(smm,options,deltat);
  setappdata(handles.figMain,'onset',onset);
  % Choose the closest onset to existing prototype
  for stimIndex = 1:length(onset)
    [~,minIndex] = min(abs(onset{stimIndex}-prototype(stimIndex)));
    prototype(stimIndex) = onset{stimIndex}(minIndex);
  end
  setappdata(handles.figMain,'prototype',prototype);
  stimIndex = get(handles.popupmenuStimulus,'Value');
  rdvs_set_popup_stacknums(handles,stimIndex)
  
% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles.output = [];
  guidata(handles.figMain,handles);
  uiresume(handles.figMain);
  
% --- Executes on button press in btnDone.
function btnDone_Callback(hObject, eventdata, handles)
% hObject    handle to btnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  options = getappdata(handles.figMain,'options');
  prototype = getappdata(handles.figMain,'prototype');
  stimuli = getappdata(handles.figMain,'stimuli');
  include = getappdata(handles.figMain,'include');
  pre = rdvs_get_deltat(handles.editPre);
  peri = rdvs_get_deltat(handles.editPeri);
  dfof_clip = str2num(get(handles.editDfofClip,'String'));
  if isempty(pre) || isempty(peri)
    errordlg('Cannot finish until the pre- and peri-stimulus stacks are properly defined');
    return
  end
  handles.output = struct(...
    'stimuli',{stimuli},...
    'include',include,...
    'prototype',prototype,...
    'prestim',pre,...
    'peristim',peri,...
    'dfof_clip',dfof_clip,...
    'base_stacknum',options.base_stacknum,...
    'isbad',options.isbad,...
    'err',options.err);
  if isfield(options,'shift')
    handles.output.shift = options.shift;
  end
  guidata(handles.figMain,handles)
  uiresume(handles.figMain);

  
% --- Executes on selection change in popupmenuStimulus.
function popupmenuStimulus_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuStimulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  stimIndex = get(hObject,'Value');
  rdvs_set_popup_stacknums(handles,stimIndex);
  include = getappdata(handles.figMain,'include');
  set(handles.checkboxInclude,'Value',include(stimIndex));

% --- Executes during object creation, after setting all properties.
function popupmenuStimulus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuStimulus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editPre_Callback(hObject, eventdata, handles)
% hObject    handle to editPre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Ensure that the result is interpretable
  rdvs_get_deltat(handles.editPre);
  rdvs_reset_stacknums_prototypes(handles);
  rdvs_plot_results(handles);

% --- Executes during object creation, after setting all properties.
function editPre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPeri_Callback(hObject, eventdata, handles)
% hObject    handle to editPeri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Ensure the result is interpretable
  rdvs_get_deltat(handles.editPeri);
  rdvs_reset_stacknums_prototypes(handles);
  rdvs_plot_results(handles)

% --- Executes during object creation, after setting all properties.
function editPeri_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPeri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenuPrototype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuPrototype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  prototype = getappdata(handles.figMain,'prototype');
  stimIndex = get(handles.popupmenuStimulus,'Value');
  prototypeIndex = get(hObject,'Value');
  onsets = get(hObject,'String');
  prototype(stimIndex) = str2double(onsets(prototypeIndex,:));
  setappdata(handles.figMain,'prototype',prototype)
  rdvs_plot_results(handles)

% --- Executes during object creation, after setting all properties.
function popupmenuPrototype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuPrototype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnPlayDfof.
function btnPlayDfof_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayDfof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  bgIndex = rdvs_get_deltat(handles.editPre);
  bg = mean(single(rdvs_get_current_stacks(handles,bgIndex)),4);
  fgIndex = rdvs_get_deltat(handles.editPeri);
  fg = mean(single(rdvs_get_current_stacks(handles,fgIndex)),4);
  dfof = fg./bg-1;
  dfofmax = str2double(get(handles.editDfofMax,'String'));
  dfof(dfof > dfofmax) = dfofmax;
  dfof(dfof < -dfofmax) = -dfofmax;
  % Normalize each
  dfof = dfof/dfofmax;
  bg = bg/max(bg(:));
  % Create the color stack
  rgb = colorize_pixels(bg,dfof);
  mplay(permute(rgb,[1 2 4 3]))


% --- Executes on button press in btnPlayCompare.
function btnPlayCompare_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlayCompare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  stk = rdvs_get_current_stacks(handles,-1);
  smm = getappdata(handles.figMain,'smm');
  options = getappdata(handles.figMain,'options');
  fixed = smm(:,:,:,options.base_stacknum);
  rgb = cat(4,fixed,stk);
  rgb = single(rgb)/single(max(rgb(:)));
  rgb(:,:,:,3) = 0;
  mplay(permute(rgb,[1 2 4 3]))
  
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

% --- Executes on button press in btnCalculate.
function btnCalculate_Callback(hObject, eventdata, handles)
% hObject    handle to btnCalculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  % Collect current settings
  smm = getappdata(handles.figMain,'smm');
  base_stacknum = getappdata(handles.figMain,'base_stacknum');
  prototype = getappdata(handles.figMain,'prototype');
  prototype_old = getappdata(handles.figMain,'prototype_old');
  onset = getappdata(handles.figMain,'onset');
  err = getappdata(handles.figMain,'err');
  shift = getappdata(handles.figMain,'shift');
  bgIndex = rdvs_get_deltat(handles.editPre);
  fgIndex = rdvs_get_deltat(handles.editPeri);
  % Get the base stack
  fixed = single(smm(:,:,:,base_stacknum));
  % For each prototype # that has changed, recalculate the shift and the
  % mismatch
  pgoptions = struct('max',length(prototype));
  tic
  for stimIndex = 1:length(prototype)
    if (toc > 3 || stimIndex == length(prototype))
      pgoptions.progress = stimIndex;
      pgoptions = progress_bar(pgoptions);
      tic
    end
    if (prototype(stimIndex) == prototype_old(stimIndex))
      continue
    end
    thisstacknum = onset{stimIndex}(trial(stimIndex));
    bg = mean(single(smm(:,:,:,thisstacknum+bgIndex)),4);
    [imr,shift(stimIndex,:)] = register_rigid(fixed,bg,struct('subpixel',false));
    thiserr = nanmean(nanmean((imr - fixed).^2,1),2);
    fg = mean(single(smm(:,:,:,thisstacknum+fgIndex)),4);
    imr = image_shift(fg,shift(stimIndex,:));
    thiserr = thiserr + nanmean(nanmean((imr - fixed).^2,1),2);
    err(stimIndex,:) = thiserr(:)';
  end
  prototype_old = prototype;
  setappdata(handles.figMain,'shift',shift);
  setappdata(handles.figMain,'err',err);
  setappdata(handles.figMain,'prototype_old',prototype_old);
  % Clear the old plots
  cla(handles.axesErr)
  cla(handles.axesShift)
  % Plot the error
  for i = 1:size(err,2)
    line(err(:,i),1:length(prototype),'Parent',handles.axesErr,'LineStyle','none','Marker','x','Color','b');
  end
  % Plot the shift, converting to physical units
  header = smm.header;
  pixel_spacing = header.pixel_spacing;
  line(shift(:,1)*pixel_spacing(1),1:length(prototype),'Parent',handles.axesShift,'LineStyle','none','Marker','x','Color','b');
  line(shift(:,2)*pixel_spacing(2),1:length(prototype),'Parent',handles.axesShift,'LineStyle','none','Marker','x','Color','g');
  line(shift(:,3)*pixel_spacing(3),1:length(prototype),'Parent',handles.axesShift,'LineStyle','none','Marker','x','Color','r');
  
    


% --- Executes on button press in checkboxInclude.
function checkboxInclude_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxInclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  include = getappdata(handles.figMain,'include');
  stimIndex = get(handles.popupmenuStimulus,'Value');
  include(stimIndex) = get(handles.checkboxInclude,'Value');
  setappdata(handles.figMain,'include',include);



function editDfofClip_Callback(hObject, eventdata, handles)
% hObject    handle to editDfofClip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDfofClip as text
%        str2double(get(hObject,'String')) returns contents of editDfofClip as a double
  n = str2num(get(hObject,'String'));
  if (length(n) ~= 2)
    errordlg('Must insert a 2-vector, [min max]');
  end

% --- Executes during object creation, after setting all properties.
function editDfofClip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDfofClip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
