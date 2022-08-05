function varargout = psf_gui(varargin)
% psf_gui: calculate the PSF by interactively clicking on beads
%
% Syntax:
%   [psf,roicenters] = psf_gui(stk)
% where
%   stk is a three-dimensional array of images of diffraction-limited beads
% and
%   psf is a three-dimensional array containing the average bead image
%   roicenters is an n_rois-by-3 matrix containing the coordinates of the
%     bead centers
%
% Click in the vicinity of any beads. The peak of each bead will be entered
% into the ROI list. The summed intensity in each ROI can be used to guess
% whether you have one bead or two (if the illumination is even). You can
% delete ROIs that you don't like.
%
% Keyboard shortcuts (only work when the figure is selected, e.g., click on
% background):
%   pageup/pagedown/uparrow/downarrow: change the frame
%   delete/backspace: delete any selected ROIs
% Selections of ROIs can be controlled using standard ctrl-click,
% shift-click.

% Copyright 2011 by Timothy E. Holy

 
% PSF_GUI M-file for psf_gui.fig
%      PSF_GUI, by itself, creates a new PSF_GUI or raises the existing
%      singleton*.
%
%      H = PSF_GUI returns the handle to a new PSF_GUI or the handle to
%      the existing singleton*.
%
%      PSF_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PSF_GUI.M with the given input arguments.
%
%      PSF_GUI('Property','Value',...) creates a new PSF_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before psf_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to psf_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help psf_gui

% Last Modified by GUIDE v2.5 17-Sep-2011 08:58:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @psf_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @psf_gui_OutputFcn, ...
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


% --- Executes just before psf_gui is made visible.
function psf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to psf_gui (see VARARGIN)

stk = varargin{1};
sz = size(stk);
setappdata(handles.figMain,'stk',stk);
setappdata(handles.figMain,'sz',sz);
roilist = zeros(0,3);
setappdata(handles.figMain,'roilist',roilist)

% Set up the display
handles.image = image(stk(:,:,1),'Parent',handles.axesImage);
set(handles.axesImage,'DataAspectRatio',[1 1 1],'TickDir','out','CLim',[0 max(stk(:))]);
set(handles.image,'CDataMapping','scaled','ButtonDownFcn',@newroi)
set(handles.figMain,'Colormap',gray(256));

% Set up keypress handling
set(handles.figMain,'KeyPressFcn',@parsekey)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes psf_gui wait for user response (see UIRESUME)
uiwait(handles.figMain);


% --- Outputs from this function are returned to the command line.
function varargout = psf_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout = handles.output;
delete(handles.figMain)


function showimage(handles)
  frameIndex = str2double(get(handles.editFrame,'String'));
  stk = getappdata(handles.figMain,'stk');
  set(handles.image,'CData',stk(:,:,frameIndex))
  showrois(handles)
  
function showrois(handles)
  hline = findobj(handles.axesImage,'type','line');
  if ~isempty(hline)
    delete(hline)
  end
  if get(handles.checkboxShowROIs,'Value')
    roilist = getappdata(handles.figMain,'roilist');
    rng = extract_range(handles);
    frameIndex = str2double(get(handles.editFrame,'String'));
    for i = 1:size(roilist,1)
      if frameIndex >= roilist(i,3)+rng(1,3) && frameIndex <= roilist(i,3)+rng(2,3)
        x = roilist(i,2)+rng(:,2);
        y = roilist(i,1)+rng(:,1);
        line(x([1 2 2 1 1]), y([1 1 2 2 1]),'Color','w')
      end
    end
  end

function parsekey(obj,evt)
  handles = guidata(obj);
  switch evt.Key
    case {'pagedown','downarrow'}
      frameIndex = str2double(get(handles.editFrame,'String'));
      frameIndex = max(1,frameIndex-1);
      set(handles.editFrame,'String',num2str(frameIndex));
      showimage(handles)
    case {'pageup','uparrow'}
      frameIndex = str2double(get(handles.editFrame,'String'));
      sz = getappdata(handles.figMain,'sz');
      frameIndex = min(sz(3),frameIndex+1);
      set(handles.editFrame,'String',num2str(frameIndex));
      showimage(handles)
    case {'delete','backspace'}
      indx = get(handles.listboxROI,'Value');
      if ~isempty(indx)
        roilist = getappdata(handles.figMain,'roilist');
        roilist(indx,:) = [];
        setappdata(handles.figMain,'roilist',roilist);
        set(handles.listboxROI,'Value',[]);
        show_roilist(handles)
      end
  end

function rng = extract_range(handles)
  rng = [str2double(get(handles.editSzY,'String')) ...
    str2double(get(handles.editSzX,'String')) ...
    str2double(get(handles.editSzZ,'String'))];
  rng = round([-1; 1]*rng);
  
function psf = calculate_psf(handles)
  fprintf('Registering ROIs...');
  roilist = getappdata(handles.figMain,'roilist');
  stk = getappdata(handles.figMain,'stk');
  if isempty(roilist)
    psf = [];
    return
  end
  % First pass, compute the psf as the average of the ROIs
  rng = extract_range(handles);
  snip = snip_rois(stk,roilist,rng);
  psf = nanmean(snip,4);
  shift = array_findpeak(psf);
  roilist = bsxfun(@plus,roilist,shift);
  % Second and later passes, register each individual bead image to the
  % psf, then re-snip and re-average
  n_registrations = 3;
  n_rois = size(roilist,1);
  for i = 1:n_registrations
    fixed_pad = register_translate_pad(psf);
    fixed_fft = register_translate_nancorr(fixed_pad);
    for j = 1:n_rois
      moving_pad = register_translate_pad(snip(:,:,:,j));
      shift = register_translate_nancorr(fixed_fft,moving_pad,struct('dx_max',rng(2,:)));
      fprintf('ROI %d, shift ',j);
      fprintf('%g ',shift);
      fprintf('\n');
      roilist(j,:) = roilist(j,:) + shift;
    end
    snip = snip_rois(stk,roilist,rng);
    psf = nanmean(snip,4);
    shift = array_findpeak(psf);
    fprintf('Mean shifted to peak, shift ',j);
    fprintf('%g ',shift);
    fprintf('\n');
    roilist = bsxfun(@plus,roilist,shift);
  end
  fprintf('done.\n');
  setappdata(handles.figMain,'roilist',roilist);
  psf = nanmean(snip,4);
  
  
function snip = snip_rois(stk,roilist,rng)
  n_rois = size(roilist,1);
  sz = size(stk);
  snip = nan([diff(rng,1,1)+1 n_rois],'single');
  for i = 1:n_rois
    thisroi = roilist(i,:);
    thisrng = bsxfun(@plus,rng,thisroi);
    x = cell(1,3);
    keepFlag = cell(1,3);
    if all(thisroi == round(thisroi))
      % All ranges are integers, so can just cut it out
      for j = 1:3
        x{j} = thisrng(1,j):thisrng(2,j);
        keepFlag{j} = x{j} > 0 & x{j} <= sz(j);
        x{j} = x{j}(keepFlag{j});
      end
      snip(keepFlag{:},i) = stk(x{:});
    else
      % Use subpixel interpolation
      % Calling qinterp_grid_inverse on the whole stack could take a long
      % time! So just snip out a slightly larger-than-necessary region and
      % invert the quadratic interpolation over this region.
      szout = diff(rng,1,1)+1;
      xfirst = zeros(1,3);
      for j = 1:3
        x{j} = round(thisrng(1,j)-3:thisrng(2,j)+3);
        keepFlag{j} = x{j} > 0 & x{j} <= sz(j);
        x{j} = x{j}(keepFlag{j});
        xfirst(j) = x{j}(1);
      end
      tmpsnip = stk(x{:});
      % Snip from the square root of the image, so that we interpolate in a way
      % that respects shot noise
      [keepFlag2,thissnip] = image_snip_qinterp(szout,thisroi-(szout-1)/2-xfirst,sqrt(tmpsnip));
      snip(keepFlag2{:},i) = thissnip.^2;
    end
  end
  
function newroi(obj,~)
  hax = get(obj,'Parent');
  pt = get(hax,'CurrentPoint');
  pt = round(pt(1,[2 1]));  % use coordinate order, not x then y
  hfig = get(hax,'Parent');
  handles = guidata(hfig);
  frameIndex = str2double(get(handles.editFrame,'String'));
  rng = extract_range(handles);
  % Snip out around the click point
  stk = getappdata(handles.figMain,'stk');
  roitmp = [pt frameIndex];
  snip = snip_rois(stk,roitmp,rng);
  % Find the peak of the ROI (no need for subpixel accuracy)
  dx = array_findpeak(snip,struct('subpixel_mode',''));
  roitmp = roitmp+dx;
  % Store the result
  roilist = getappdata(handles.figMain,'roilist');
  roilist(end+1,:) = roitmp;
%   if ~isempty(roilist)
%     % Register the new ROI to the first one
%     roitmp = [roilist(1,:); pt frameIndex];
%     stk = getappdata(handles.figMain,'stk');
%     rng = extract_range(handles);
%     snip = snip_rois(stk,roitmp,rng);
%     fixedfft = register_translate_nancorr(snip(:,:,:,1));
%     shift = register_translate_nancorr(fixedfft,snip(:,:,:,2),struct('dx_max',rng(2,:)));
%     pt = pt + round(shift(1:2));
%     frameIndex = frameIndex + round(shift(3));
%   end
%   roilist(end+1,:) = [pt frameIndex];
  setappdata(handles.figMain,'roilist',roilist);
  show_roilist(handles)
  
function show_roilist(handles)
  roilist = getappdata(handles.figMain,'roilist');
  stk = getappdata(handles.figMain,'stk');
  rng = extract_range(handles);
  snip = snip_rois(stk,roilist,rng);
  I = nansum(nansum(nansum(snip,1),2),3);
  roitmp = [roilist I(:)/mean(I(:))];
  set(handles.listboxROI,'String',num2str(roitmp));
  
  
% --- Executes on selection change in listboxROI.
function listboxROI_Callback(hObject, eventdata, handles)
% hObject    handle to listboxROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxROI contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxROI


% --- Executes during object creation, after setting all properties.
function listboxROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFrame_Callback(hObject, eventdata, handles)
% hObject    handle to editFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrame as text
%        str2double(get(hObject,'String')) returns contents of editFrame as a double
sz = getappdata(handles.figMain,'sz');
frameIndex = str2double(get(handles.editFrame,'String'));
frameIndex = max(1,min(frameIndex,sz(3)));
set(handles.editFrame,'String',num2str(frameIndex));
showimage(handles)


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



function editSzY_Callback(hObject, eventdata, handles)
% hObject    handle to editSzY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSzY as text
%        str2double(get(hObject,'String')) returns contents of editSzY as a double
showrois(handles)
show_roilist(handles)


% --- Executes during object creation, after setting all properties.
function editSzY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSzY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSzX_Callback(hObject, eventdata, handles)
% hObject    handle to editSzX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSzX as text
%        str2double(get(hObject,'String')) returns contents of editSzX as a double
showrois(handles)
show_roilist(handles)


% --- Executes during object creation, after setting all properties.
function editSzX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSzX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSzZ_Callback(hObject, eventdata, handles)
% hObject    handle to editSzZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSzZ as text
%        str2double(get(hObject,'String')) returns contents of editSzZ as a double
showrois(handles)
show_roilist(handles)

% --- Executes during object creation, after setting all properties.
function editSzZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSzZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles.output = {calculate_psf(handles)};
  roilist = getappdata(handles.figMain,'roilist');
  handles.output{2} = roilist;
  guidata(handles.figMain,handles)
  uiresume(handles.figMain)


% --- Executes on button press in pushbuttonContrast.
function pushbuttonContrast_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  cl = get(handles.axesImage,'CLim');
  frameIndex = str2double(get(handles.editFrame,'String'));
  stk = getappdata(handles.figMain,'stk');
  cl = imrangegui(stk(:,:,frameIndex),cl);
  if ~isempty(cl)
    set(handles.axesImage,'CLim',cl)
  end
  showimage(handles)


% --- Executes on button press in checkboxShowROIs.
function checkboxShowROIs_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxShowROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxShowROIs
