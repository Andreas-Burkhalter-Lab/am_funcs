function varargout = threshold_stack_gui(varargin)
% THRESHOLD_STACK_GUI M-file for threshold_stack_gui.fig
%      THRESHOLD_STACK_GUI, by itself, creates a new THRESHOLD_STACK_GUI or raises the existing
%      singleton*.
%
%      H = THRESHOLD_STACK_GUI returns the handle to a new THRESHOLD_STACK_GUI or the handle to
%      the existing singleton*.
%
%      THRESHOLD_STACK_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHOLD_STACK_GUI.M with the given input arguments.
%
%      THRESHOLD_STACK_GUI('Property','Value',...) creates a new THRESHOLD_STACK_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before threshold_stack_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to threshold_stack_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help threshold_stack_gui

% Last Modified by GUIDE v2.5 27-Mar-2012 12:25:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @threshold_stack_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @threshold_stack_gui_OutputFcn, ...
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


% --- Executes just before threshold_stack_gui is made visible.
function threshold_stack_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to threshold_stack_gui (see VARARGIN)

% initialize appdata
fig = handles.figure1;
img = varargin{1};
ops = varargin{2};
s_lab = varargin{3};

% set up "max project" view
ops.stim_to_use(end+1)= length(s_lab)+1; % so as (hopefully) not to muddle with another stim :-)
s_lab{end+1} = 'max_project';
set(handles.stimlabel_popup,'string',s_lab(ops.stim_to_use));
img.dfof(:,:,:,end+1) = max(img.dfof,[],4);
img.zscore(:,:,:,end+1) = max(img.zscore,[],4);

this_stim = ops.stim_to_use(end);
framenum = ceil(size(img.dfof,3)/2);
this_label = s_lab(this_stim);
this_frame = img.(ops.metric)(:,:,framenum,find(ops.stim_to_use==this_stim));
dfof_thresh = repmat(ops.dfof_thresh,[size(img.dfof,3),length(s_lab)]);
zscore_thresh = repmat(ops.zscore_thresh,[size(img.zscore,3),length(s_lab)]);
ops = default(ops,'dfof_thresh',dfof_thresh);
ops = default(ops,'zscore_thresh',zscore_thresh);
ops = default(ops,'edge_mask',[0 0 0 0]);

setappdata(fig,'img',img);
setappdata(fig,'ops',ops);
setappdata(fig,'s_lab',s_lab);
setappdata(fig,'this_stim',this_stim);
setappdata(fig,'framenum',framenum);
setappdata(fig,'this_frame',this_frame);
setappdata(fig,'handles',handles);
setappdata(fig,'dfof_thresh',ops.dfof_thresh);
setappdata(fig,'zscore_thresh',ops.zscore_thresh);
setappdata(fig,'displayed',ops.metric);
setappdata(fig,'normalized',false);
setappdata(fig,'edge_mask',ops.edge_mask);


% set up the image in stack_axis
axes(handles.stack_axis);
handles.image_axis = imagesc(img.(ops.metric)(:,:,framenum,find(ops.stim_to_use ==this_stim))); set(gca,'visible','off');
colormap(gray);
guidata(hObject,handles);
dclims = [-0.03 0.03];
setappdata(fig,'dclims',dclims);
zclims = [-0.1 1];
setappdata(fig,'zclims',zclims);
if ~isempty(strmatch(ops.metric,'dfof'))
	set(handles.dfof_displayed_radio,'value',1);
	UpdateClims(hObject,eventdata);
elseif ~isempty(strmatch(ops.metric,'zscore'))
	set(handles.zscore_displayed_radio,'value',1);
	UpdateClims(hObject,eventdata);
end

% initializing static text boxes
imsz = ops.max_size(1:2);
set(handles.imsz_text,'string',['actual: ' num2str(imsz) ' pix']);
set(handles.stack_axis,'units','pixels');
axpos = get(handles.stack_axis,'position');
set(handles.stack_axis,'units','normalized');
set(handles.edit6, 'string', [num2str(ops.edge_mask(1))]);
set(handles.edit7, 'string', [num2str(ops.edge_mask(2))]);
set(handles.edit8, 'string', [num2str(ops.edge_mask(3))]);
set(handles.edit9, 'string', [num2str(ops.edge_mask(4))]);
if any(axpos(3:4)<imsz)
	set(handles.imdisp_text,'string',['displayed: ' num2str(axpos(3:4)) ' pix']);
else
	set(handles.imdisp_text,'string',['displayed: ' num2str(imsz) ' pix']);  
end
set(handles.filt_text,'string',['filter: [' num2str(ops.sigma) ']']);
UpdateLabel(hObject,eventdata);

% set radio button to default;

set(fig,'KeyReleaseFcn',@ChangeFrame,'BusyAction','cancel')

% Choose default command line output for threshold_stack_gui
handles.output = {[]};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes threshold_stack_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = threshold_stack_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
warning('off','all');
fig = hObject;
uiwait(fig);
if ~ishandle(hObject)
	varargout = {[]};
else
	ops = getappdata(fig,'ops');
    maskEdges = getappdata(fig, 'maskEdges');  %Deleter after chedking for match
	dfof_thresh = getappdata(fig,'dfof_thresh');
	zscore_thresh = getappdata(fig,'zscore_thresh');
    edge_mask = getappdata(fig,'edge_mask');
    ops.dfof_thresh = dfof_thresh;
	ops.zscore_thresh = zscore_thresh;
    ops.edge_mask = edge_mask;
    ops.maskEdges = maskEdges; %Delete after checking for match
	varargout = {ops};
	close(fig);
end
warning('on');


% --- Executes on button press in done_button.
function done_button_Callback(hObject, eventdata, handles)
% hObject    handle to done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get(hObject,'parent');
uiresume(fig);


% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(get(hObject,'parent'),'dfof_thresh',[]);
setappdata(get(hObject,'parent'),'zscore_thresh',[]);
close(get(hObject,'parent'));

% 
function CloseFig

function dfof_thresh_box_Callback(hObject, eventdata, handles)
% hObject    handle to dfof_thresh_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = get_parent_fig(hObject);
ops = getappdata(fig,'ops');
this_stim = getappdata(fig,'this_stim');
framenum = getappdata(fig,'framenum');
si = find(ops.stim_to_use == this_stim);
dfof_thresh = getappdata(fig,'dfof_thresh');
dfof_thresh(framenum,si) = str2num(get(hObject,'string'));
setappdata(fig,'dfof_thresh',dfof_thresh);
UpdateThresh(fig,eventdata);

% --- Executes during object creation, after setting all properties.
function dfof_thresh_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfof_thresh_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function zscore_thresh_box_Callback(hObject, eventdata, handles)
% hObject    handle to zscore_thresh_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = get_parent_fig(hObject);
ops = getappdata(fig,'ops');
this_stim = getappdata(fig,'this_stim');
framenum = getappdata(fig,'framenum');
si = find(ops.stim_to_use == this_stim);
zscore_thresh = getappdata(fig,'zscore_thresh');
zscore_thresh(framenum,si) = str2num(get(hObject,'string'));
setappdata(fig,'zscore_thresh',zscore_thresh);
UpdateThresh(fig,eventdata);


% --- Executes during object creation, after setting all properties.
function zscore_thresh_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zscore_thresh_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function clim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to clim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = get(hObject,'parent');
displayed = getappdata(fig,'displayed');
clim = str2num(get(hObject,'string'));
if ~isempty(strmatch(displayed,'dfof','exact'))
	if length(clim)==2 && clim(1) < clim(2)
		setappdata(fig,'dclims',clim);
	else
		set(hObject,'string',num2str(getappdata(fig,'dclims')));
	end
elseif ~isempty(strmatch(displayed,'zscore','exact'))
		if length(clim)==2 && clim(1) < clim(2)
		setappdata(fig,'zclims',clim);
	else
		set(hObject,'string',num2str(getappdata(fig,'zclims')));
	end
end
UpdateClims(fig,eventdata);

% --- Executes during object creation, after setting all properties.
function clim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stimlabel_popup.
function stimlabel_popup_Callback(hObject, eventdata, handles)
% hObject    handle to stimlabel_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = get(hObject,'parent');
ops = getappdata(fig,'ops');
newval = get(hObject,'value'); 
setappdata(fig,'this_stim',ops.stim_to_use(newval));
UpdateLabel(fig,eventdata);


% --- Executes during object creation, after setting all properties.
function stimlabel_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimlabel_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function frame_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimlabel_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function frame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = get(hObject,'parent');
fh = guihandles(fig);
framenum = getappdata(fig,'framenum');
newframe = str2num(get(hObject,'string'));
if newframe ~= framenum
	img = getappdata(fig,'img');
	if newframe < size(img.dfof,3) && newframe > 0
		setappdata(fig,'framenum',newframe);
		UpdateFrame(fig,eventdata);
	end
end

% --- Executes during object creation, after setting all properties.
function frame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% BEGIN NON-UITOOL FUNCTIONS
% --- Executes when user presses a key
function ChangeFrame(hObject,eventdata)
fig = hObject;
ops = getappdata(fig,'ops');
img = getappdata(fig,'img');
sz = size(img.dfof);
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
switch eventdata.Key
	case 'uparrow'
		cur = find(ops.stim_to_use==this_stim);
		if length(ops.stim_to_use) > cur
			setappdata(fig,'this_stim',ops.stim_to_use(cur+1));
			UpdateLabel(hObject,eventdata)
		end
	case 'downarrow'
		cur = find(ops.stim_to_use==this_stim);
		if cur > 1
			setappdata(fig,'this_stim',ops.stim_to_use(cur-1));
			UpdateLabel(hObject,eventdata);
		end
	case 'leftarrow'
		if framenum > 1
			setappdata(fig,'framenum',framenum-1);
			UpdateThresh(fig,eventdata);
		end
	case 'rightarrow'
		if framenum < sz(3)
			setappdata(fig,'framenum',framenum+1);
			UpdateThresh(fig,eventdata);
		end
end

% --- Updates the displayed stimulus label
function UpdateLabel(hObject,eventdata)
fig = hObject;
fh = guihandles(fig);
ops = getappdata(fig,'ops');
s_lab = getappdata(fig,'s_lab');
this_stim  = getappdata(fig,'this_stim');
set(fh.stimlabel_popup,'value',find(ops.stim_to_use==this_stim));
UpdateThresh(fig,eventdata);

% --- Updates the image clims
function UpdateClims(hObject,eventdata)
fig = hObject;
displayed = getappdata(fig,'displayed');
fh = guihandles(fig);
ax = imgca(fig);
cledit = fh.clim_edit;
dclims = getappdata(fig,'dclims');
zclims = getappdata(fig,'zclims');
fh.stack_axis = ax;
guidata(fig,fh);
if ~isempty(strmatch(displayed,'dfof','exact'))
	set(ax,'clim', dclims);
	set(cledit,'string',num2str(dclims));
elseif ~isempty(strmatch(displayed,'zscore','exact'))
	set(ax,'clim', zclims);
	set(cledit,'string',num2str(zclims));
end


% --- Updates the image threshold
function UpdateThresh(hObject,eventdata)
fig = hObject;
fh = guihandles(fig);
ops = getappdata(fig,'ops');
img = getappdata(fig,'img');
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
idx = find(ops.stim_to_use == this_stim);
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
metric = ops.metric;
old_displayed = getappdata(fig,'displayed');
setappdata(fig,'init',true);
if isempty(dfof_thresh)
	input = questdlg('No thresholds were found, do you want to begin by calculating per-frame, per-stim auto-thresholds?', 'Auto-threshold?','Yes','No','Yes');
	if ~isempty(strmatch(input,'Yes'))		
		setappdata(fig,'displayed','dfof');
		autoAllBtn_Callback(hObject,eventdata);
	else
		dfof_thresh = zeros([size(img.dfof,3) length(ops.stim_to_use)])+0.02;
		zscore_thresh = zeros([size(img.dfof,3) length(ops.stim_to_use)])+1;
		setappdata(fig,'dfof_thresh',dfof_thresh);
		setappdata(fig,'zscore_thresh',zscore_thresh);
	end
end
if isempty(zscore_thresh)
		setappdata(fig,'displayed','zscore');
		autoAllBtn_Callback(hObject,eventdata);
% 
% 		zscore_thresh = repmat(0.5,[size(img.(ops.metric),3),length(ops.stim_to_use)]);
% 		setappdata(fig,'zscore_thresh',zscore_thresh);
end
setappdata(fig,'init',false);
setappdata(fig,'displayed',old_displayed);
dth = fh.dfof_thresh_box;
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
if idx == length(ops.stim_to_use)
	set(dth,'string',['max - ' num2str(max(dfof_thresh(framenum,1:end-1)))]);
else
	set(dth,'string',num2str(dfof_thresh(framenum,idx)));
end
zth = fh.zscore_thresh_box;
if idx == length(ops.stim_to_use)
	set(zth,'string',['max- ' num2str(max(zscore_thresh(framenum,1:end-1)))]);
else
	set(zth,'string',num2str(zscore_thresh(framenum,idx)));
end
UpdateFrame(hObject,eventdata);

% --- Updates the displayed frame
function UpdateFrame(hObject,eventdata)
fig = hObject;
fh = guihandles(fig);
ops = getappdata(fig,'ops');
edge_mask = getappdata(fig,'edge_mask');
im_sz = ops.max_size;
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
idx = find(ops.stim_to_use == this_stim);
img = getappdata(fig,'img');
displayed = getappdata(fig,'displayed');
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');

% plot img & update stats
h = imgca(fig);
thisimg = img.(displayed)(:,:,framenum,idx);
totalpix = sum(~isnan(thisimg(:)));
maskEdges = false(im_sz(1),im_sz(2));
maskEdges((1+edge_mask(2)):(im_sz(1)-edge_mask(1)),...
    (1+edge_mask(3)):(im_sz(2)-edge_mask(4)))=true; % note plot order, y
setappdata(fig,'maskEdges',maskEdges); %Temp line to check mask in roi_by_imflow
if idx == length(ops.stim_to_use)
	sz = length(ops.stim_to_use);
	imdims = size(img.dfof);
	for i = 1:sz-1	
		dmulti = img.dfof(:,:,framenum,i) > dfof_thresh(framenum,i);
		zmulti = img.zscore(:,:,framenum,i) > zscore_thresh(framenum,i);
		if ~exist('multi','var')
			multi = dmulti.*zmulti;
		else
			multi = multi | dmulti.*zmulti;
        end        
    end

% 		tmp = img.dfof(:,:,framenum,i); dlin(i,:) = tmp(:);
% 		tmp = img.zscore(:,:,framenum,i); zlin(i,:) = tmp(:);
% 	end
% 	dmulti = max(bsxfun(@gt,dlin,dfof_thresh(framenum,1:sz-1)'),[],1);
% 	dmulti = reshape(dmulti,imdims(1:2));
% 	zmulti = max(bsxfun(@gt,zlin,zscore_thresh(framenum,1:sz-1)'),[],1);
% 	zmulti = reshape(zmulti,imdims(1:2));
else
	dmulti = (img.dfof(:,:,framenum,idx)>dfof_thresh(framenum,idx));
	zmulti = (img.zscore(:,:,framenum,idx)>zscore_thresh(framenum,idx));
	multi = dmulti.*zmulti;
end
multi = multi.*maskEdges;
supra = sum(multi(:));
pct = supra/totalpix*100;
set(fh.numpix_text,'string',{'# pixels:', num2str(supra)});
set(fh.pctpix_text,'string',{'% pixels:', num2str(pct,4)});
multi = multi+0.4; multi(multi>1) = 1;
set(findobj(h,'type','image'),'cdata',thisimg.*multi);
set(fh.frame_edit,'string',[num2str(framenum)]);
fh.stack_axis = h;
colormap(gray);



guidata(fig,fh);


% --- Executes on button press in dfof_displayed_radio.
function dfof_displayed_radio_Callback(hObject, eventdata, handles)
% hObject    handle to dfof_displayed_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
fh = guihandles(fig);
displayed = getappdata(fig,'displayed');
if isempty(strmatch(displayed,'dfof','exact'))
	set(hObject,'value',1);
	set(fh.zscore_displayed_radio,'value',0);
	setappdata(fig,'displayed','dfof');
	UpdateClims(fig,eventdata);
	UpdateFrame(fig,eventdata);
else
	set(hObject,'value',1);
	set(fh.zscore_displayed_radio,'value',0);
end


% --- Executes on button press in zscore_displayed_radio.
function zscore_displayed_radio_Callback(hObject, eventdata, handles)
% hObject    handle to zscore_displayed_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
fh = guihandles(fig);
displayed = getappdata(fig,'displayed');
if isempty(strmatch(displayed,'zscore','exact'))
	set(hObject,'value',1);
	set(fh.dfof_displayed_radio,'value',0);
	setappdata(fig,'displayed','zscore');
	UpdateClims(fig,eventdata);
	UpdateFrame(fig,eventdata);
else
	set(hObject,'value',1);
	set(fh.dfof_displayed_radio,'value',0);
end


% --- Executes on button press in frameApplyBtn.
function frameApplyBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
ops = getappdata(fig,'ops');
idx = find(ops.stim_to_use == this_stim);
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
input = questdlg('Really apply this to all frames?','Warning: overwrite?','Yes','No','No');
if ~isempty(strmatch(input,'Yes'))
	dfof_thresh(:,idx) = dfof_thresh(framenum,idx);
	zscore_thresh(:,idx) = zscore_thresh(framenum,idx);
	setappdata(fig,'dfof_thresh',dfof_thresh);
	setappdata(fig,'zscore_thresh',zscore_thresh);
	UpdateFrame(fig,eventdata);
else
	return;
end


% --- Executes on button press in stimApplyBtn.
function stimApplyBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
ops = getappdata(fig,'ops');
idx = find(ops.stim_to_use == this_stim);
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
input = questdlg('Really apply this to all stims?','Warning: overwrite?','Yes','No','No');
if ~isempty(strmatch(input,'Yes'))
	dfof_thresh(framenum,:) = dfof_thresh(framenum,idx);
	zscore_thresh(framenum,:) = zscore_thresh(framenum,idx);setappdata(fig,'dfof_thresh',dfof_thresh);
	setappdata(fig,'zscore_thresh',zscore_thresh);
	UpdateFrame(fig,eventdata);
else
	return;
end


% --- Executes on button press in rangeApplyBtn.
function rangeApplyBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
ops = getappdata(fig,'ops');
idx = find(ops.stim_to_use == this_stim);
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
this_d_thresh = dfof_thresh(framenum,idx);
this_z_thresh = zscore_thresh(framenum,idx);
input = inputdlg({'start:','end'},'Specify Range',[1; 1],{num2str(framenum), num2str(framenum)});
range(1)=str2num(input{1});range(2)=str2num(input{2});
if all(range>0) && all(range<size(dfof_thresh,1)) && range(1) <= range(2)
	dfof_thresh(range(1):range(2),idx) = dfof_thresh(framenum,idx);
	zscore_thresh(range(1):range(2),idx) = zscore_thresh(framenum,idx);
	setappdata(fig,'dfof_thresh',dfof_thresh);
	setappdata(fig,'zscore_thresh',zscore_thresh);
	UpdateFrame(fig,eventdata);
end

% --- Executes on button press in allApplyBtn.
function allApplyBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
framenum = getappdata(fig,'framenum');
this_stim = getappdata(fig,'this_stim');
ops = getappdata(fig,'ops');
idx = find(ops.stim_to_use == this_stim);
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
this_d_thresh = dfof_thresh(framenum,idx);
this_z_thresh = zscore_thresh(framenum,idx);
input = questdlg('Really apply this to all frames, all stim?','Warning: overwrite?','Yes','No','No');
if ~isempty(strmatch(input,'Yes'))
	dfof_thresh(:,:) = dfof_thresh(framenum,idx);
	zscore_thresh(:,:) = zscore_thresh(framenum,idx);
	setappdata(fig,'dfof_thresh',dfof_thresh);
	setappdata(fig,'zscore_thresh',zscore_thresh);
	UpdateFrame(fig,eventdata);
else
	return;
end


% --- Executes on button press in saveThreshBtn.
function saveThreshBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
filename = uiputfile('*.roithresh');
save(filename,'dfof_thresh','zscore_thresh','-mat');
 
% --- Executes on button press in loadThreshBtn.
function loadThreshBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
filename = uigetfile('*.roithresh');
load(filename,'-mat');
if ~exist('dfof_thresh','var') || ~exist('zscore_thresh','var')
	warndlg('The suggested .roithresh file was either the wrong type of file or is corrupt');
end
setappdata(fig,'dfof_thresh',dfof_thresh);
setappdata(fig,'zscore_thresh',zscore_thresh);
UpdateThresh(fig,eventdata);
UpdateFrame(fig,eventdata);


% --- Executes on button press in rawValRadio.
function rawValRadio_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
handles = guihandles(fig);
set(handles.rawValRadio,'value',1);
set(handles.normValRadio,'value',0);
normalized = getappdata(fig,'normalized');
if normalized
	dfof_thresh = getappdata(fig,'dfof_thresh');
	newdfof = getappdata(fig,'dfof_thresh_bkup');
	setappdata(fig,'dfof_thresh',newdfof);
	setappdata(fig,'dfof_thresh_bkup',dfof_thresh);
	zscore_thresh = getappdata(fig,'zscore_thresh');
	newzscore = getappdata(fig,'zscore_thresh_bkup');
	setappdata(fig,'zscore_thresh',newzscore);
	setappdata(fig,'zscore_thresh_bkup',zscore_thresh);
	img = getappdata(fig,'img');
	newim = getappdata(fig,'imbkup');
	setappdata(fig,'imbkup',img);
	setappdata(fig,'img',newim);
	setappdata(fig,'normalized',false);
	UpdateThresh(fig,eventdata);
else
	return;
end
	


% --- Executes on button press in normValRadio.
function normValRadio_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
handles = guihandles(fig);
set(handles.rawValRadio,'value',0);
set(handles.normValRadio,'value',1);
img = getappdata(fig,'img');
dfof_thresh = getappdata(fig,'dfof_thresh');
zscore_thresh = getappdata(fig,'zscore_thresh');
if isempty(getappdata(fig,'imbkup')) && ~getappdata(fig,'normalized');
	setappdata(fig,'imbkup',img);
	setappdata(fig,'normalized',true);
	setappdata(fig,'zscore_thresh_bkup',zscore_thresh)
	setappdata(fig,'dfof_thresh_bkup',dfof_thresh)
else
	dfof_thresh = getappdata(fig,'dfof_thresh');
	newdfof = getappdata(fig,'dfof_thresh_bkup');
	setappdata(fig,'dfof_thresh',newdfof);
	setappdata(fig,'dfof_thresh_bkup',dfof_thresh);
	zscore_thresh = getappdata(fig,'zscore_thresh');
	newzscore = getappdata(fig,'zscore_thresh_bkup');
	setappdata(fig,'zscore_thresh',newzscore);
	setappdata(fig,'zscore_thresh_bkup',zscore_thresh);
	img = getappdata(fig,'img');
	newim = getappdata(fig,'imbkup');
	setappdata(fig,'imbkup',img);
	setappdata(fig,'img',newim);
	setappdata(fig,'normalized',true);
	UpdateThresh(fig,eventdata);
	return;
end

imsz = size(img.dfof); imsz(end)=imsz(end)-1;
tmpdf = single(zeros(imsz));
tmpdz = single(zeros(imsz));
for i = 1:size(img.dfof,4)-1
	tmpd = img.dfof(:,:,:,i);
	tmpz = img.zscore(:,:,:,i);
	dfrange{i} = eval([num2str(min(tmpd(:))-0.001) ':0.001:' num2str(max(tmpd(:))+0.001)]);
	minn.dfof(i) = dfrange{i}(end)+0.001;
	maxn.dfof(i) = dfrange{i}(end)-0.001;
	zrange{i} = eval([num2str(min(tmpz(:))-0.001) ':0.001:' num2str(max(tmpz(:))+0.001)]);
	minn.zscore(i) = zrange{i}(end)+0.001;
	maxn.zscore(i) = zrange{i}(end)-0.001;
	hcdf{i} = histc(tmpd(:),dfrange{i});
	hcdz{i} = histc(tmpz(:),zrange{i});
 	td = fit_gaussian(dfrange{i},hcdf{i});
	tz = fit_gaussian(zrange{i},hcdz{i});
% 	mu.dfof(i) = tmp.center;
	hm.dfof(i) = dfrange{i}(hcdf{i}==max(hcdf{i}));
	hm.zscore(i) = zrange{i}(hcdz{i}==max(hcdz{i}));
	sigma.dfof(i) = td.sigma;
	sigma.zscore(i) = tz.sigma;
	tmpdf(:,:,:,i) = tmpd;%/sigma.dfof(i);
	tmpdz(:,:,:,i) = tmpz;
end
for i = 1:size(img.dfof,4)-1
	tmpdf(:,:,:,i) = tmpdf(:,:,:,i)-hm.dfof(i);
	tmpdz(:,:,:,i) = tmpdz(:,:,:,i)-hm.zscore(i);
	dfdivfac(i) = max([maxn.dfof(i) median(maxn.dfof)])/min([sigma.dfof(i) median(sigma.dfof)]);
	zdivfac(i) = max([maxn.zscore(i) median(maxn.zscore)])/min([sigma.zscore(i) median(sigma.zscore)]);
end
dfdivfac = dfdivfac/median(dfdivfac);
zdivfac = zdivfac/median(zdivfac);
for i = 1:size(img.dfof,4)-1
	tmpdf(:,:,:,i) = tmpdf(:,:,:,i)/dfdivfac(i);
	tmpdz(:,:,:,i) = tmpdz(:,:,:,i)/zdivfac(i);
end

tmpdf(:,:,:,end+1) = max(tmpdf,[],4);
tmpdz(:,:,:,end+1) = max(tmpdz,[],4);

img.dfof = tmpdf;
img.zscore = tmpdz;
setappdata(fig,'img',img);
setappdata(fig,'normalized',true);
UpdateThresh(fig,eventdata);

% --- Executes on button press in autoOneBtn.
function autoOneBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
handles = guihandles(fig);
img = getappdata(fig,'img');
framenum = getappdata(fig,'framenum');
ops = getappdata(fig,'ops');
this_stim = getappdata(fig,'this_stim');
displayed = getappdata(fig,'displayed');
idx = find(ops.stim_to_use == this_stim);
if idx == length(ops.stim_to_use) 
	% try to pick for all stim on this frame
	i = 1:length(ops.stim_to_use)-1;
else
	i = idx;
end
for i = i
	thisim = img.(displayed)(:,:,framenum,i);
	imrange = eval([num2str(min(thisim(:))-0.001) ':0.001:' num2str(max(thisim(:))+0.001)]);
	hi = histc(thisim(:),imrange);
	g = fit_gaussian(imrange,hi);
	mu = g.center;
	sigma = g.sigma;
	tmp = mu+6*sigma;
	if ~isempty(strmatch(displayed,'dfof','exact'))
		dfof_thresh = getappdata(fig,'dfof_thresh');
		if tmp < 0.01
			tmp = 0.01;
		end
		dfof_thresh(framenum,i)=tmp;
		setappdata(fig,'dfof_thresh',dfof_thresh);
	elseif ~isempty(strmatch(displayed,'zscore','exact'))
		zscore_thresh = getappdata(fig,'zscore_thresh');
		if tmp < 0.1
			tmp = 0.1;
		end
		zscore_thresh(framenum,i)=tmp;
		setappdata(fig,'zscore_thresh',zscore_thresh);
	end
end
UpdateThresh(fig,eventdata);



% --- Executes on button press in autoAllBtn.
function autoAllBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
handles = guihandles(fig);
img = getappdata(fig,'img');
framenum = getappdata(fig,'framenum');
ops = getappdata(fig,'ops');
this_stim = getappdata(fig,'this_stim');
displayed = getappdata(fig,'displayed');
idx = find(ops.stim_to_use == this_stim);
if idx == length(ops.stim_to_use) 
	% try to pick for all stim on this frame
	i = 1:length(ops.stim_to_use)-1;
else
	i = idx;
end
N = (length(ops.stim_to_use)-1)*size(img.dfof,3);
prog  = progress(struct('progress',0,'max',N));
count = 0;
if getappdata(fig,'init');
	if isempty(getappdata(fig,'dfof_thresh'))
		setappdata(fig,'dfof_thresh',zeros([size(img.dfof,3) length(ops.stim_to_use)]));
	end
  if isempty(getappdata(fig,'zscore_thresh'))
		setappdata(fig,'zscore_thresh',zeros([size(img.dfof,3) length(ops.stim_to_use)]));
	end
end
for i = i
	for j = 1:size(img.dfof,3)
		thisim = img.(displayed)(:,:,j,i);
		imrange = eval([num2str(min(thisim(:))-0.001) ':0.001:' num2str(max(thisim(:))+0.001)]);
		hi = histc(thisim(:),imrange);
		g = fit_gaussian(imrange,hi);
		mu = g.center;
		sigma = g.sigma;
		tmp = mu+6*sigma;
		if ~isempty(strmatch(displayed,'dfof','exact'))
			dfof_thresh = getappdata(fig,'dfof_thresh');
			if tmp < 0.01
				tmp = 0.01;
			end
			dfof_thresh(j,i)=tmp;
			setappdata(fig,'dfof_thresh',dfof_thresh);
		elseif ~isempty(strmatch(displayed,'zscore','exact'))
			zscore_thresh = getappdata(fig,'zscore_thresh');
			if tmp < 0.1
				tmp = 0.1;
			end
			zscore_thresh(j,i)=tmp;
			setappdata(fig,'zscore_thresh',zscore_thresh);
		end
		count = count+1;
		prog.progress = count;
		prog = progress(prog);
	end
end

prog.progress = -1;
prog = progress(prog); if ishandle(prog.handle); delete(prog.handle); end; clear prog;

if ~getappdata(fig,'init')
	UpdateThresh(fig,eventdata);
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
edge_mask = getappdata(fig,'edge_mask');
edge_mask(1) = str2num(get(hObject,'String'));
setappdata(fig,'edge_mask',edge_mask);
UpdateFrame(fig,eventdata)
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
edge_mask = getappdata(fig,'edge_mask');
edge_mask(2) = str2num(get(hObject,'String'));
setappdata(fig,'edge_mask',edge_mask);
UpdateFrame(fig,eventdata)
% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
edge_mask = getappdata(fig,'edge_mask');
edge_mask(3) = str2num(get(hObject,'String'));
setappdata(fig,'edge_mask',edge_mask);
UpdateFrame(fig,eventdata)
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
edge_mask = getappdata(fig,'edge_mask');
edge_mask(4) = str2num(get(hObject,'String'));
setappdata(fig,'edge_mask',edge_mask);
UpdateFrame(fig,eventdata)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
     