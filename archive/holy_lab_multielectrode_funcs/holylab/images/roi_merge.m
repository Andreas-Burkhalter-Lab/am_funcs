function varargout = roi_merge(varargin)
% ROI_MERGE allows interactive merging of ROIs selected by roi_by_imflow
%  Syntax: roi_merge() - user is asked to select files using uigetfile and default options chosen
%          roi_merge(options) - user can supply an options structure with predefined fields
%
%  Inputs: 
%          options: a structure with the following fields
%                   .imfile: a string with the name of an existing .imagine file
%                   .smm: a smm-object (alternative to .imfile for multi-file experiments)
%                   .roifile: a string with the name of an existing .roi file
%                   .intenfile: a string with the name of an existing .intensity file
%                   .displayfile: a string with the name of the .mat file with image info (i.e. temp file from roi_by_imflow)
%
%  Outputs: not applicable/useful for this (GUIDE-based ROI default)
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2011 Julian P Meeks (Timothy Holy Laboratory)
% 
% Revision History
% 2011_08_29: wrote it (JPM)
%  

% Last Modified by GUIDE v2.5 09-Nov-2011 17:58:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @roi_merge_OpeningFcn, ...
                   'gui_OutputFcn',  @roi_merge_OutputFcn, ...
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


% --- Executes just before roi_merge is made visible.
function roi_merge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% parse the inputs
if isobject(varargin{1})
	options.imfile = [];
	smm = varargin{1};
elseif length(varargin)< 1 || ~isstruct(varargin{1})
	options.imfile = uigetfile('*.imagine','Select a .imagine file');
	smm = stackmm(options.imfile);
else
	options = varargin{1};
	if ~isfield(options,'imfile') && ~isfield(options,'smm')
		options.imfile = uigetfile('*.imagine','Select a .imagine file');
		smm = stackmm(options.imfile);
	end
	if ~isfield(options,'smm') && isfield(options,'imfile')
		smm = stackmm(options.imfile);
	end
end
	
if ~isfield(options,'roifile')
	options.roifile = uigetfile('*.roidef', 'Select a .roidef file');
end
if ~isfield(options,'intenfile')
		options.intenfile = uigetfile('*.intensity', 'Select a .intensity file');
end

% begin setting up application data
fig = hObject;
setappdata(fig,'smm',smm);
tmp = load(options.roifile,'-mat');
% IMPOSE A RULE INTERNALLY: if "tags" was passed as a numeric index, convert it to a string index
% RATIONALE: if, at some point, we want to use something more descriptive for a glomerulus id (i.e. Tom, Tina, etc.)
%            it will be great for this to be implemented without a re-write of the software.  
if all(isnumeric([tmp.roi_defs.label]))
	for i = 1:length(tmp.roi_defs)
		tmp.roi_defs(i).label = num2str(tmp.roi_defs(i).label);
	end
end

setappdata(fig,'roi',tmp);
tmp = load(options.intenfile,'-mat');
setappdata(fig,'inten',tmp);
setappdata(fig,'init_done',false);
setappdata(fig,'tempchar',[]);
setappdata(fig,'needsSave',false);
setappdata(fig,'calcNeighbors',true);

if ~isfield(options,'displayfile')
	ans = questdlg('Do you want to select a display file?','Use pre-calculated display file?','Yes','No','Yes');
	if ~isempty(strmatch(ans,'Yes','exact'))
		options.displayfile = uigetfile('*.mat');
		tmp = load(options.displayfile,'-mat');
		if ~isfield(tmp,'im2flow')
			warning('supplied displayfile was not in useful format.  ignored');
			sz = smm.size;
			img = single(smm(:,:,:,floor(sz(4)/2)));
			img(img==0)=NaN;
			setappdata(fig,'img',img);
			setappdata(fig,'displayfile',[]);
		else
			flow_metric = tmp.options.flow_metric;
			i2f = tmp.im2flow.(flow_metric);
			img = prepDisplayImg(fig,i2f);
			img(img==0)=NaN;
			setappdata(fig,'im2flow',tmp.im2flow);
			setappdata(fig,'img',img);
			setappdata(fig,'displayfile',options.displayfile);
		end
	else
		sz = smm.size;
		img = single(smm(:,:,:,floor(sz(4)/2)));
		img(img==0)=NaN;
		setappdata(fig,'img',img);
		setappdata(fig,'displayfile',[]);
	end
else
	tmp = load(options.displayfile,'-mat');
	if ~isfield(tmp,'im2flow')
		warning('supplied displayfile was not in useful format.  ignored');
		sz = smm.size;
		img = single(smm(:,:,:,floor(sz(4)/2)));
		img(img==0) = NaN;
		setappdata(fig,'img',img);
		setappdata(fig,'displayfile',[]);
	else
		flow_metric = tmp.options.flow_metric;
		i2f = tmp.im2flow.(flow_metric);
		img = prepDisplayImg(fig,i2f);
		img(img==0) = NaN;
		setappdata(fig,'im2flow',tmp.im2flow);
		setappdata(fig,'img',img);
		setappdata(fig,'displayfile',options.displayfile);
	end
end
clear tmp;

% initialize browsing variables
tmp = getappdata(fig,'roi');
setappdata(fig,'curRoi',tmp.roi_defs(1).label);
setappdata(fig,'curFrame',tmp.roi_defs(1).centerInPixels(3));
setappdata(fig,'mergeOrders',struct);
setappdata(fig,'selectedRois',{''});
setappdata(fig,'imgHandle',handles.imgAxis);
setappdata(fig,'stimData',struct);
setappdata(fig,'closest',[]);

% call initialization functions
calcDistances(fig);
calcIntervals(fig);
updateCurRoiList(fig);
updateFrame(fig);

setappdata(fig,'init_done',true);

% Choose default command line output for roi_merge
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = roi_merge_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

% --- Executes on button press in LoadBtn.
function LoadBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
handles = guihandles(fig);
selection = uigetfile('*.roidef','Load ROIdef file');
if ~isempty(selection)
	tmp = load(selection,'-mat');
end
setappdata(fig,'needsSave',false);
setappdata(fig,'mergeOrders',struct);
if all(isnumeric([tmp.roi_defs.label]))
	for i = 1:length(tmp.roi_defs)
		tmp.roi_defs(i).label = num2str(tmp.roi_defs(i).label);
	end
end
setappdata(fig,'roi',tmp);
[tags] = {tmp.roi_defs.label};
setappdata(fig,'tags',tags);
setappdata(fig,'curRoi',tmp.roi_defs(1).label);
setappdata(fig,'selectedRois',{''});
setappdata(fig,'calcNeighbors',true);
set(handles.selRoiList,'string',{''});
set(handles.curRoiList,'string',tags);
set([handles.curRoiList handles.selRoiList],'value',1);
calcIntervals(fig);
updateFrame(fig);

% --- Executes during object creation, after setting all properties.
function curRoiList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in selRoiList.
function selRoiList_Callback(hObject, eventdata, handles)
% nothing here (non-interactive)

% --- Executes during object creation, after setting all properties.
function selRoiList_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in exitBtn.
function exitBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
if fig == 0
	fig = hObject;
end
if getappdata(fig,'needsSave')
	response = questdlg('You have not saved changes, what do you want to do?', 'Warning: unsaved data', 'Return to GUI','Close without saving','Return to GUI');
	if isempty(strmatch(response,'Return to GUI','exact'))
		delete(fig);
	end
else
	delete(fig);
end

% --- Executes on mouse press over axes background.
function ldAx1_ButtonDownFcn(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
roi = getappdata(fig,'roi');
mergeOrders = getappdata(fig,'mergeOrders');
handles = guihandles(fig);
this = getappdata(fig,'ldAx1Idx');
if isempty(this)
	return;
end
thisstr = num2str(roi.roi_defs(this).label);
found = strmatch(thisstr,get(handles.curRoiList,'string'),'exact');
if ~isempty(found)
	set(handles.curRoiList,'value',found);
	curRoiList_Callback(handles.curRoiList,[],[]);
else
	for i = 1:length(mergeOrders)
		found = ~isempty(strmatch(thisstr,mergeOrders(i).from,'exact'));
		if found
			found = mergeOrders(i).to;
			found = strmatch(found,get(handles.curRoiList,'string'),'exact');
			set(handles.curRoiList,'value',found);
			curRoiList_Callback(handles.curRoiList,[],[]);
		end
	end
end

% --- Executes on mouse press over axes background.
function ldAx2_ButtonDownFcn(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
roi = getappdata(fig,'roi');
mergeOrders = getappdata(fig,'mergeOrders');
handles = guihandles(fig);
this = getappdata(fig,'ldAx2Idx');
if isempty(this)
	return;
end
thisstr = num2str(roi.roi_defs(this).label);
found = strmatch(thisstr,get(handles.curRoiList,'string'),'exact');
if ~isempty(found)
	set(handles.curRoiList,'value',found);
	curRoiList_Callback(handles.curRoiList,[],[]);
else
	for i = 1:length(mergeOrders)
		found = ~isempty(strmatch(thisstr,mergeOrders(i).from,'exact'));
		if found
			found = mergeOrders(i).to;
			found = strmatch(found,get(handles.curRoiList,'string'),'exact');
			set(handles.curRoiList,'value',found);
			curRoiList_Callback(handles.curRoiList,[],[]);
		end
	end
end

% --- Executes on mouse press over axes background.
function ldAx3_ButtonDownFcn(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
roi = getappdata(fig,'roi');
mergeOrders = getappdata(fig,'mergeOrders');
handles = guihandles(fig);
this = getappdata(fig,'ldAx3Idx');
if isempty(this)
	return;
end
thisstr = num2str(roi.roi_defs(this).label);
found = strmatch(thisstr,get(handles.curRoiList,'string'),'exact');
if ~isempty(found)
	set(handles.curRoiList,'value',found);
	curRoiList_Callback(handles.curRoiList,[],[]);
else
	for i = 1:length(mergeOrders)
		found = ~isempty(strmatch(thisstr,mergeOrders(i).from,'exact'));
		if found
			found = mergeOrders(i).to;
			found = strmatch(found,get(handles.curRoiList,'string'),'exact');
			set(handles.curRoiList,'value',found);
			curRoiList_Callback(handles.curRoiList,[],[]);
		end
	end
end

% --- Executes on button press in AddBtn.
function AddBtn_Callback(hObject, eventdata, handles)
h = handles.curRoiList;
ed.Character = '+';
curRoiList_KeyPressFcn(h,ed,handles);

% --- Executes on button press in RemBtn.
function RemBtn_Callback(hObject, eventdata, handles)
h = handles.selRoiList;
ed.Key = 'delete';
selRoiList_KeyPressFcn(h,ed,handles);

% --- Executes on selection change in curRoiList.
function curRoiList_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
taglocal = get(hObject,'string');
cur_selected = get(hObject,'value');
match = taglocal{cur_selected};
setappdata(fig,'curRoi',match);
selectedRois = getappdata(fig, 'selectedRois');

updateFrame(fig);
if isempty(selectedRois{1})
	 setappdata(fig,'calcNeighbors',true);
   findCtr(fig);
end

% --- Executes on button press in findCenterBtn.
function findCenterBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
findCtr(fig);

% --- Executes during object creation, after setting all properties.
function curFrame_edit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function curFrame_edit_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
oldFrame = getappdata(fig,'curFrame');
smm = getappdata(fig,'smm');
sz = smm.size;
cur_string = get(hObject,'string');
cur_num = str2num(cur_string);
if cur_num <= sz(3) && cur_num >= 1
	setappdata(fig,'curFrame',cur_num)
	updateFrame(fig);
else
	set(hObject,'string',num2str(oldFrame));
end

% --- Executes when user attempts to  figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
exitBtn_Callback(hObject,eventdata,handles);

% --- Executes on key release with focus on figure1 and none of its controls.
function figure1_KeyReleaseFcn(hObject, eventdata, handles)

fig = hObject;
smm = getappdata(fig,'smm');
sz = smm.size;
switch eventdata.Key
	case 'rightarrow'
		curFrame = getappdata(fig,'curFrame');
		if curFrame+1 <= sz(3)
			setappdata(fig,'curFrame',curFrame+1);
			setappdata(fig,'tempchar',[]);
			updateFrame(fig);
		end
	case 'leftarrow'
		curFrame = getappdata(fig,'curFrame');
		if curFrame-1 >= 1
			setappdata(fig,'curFrame',curFrame-1)
			setappdata(fig,'tempchar',[]);
			updateFrame(fig);
		end
	case 'return'
		curFrame = str2num(getappdata(fig,'tempchar'));
		if curFrame <= sz(3) && curFrame >= 1
			setappdata(fig,'curFrame',curFrame);
			updateFrame(fig);
		end
		setappdata(fig,'tempchar',[]);
	case 'c'
		findCtr(fig);
	otherwise
		tmp = cell(1);
		for i = 1:10
			tmp{i} = num2str(i-1);
		end
	  switch eventdata.Character
			case tmp
				tempchar = getappdata(fig,'tempchar');
				tempchar = [tempchar eventdata.Character];
				setappdata(fig,'tempchar',tempchar);
			otherwise
				setappdata(fig,'tempchar',[]);
		end	
end

% --- Executes on button press in selRoiClearBtn.
function selRoiClearBtn_Callback(hObject, eventdata, handles)
% hObject    handle to selRoiClearBtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig = get_parent_fig(hObject);
setappdata(fig,'selectedRois',{''});
if isempty(getappdata(fig,'curRoi'))
	curList = getappdata(fig,'curList');
	setappdata(fig,'curRoi',curList{1});
end
setappdata(fig,'calcNeighbors',true);
updateSelRoiList(fig);
updateFrame(fig);


% --- Executes on button press in MergeBtn.
function MergeBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
tags = getappdata(fig,'tags');
roi = getappdata(fig,'roi');
selectedRois = getappdata(fig,'selectedRois');
if isempty(selectedRois{1}) 
	return;
else
	for i = 1:length(selectedRois)
		[ismerged(i),prevSel{i},mergeIndex(i)] = isMerged(fig,selectedRois{i});
		if ~isnan(mergeIndex(i))
			mergeOrders = getappdata(fig,'mergeOrders');
			mergeOrders(mergeIndex(i)) = [];
			setappdata(fig,'mergeOrders',mergeOrders);
		end
	end
end

mergeOrders = getappdata(fig,'mergeOrders');
[~,tmpi]=sort(str2double(selectedRois));
if length(fieldnames(mergeOrders)) < 1 && length(selectedRois) > 1
	mergeOrders.to = selectedRois{tmpi(1)};
	mergeOrders.from = unique(selectedRois(tmpi(2:end)));
	if ~isempty(strmatch(mergeOrders.to,mergeOrders.from,'exact'))
		mergeOrders.from(strmatch(mergeOrders.to,mergeOrders.from,'exact')) = [];
	end
	setappdata(fig,'mergeOrders',mergeOrders);
elseif length(selectedRois) > 1
	mergeOrders(end+1).to = selectedRois{tmpi(1)};
	if isempty(selectedRois(tmpi(2:end)))
		mergeOrders(end) = [];
		setappdata(fig,'mergeOrders',mergeOrders);
	else
		mergeOrders(end).from = unique(selectedRois(tmpi(2:end)));
		if ~isempty(strmatch(mergeOrders(end).to,mergeOrders(end).from,'exact'))
			mergeOrders(end).from(strmatch(mergeOrders(end).to,mergeOrders(end).from,'exact')) = [];
		end
		setappdata(fig,'mergeOrders',mergeOrders);
    end
else
    return
end

setappdata(fig,'curRoi',mergeOrders(end).to);
setappdata(fig,'calcNeighbors',true);
updateCurRoiList(fig);
updateFrame(fig);
setappdata(fig,'needsSave',true);


% --- Executes on button press in compareBtn.
function compareBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
tags = getappdata(fig,'tags');
inten = getappdata(fig,'inten');
roi = getappdata(fig,'roi');
stimData = getappdata(fig,'stimData');
curRoi = getappdata(fig,'curRoi');

ops = struct;
ops.roitags = {curRoi};
ops.step = false;
f1 = intensity_stepper(inten,roi,ops);
ops.roitags = getappdata(fig,'selectedRois');
if ~isempty(ops.roitags{1})
	f2 = intensity_stepper(inten,roi,ops);
	pos = get(f1,'position'); pos(1) = pos(1)+601;
	set(f2,'position', pos);
end

% --- Executes on button press in selCurrentBtn.
function selCurrentBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
tags = getappdata(fig,'tags');
roi = getappdata(fig,'roi');
handles = guihandles(fig);
val = get(handles.selRoiList,'value');
str = get(handles.selRoiList,'string');
if iscell(str)
	this = str{val};
else
	this = str;
end
if isempty(this)
	return;
end
[ismerged groupi] = isMerged(fig,this);
if ismerged
	tg = roi.roi_defs(groupi(1)).label;
else
	tg = this;
end
setappdata(fig,'curRoi',tg);
setappdata(fig,'calcNeighbors',true);
id2 = strmatch(tg,get(handles.curRoiList,'string'),'exact');
set(handles.curRoiList,'value',id2);
updateFrame(fig);

% --- Executes on button press in selCenter.
function selCenter_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
handles = guihandles(fig);
val = get(handles.selRoiList,'value');
str = get(handles.selRoiList,'string');
if iscell(str)
	if isempty(str{1})
		return;
	end
	this = str{val};
else
	this = str;
end
findCtr(fig,this);

% --- Executes on button press in SaveBtn.
function SaveBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
roi = getappdata(fig,'roi');
smm = getappdata(fig,'smm');
header = roi.header;
pixel_spacing = header.pixel_spacing;
sz = smm.size;
mergeOrders = getappdata(fig,'mergeOrders');

output_file = uiputfile('*.roidef','Please select an output file for the new roi definitions');
if isempty(output_file)
	return;
end

prog  = progress(struct('progress',0,'max',length(roi.roi_defs)-length([mergeOrders.from])));
prog = progress(prog);

pause(0.1);

count = 1;
mergedfrom = [mergeOrders.from];
mergedto = {mergeOrders.to};
roicp = roi.roi_defs;
im2flow = getappdata(fig,'im2flow');
if ~isempty(im2flow)
		weight_metric = questdlg('Which metric should be used for calculating weights?','Choose weight metric','zscore','dfof',weight_metric);
		img = array_prolong(im2flow.(weight_metric),[sz(1:3) size(im2flow.(weight_metric),4)]);
else
	displayfile = getappdata(fig,'displayfile');
	if exist(displayfile,'file')
		tmpim = load(displayfile,'-mat');
		weight_metric = tmpim.options.weight_metric;
		weight_metric = questdlg('Which metric should be used for calculating weights?','Choose weight metric','zscore','dfof',weight_metric);
		img = array_prolong(tmpim.im2flow.(weight_metric),[sz(1:3) size(tmpim.im2flow.(weight_metric),4)]); clear tmpim;
	else
		displayfile = uigetfile('*roi_by_imflow_temp*', 'Choose a "display" file for assigning weights');
		tmpim = load(displayfile,'-mat');
		weight_metric = tmpim.options.weight_metric;
		weight_metric = questdlg('Which metric should be used for calculating weights?','Choose weight metric','zscore','dfof',weight_metric);
		img = array_prolong(tmpim.im2flow.(weight_metric),[sz(1:3) size(tmpim.im2flow.(weight_metric),4)]); clear tmpim;
	end
end
roi_defs = roicp; roi_defs(:) = [];
for i = 1:length(roicp)
	lbl = roicp(i).label;
	if ~isstr(lbl)
		lbl = num2str(lbl);
	end
	isto = strmatch(lbl,mergedto,'exact');
	isfrom = strmatch(lbl,mergedfrom,'exact');
	if isto
		[~,idx]=isMerged(fig,lbl);
		roi_defs(count).pixels = [];
		for j = idx
			roi_defs(count).pixels = [roi_defs(count).pixels; roicp(j).pixels];
		end
		if isnumeric(roicp(count).label)
			roi_defs(count).label = count; % if integer, just assign new integer
		else
			roi_defs(count).label = roicp(i).label; % if a hand-made label, keep the "to" label
		end
		
		roi_defs(count).centerInPixels = round(mean(roi_defs(count).pixels,1));
		roi_defs(count).centerInUm = roi_defs(count).centerInPixels.*pixel_spacing;
		tmp = zeros(sz(1:3));
		tmp(sub2ind(sz(1:3),roi_defs(count).pixels(:,1),roi_defs(count).pixels(:,2),roi_defs(count).pixels(:,3))) = 1;
		perim = bwperim(tmp,4);
		vtx = [];
		[vtx(:,1), vtx(:,2) vtx(:,3)] = ind2sub(sz(1:3), find(perim));
		[s, si] = sort(vtx(:,3),'ascend');
		vtx = vtx(si,:);
		for i = 1:sz(3)
			roi_defs(count).vtxInPixels{i} = [];
			these = find(s==i);
			if ~isempty(these)
				thesevtx = vtx(these,:);
				while ~isempty(thesevtx)
					this = thesevtx(1,:);
					roi_defs(count).vtxInPixels{i} = [roi_defs(count).vtxInPixels{i}; this];
					thesevtx(1,:) = [];
					dist = sqrdist(this', thesevtx');
					[~, sd] = sort(dist,'ascend');
					thesevtx = thesevtx(sd,:);
				end
			end
			if ~isempty(roi_defs(count).vtxInPixels{i})
				roi_defs(count).vtxInUm{i}(:,1) = roi_defs(count).vtxInPixels{i}(:,1)*pixel_spacing(1);
				roi_defs(count).vtxInUm{i}(:,2) = roi_defs(count).vtxInPixels{i}(:,2)*pixel_spacing(2);
				roi_defs(count).vtxInUm{i}(:,3) = roi_defs(count).vtxInPixels{i}(:,3)*pixel_spacing(3);
			else
				roi_defs(count).vtxInUm{i} = [];
			end
		end
		roi_defs(count).roiVolumeInUm3 = size(roi_defs(count).pixels,1)*prod(pixel_spacing);
		if exist(displayfile,'file')
			snip = [];
			pixi = sub2ind(sz(1:3),roi_defs(count).pixels(:,1),roi_defs(count).pixels(:,2),roi_defs(count).pixels(:,3));
			for j = 1:size(img,4)
				cut = img(:,:,:,j);
				snip(:,j) = cut(pixi);
			end
			rmean = nanmean(snip,1);
			c = sum(bsxfun(@times,snip,rmean),2);
			c(c<0)=0;c(isnan(c)) = 0;
			roi_defs(count).weight = c/sum(c);
		end
		count = count+1;
	elseif isfrom
		continue;
	else
		roi_defs(count) = roicp(i);
		if isnumeric(roicp(i).label)
			roi_defs(count).label = count;
		end
		count = count+1;
    end
    prog.progress=count-1;
    progress(prog);
end

pixelPerUm = roi.pixelPerUm;
save(output_file,'roi_defs','header','pixelPerUm');

prog.progress = -1;
prog = progress(prog); if ishandle(prog.handle); delete(prog.handle); end; clear prog;

setappdata(fig,'needsSave',false);

% --- Custom functions below ---
function updateCurRoiList(fig)
roi = getappdata(fig,'roi');
handles = guihandles(fig);

if ~getappdata(fig,'init_done')
	lst = {roi.roi_defs.label};
	set(handles.curRoiList,'string',lst);
	set(handles.curRoiList,'value',1);
	setappdata(fig,'tags',lst);
	setappdata(fig,'curList',lst);
else
    [ismerged,groupi] = isMerged(fig);
    if ~ismerged
        groupi = getIdx(fig,getappdata(fig,'curRoi'));
    end
	tagcp = getappdata(fig,'tags');
	mergeOrders = getappdata(fig,'mergeOrders');
	merged = [mergeOrders.from];
	found = [];
	for i = 1:length(merged)
		found = [found strmatch(merged{i},tagcp,'exact')];
	end
	tagcp(found) = [];
	[~,st] = sort(groupi);
	curRoi = roi.roi_defs(groupi(st(1))).label;
	setappdata(fig,'curRoi',curRoi);
	set(handles.curRoiList,'string',tagcp);
	setappdata(fig,'curList',tagcp);
	if isnumeric(curRoi)
		lbli = strmatch(num2str(curRoi),get(handles.curRoiList,'string'),'exact');
	else
		lbli = strmatch(curRoi,get(handles.curRoiList,'string'),'exact');
	end
	set(handles.curRoiList,'value',lbli);
end


% --- helper function to set SelRoiList string & value
function updateSelRoiList(fig)

handles = guihandles(fig);
selectedRois = getappdata(fig,'selectedRois');
set(handles.selRoiList,'string',selectedRois);
set(handles.selRoiList,'value',1);


% --- Executes on button press in deleteCurBtn.
function deleteCurBtn_Callback(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
roi = getappdata(fig,'roi');
curRoi_val = get(handles.curRoiList,'value');
tmp = get(handles.curRoiList,'string');
mergeOrders = getappdata(fig,'mergeOrders');
if isempty(curRoi_val)
    curRoi_str = [];
else
    curRoi_str = tmp{curRoi_val};
end

[ismerged groupi mergeIndex] = isMerged(fig,curRoi_str);

if ~isempty(curRoi_str)
	if ismerged
		mergeOrders(mergeIndex) = [];
		if length(mergeOrders) == 1 && isempty(mergeOrders(1).from)
			mergeOrders = struct;
		end
		for i = 1:length(groupi)
			thislabel = roi.roi_defs(groupi(i)).label;
			if isempty(fieldnames(mergeOrders(1)))
				mergeOrders(1).to = '';
				mergeOrders(1).from = {thislabel};
			else
				mergeOrders(end+1).to = '';
				mergeOrders(end).from = {thislabel};
			end
		end
	else
		if isempty(fieldnames(mergeOrders(1)))
			mergeOrders(1).to = '';
			mergeOrders(1).from = {curRoi_str};
		else
			mergeOrders(end+1).to = '';
			mergeOrders(end).from = {curRoi_str};
		end
	end
	setappdata(fig,'mergeOrders',mergeOrders);
else
	return;
end

curList = getappdata(fig,'curList');
found = strmatch(curRoi_str,curList);
if ~isempty(found)
	curList(found) = [];
end
setappdata(fig,'curList',curList);
setappdata(fig,'needsSave',true);
roiIdx = getIdx(fig,curRoi_str);
if roiIdx-1 > 0 && ~isempty(strmatch(roi.roi_defs(roiIdx-1).label,curList))
	setappdata(fig,'curRoi',roi.roi_defs(roiIdx-1).label);
elseif roiIdx+1 < length(roi.roi_defs) && ~isempty(strmatch(roi.roi_defs(roiIdx+1).label,curList))
	setappdata(fig,'curRoi',roi.roi_defs(roiIdx+1).label);
else
	found = [];
	rois = 1:length(roi.roi_defs); rois(roiIdx) = [];
	while isempty(found)
		[~,roiIdx] = mindist(roiIdx,rois); roiIdx = rois(roiIdx);
		found = strmatch(roi.roi_defs(roiIdx).label,curList);
		if isempty(found)
			rois(roiIdx) = [];
		end
	end
	setappdata(fig,'curRoi',roi.roi_defs(roiIdx).label);
end
setappdata(fig,'calcNeighbors',true);
updateCurRoiList(fig);
updateFrame(fig);
findCtr(fig);


% --- Executes on key press with focus on curRoiList and none of its controls.
function curRoiList_KeyPressFcn(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
switch eventdata.Character
	case '+'
		roi = getappdata(fig,'roi');
		curRoi = getappdata(fig,'curRoi');
        [ismerged, groupi] = isMerged(fig);
        if ismerged
            roiIdx = groupi;
        else
            roiIdx = getIdx(fig,curRoi);
        end
% 		roiIdx = getIdx(fig,curGroup);
		curLbl = {roi.roi_defs(roiIdx).label}; 
		for i = 1:length(curLbl)
			curLbl{i} = num2str(curLbl{i});
		end
		selectedRois = getappdata(fig,'selectedRois');
		if length(curLbl)<2
			if isempty(strmatch(curLbl,selectedRois,'exact'))
				if isempty(selectedRois{1})
					selectedRois = curLbl;
				else
					selectedRois = [selectedRois curLbl];
				end
			end
		elseif length(selectedRois) < 2
			if isempty(strmatch(selectedRois,curLbl,'exact'))
				if isempty(selectedRois{1})
					selectedRois = curLbl;
				else
					selectedRois = [selectedRois curLbl];
				end
			end
		else
			sm = length(curLbl) < length(selectedRois);
			if sm
				overlap = indexainb(curLbl,selectedRois);
				if ~isempty(overlap)
					selectedRois(overlap) = [];
				else
					selectedRois = [selectedRois curLbl];
				end
			else
				overlap =  indexainb(selectedRois,curLbl);
				if ~isempty(overlap)
					selectedRois(overlap) = [];
				else
					selectedRois = [selectedRois curLbl];
				end
			end
		end
		set(handles.selRoiList,'string',selectedRois);
		set(handles.selRoiList,'value',length(selectedRois));
		setappdata(fig,'selectedRois',selectedRois);
		setappdata(fig,'calcNeighbors',true);
		updateFrame(fig);
	otherwise
		return;
end


% --- Executes on key press with focus on selRoiList and none of its controls.
function selRoiList_KeyPressFcn(hObject, eventdata, handles)

fig = get_parent_fig(hObject);
if ~isempty(strmatch(eventdata.Key,'delete'))
	eventdata.Character = '-';
end
switch eventdata.Character
	case '-'
		roi = getappdata(fig,'roi');
		str = get(hObject,'string');
		selected_roi = get(hObject,'value');
		if ~isempty(str)
			curLbl = str{selected_roi};
		else
			return
		end
		if ~isempty(curLbl)
			[ismerged,curSelGrp] = isMerged(fig,curLbl);
			selectedRois = getappdata(fig,'selectedRois');
			if ~isempty(strmatch(curLbl,selectedRois,'exact'))
				if length(selectedRois)<2
					selectedRois = {''};
				else
					selectedRois(strmatch(curLbl,selectedRois,'exact')) = [];
				end
				set(handles.selRoiList,'value',length(selectedRois));
				set(handles.selRoiList,'string',selectedRois);
			end
			setappdata(fig,'selectedRois',selectedRois);
			setappdata(fig,'calcNeighbors',true);
			updateFrame(fig);
		end
end

% --- Executes on button press in pcaSplitBtn.
function pcaSplitBtn_Callback(hObject, eventdata, handles)
fig = get_parent_fig(hObject);
% im2flow = getappdata(fig,'im2flow');
% selectedRois = getappdata(fig,'selectedRois');
% roi = getappdata(fig,'roi');
% smm = getappdata(fig,'smm');
% smmh = smm.header;
% smm_size = smm.size;
% 
% if isempty(selectedRois{1})
% 	return;
% end
% default_display = questdlg('Which metric do you want to use as the basis for PCA calculations?','Choose metric for PCA', 'dfof','zscore','dfof');
% idx = getIdx(fig,selectedRois);
% 
% if any(smm_size(1:3)>size(im2flow.(default_display)(:,:,:,1)))
% 	temp = zeros(smm_size(1:3));
% 	thesepix = [roi.roi_defs(idx).pixels];
% 	thesepix = sub2ind(smm_size(1:3),thesepix(:,1),thesepix(:,2),thesepix(:,3));
% 	thesepix = sort(thesepix);
% 	temp(thesepix)=1;
% 	fromsz = smm_size(1:3); tosz = size(im2flow.(default_display)(:,:,:,1));
% 	red = [0 0 0];
% 	while any(fromsz>tosz)
% 		red = red + fromsz>tosz;
% 		fromsz = fromsz./((red>0)+1);
% 	end
% 	temp = array_restrict(temp,red);
% 	thesepix = find(temp>0.5);
% % 	temp(thesepix) = 1;
% else
% 	thesepix = sort([roi.roi_defs(idx).pixels]);
% 	thesepix = sub2ind(size(im2flow.(default_display)(:,:,:,1)),thesepix(:,1),thesepix(:,2),thesepix(:,3));
% end
% 
% for i = 1:size(im2flow.(default_display),4)
% 	temp = im2flow.(default_display)(:,:,:,i);
% 	pixvals(:,i) = temp(thesepix);
% end
% [pcomps sqrt_lambda]= pca(pixvals);
% eigen_sum = sum(sqrt_lambda);
% eigen_normd = sqrt_lambda/eigen_sum;
% for i = 1:length(eigen_normd)-1
% 	eigen_ratio(i) = eigen_normd(i)/(eigen_normd(i)-eigen_normd(i+1));
% end	
% 
% % set up pop-up fig
% units = get(fig,'units');
% set(fig,'units','pixels');
% curpos = get(fig,'position');
% set(fig,'units',units);
% htem = figure('units','pixels','name','Split GUI','position',[curpos(1)-500 curpos(2) 600 400]);
% ax = axes('parent',htem);



% --- this is the core rendering/highlighting function ---
function updateFrame(fig)
handles = guihandles(fig);
img = getappdata(fig,'img');
curFrame = getappdata(fig,'curFrame');
curRoi = getappdata(fig,'curRoi');
selectedRois = getappdata(fig,'selectedRois');
mergeOrders = getappdata(fig,'mergeOrders');
tags = getappdata(fig,'tags');
roi = getappdata(fig,'roi');
handles.imgAxis = getappdata(fig,'imgHandle');

% -- if needing to initialize --
if ~getappdata(fig,'init_done')
	imagesc(img(:,:,curFrame),'parent',handles.imgAxis);
	set(handles.imgAxis,'visible','off','xtick',[],'ytick',[],'color',[0 0 0]);
	cmap = colorize_asymmetric([min(img(:)) max(img(:))]);
	cmap(1,:) = 0;
	colormap(handles.imgAxis,cmap);
	roiIdx = getIdx(fig,curRoi);
	roi_vtx = roi.roi_defs(roiIdx).vtxInPixels{curFrame}(:,1:2);
	roi_vtx = rot90(roi_vtx,2);  % NOTE: the vertices must be inverted to match image!!!
	hpatch = patch(roi_vtx(:,1), roi_vtx(:,2), [0.0 0.4 0.0],'parent',handles.imgAxis);
	set(hpatch,'buttondownfcn', @patchFocus);
	setappdata(hpatch,'thisRoi',curRoi);
	set(hpatch,'hittest','on');
	if length(roi_vtx) > 50
		set(hpatch,'linestyle','none','facealpha',0.3);
	end
	setappdata(fig,'curPatch',hpatch);
	findNeighbors(fig);
	return;
end

% -- section for moving between frames --
delete(get(handles.imgAxis,'children'));
imagesc(img(:,:,curFrame),'parent',handles.imgAxis);
set(handles.imgAxis,'visible','off','tag','imgAxis');
set(handles.curFrame_edit,'string',num2str(curFrame));

% -- section for putting a patch on "current" rois --
[ismerged, groupi] = isMerged(fig);
if ismerged > 0
	if ismerged > 1
		error('there are multiple merge orders for the same ROI. this shouldn''t have happened.  fix the GUI!');
	end
	roi_vtx = cell(1);
	maxsz = 0;
	for i = 1:length(groupi)
		tmp = roi.roi_defs(groupi(i)).vtxInPixels{curFrame};
		if ~isempty(tmp)
			tmp = rot90(tmp(:,1:2),2);
			maxsz = max([maxsz size(tmp,1)]); 
			roi_vtx{i} = tmp;
		end
	end
	if maxsz > 0
		tmp = NaNs([maxsz,2,i]);
		for i = 1:length(roi_vtx);
			thissz = size(roi_vtx{i},1);
			if thissz > 0
				tmp(1:thissz,:,i) = roi_vtx{i};
				tmp(thissz:end,1,i) = tmp(thissz,1,i);
				tmp(thissz:end,2,i) = tmp(thissz,2,i);
			end
		end
		roi_vtx = tmp;
	else
		roi_vtx = 0;
	end

	if nansum(roi_vtx(:))~=0
		patchMyRoi(fig,roi_vtx,[0.0 0.9 0.0], 0.3,roi.roi_defs(groupi(1)).label);
	else
		setappdata(fig,'curPatch',[]);
    end
else
	roiIdx = getIdx(fig,curRoi);
	roi_vtx = roi.roi_defs(roiIdx).vtxInPixels{curFrame};
	if ~isempty(roi_vtx)
		roi_vtx = roi_vtx(:,1:2);
		roi_vtx = rot90(roi_vtx,2);
		cpatch = patchMyRoi(fig,roi_vtx,[0.0 0.9 0.0], 0.3,roiIdx);
		setappdata(fig,'curPatch',cpatch);
	else
		setappdata(fig,'curPatch',[]);
    end
end

% -- call findNeighbors --
% -- avoid this step (slow) if no changes have been made to ROI selection
if getappdata(fig,'calcNeighbors')
    findNeighbors(fig);
    setappdata(fig,'calcNeighbors',false);
end
    
% -- section for putting a patch on "selected" rois --
if ~isempty(selectedRois{1})
	selIndex = [];
	for i = 1:length(selectedRois)
		selIndex = [selIndex strmatch(selectedRois{i},tags,'exact')];
	end

	roi_vtx = cell(1);
	maxsz = 0;
	for i = 1:length(selIndex)
		tmp = roi.roi_defs(selIndex(i)).vtxInPixels{curFrame};
		if ~isempty(tmp)
			tmp = rot90(tmp(:,1:2),2);
			maxsz = max([maxsz size(tmp,1)]); 
			roi_vtx{i} = tmp;
		end
	end
	if maxsz > 0
		tmp = NaNs([maxsz,2,i]);
		for i = 1:length(roi_vtx);
			thissz = size(roi_vtx{i},1);
			if thissz > 0
				tmp(1:thissz,:,i) = roi_vtx{i};
				tmp(thissz:end,1,i) = tmp(thissz,1,i);
				tmp(thissz:end,2,i) = tmp(thissz,2,i);
			end
		end
		roi_vtx = tmp;
	else
		roi_vtx = 0;
	end

	if nansum(roi_vtx(:))~=0
		patchMyRoi(fig,roi_vtx,[0.4 0.4 0.4], 0.3,selIndex(1));
	else
		setappdata(fig,'selPatch',[]);
	end
else
	setappdata(fig,'selPatch',[]);
end


% -- put a patch down with the specified color, alpha)
function hpatch = patchMyRoi(fig,roi_vtx,color,alpha,roi)

hpatch = patch(squeeze(roi_vtx(:,1,:)), squeeze(roi_vtx(:,2,:)), color,'parent',getappdata(fig,'imgHandle'));
if length(roi_vtx) > 20
	set(hpatch,'linestyle','none','facealpha',alpha);
else
	set(hpatch,'edgecolor',color, 'linewidth',4);
end
setappdata(hpatch,'thisRoi',roi);
set(hpatch,'buttondownfcn',@patchFocus);
set(hpatch,'hittest','on');


% --- returns the indices of ROIs set to be merged with the supplied ROI
function [ismerged, groupi, mergeIndex] = isMerged(fig,varargin)
 handles = guihandles(fig);
tags = getappdata(fig,'tags');% get(handles.curRoiList,'string');
roi = getappdata(fig,'roi');
mergeOrders = getappdata(fig,'mergeOrders');
if isempty(varargin)
	curRoi = getappdata(fig,'curRoi');
elseif ischar(varargin{1})
	curRoi = varargin{1};
end

roiIdx = getIdx(fig,curRoi);

if length(fieldnames(mergeOrders))<1
	ismerged = 0;
else
	ismerged = sum(find(strmatch(num2str(roi.roi_defs(roiIdx).label),[mergeOrders.to mergeOrders.from],'exact')));
end

if ismerged < 1
	groupi = NaN;
	mergeIndex = NaN;
	return;
end
for i = 1:length(mergeOrders)
	these = [mergeOrders(i).to mergeOrders(i).from];
	foundit = ~isempty(strmatch(num2str(roi.roi_defs(roiIdx).label),these,'exact'));
	if foundit
		mergeIndex = i;
		for j = 1:length(these)
			groupi(j) = strmatch(these{j},tags,'exact'); %convert to index
		end
		break;
	end
end


% --- the main function for assigning "closest neighbors" ---
function findNeighbors(fig)
handles = guihandles(fig);
roi = getappdata(fig,'roi');
pixel_spacing = roi.header.pixel_spacing;
inten = getappdata(fig,'inten');
mergeOrders = getappdata(fig,'mergeOrders');
curRoi = getappdata(fig,'curRoi');
selectedRois = getappdata(fig,'selectedRois');
stimData = getappdata(fig,'stimData');
	
[ismerged groupi] = isMerged(fig,curRoi);
if ismerged
	roiIdx = groupi;
else
	roiIdx = getIdx(fig,curRoi);
end

if ~isempty(fieldnames(mergeOrders))
	merged_or_deleted = [mergeOrders.from];
	mdIdx = getIdx(fig,merged_or_deleted);
else
	mdIdx = [];
end

if ~isempty(getappdata(fig,'distmatrix'))
	distmatrix = getappdata(fig,'distmatrix');
	cutoff = getappdata(fig,'cutoff');
	if ~isempty(fieldnames(mergeOrders))
		for i = 1:length(mergeOrders)
			to = mergeOrders(i).to;
			from = mergeOrders(i).from;
			if isempty(to)
				distmatrix(getIdx(fig,from),:) = NaN;
				distmatrix(:,getIdx(fig,from)) = NaN;
			else
				tofrom = [{to} from];
				groupmin = min(distmatrix(getIdx(fig,tofrom),roiIdx));
				distmatrix(getIdx(fig,to),roiIdx) = groupmin;
				distmatrix(roiIdx,getIdx(fig,to)) = groupmin;
				distmatrix(getIdx(fig,from),roiIdx) = NaN;
				distmatrix(roiIdx,getIdx(fig,from)) = NaN;
			end
		end
	end
	curRoiDists = min(distmatrix(roiIdx,:),[],1);
	curRoiDists(curRoiDists==0)=NaN; % self-pointing
	[closest, i_closest] = sort(curRoiDists);
	i_closest(isnan(closest))=[];
	closest(isnan(closest))=[];	
	i_closest(closest>cutoff) = [];
	closest(closest>cutoff)=[];
else
	% HACK - hard-code max & of landmarks
	nland = 5e2;
	% find physically-closest ROIs (top 10% ?)
	temp = cat(1,roi.roi_defs(roiIdx).pixels);
	ncur = size(temp,1);
	cur_pix = temp(randsample(1:ncur,min([nland max([ceil(ncur/10) min([nland ncur])])])),:);  % subsample for speed  HACK hard-coded
	cpix(:,1) = cur_pix(:,1)*pixel_spacing(1);
	cpix(:,2) = cur_pix(:,2)*pixel_spacing(2);
	cpix(:,3) = cur_pix(:,3)*pixel_spacing(3);
	clear temp;
	tocheck = 1:length(roi.roi_defs);
	tocheck(findainb(unique([roiIdx mdIdx]),tocheck))=[];

	minPhys = NaNs([1 length(roi.roi_defs)]);
	%tic;
	for i = tocheck
		[ismerged groupi] = isMerged(fig,roi.roi_defs(i).label);
		if ismerged
			thesei = groupi;
		else
			thesei = i;
		end
		thesepix = cat(1,roi.roi_defs(thesei).pixels);
		thislen = size(thesepix,1);
		thesepix = thesepix(randsample(1:thislen,min([nland max([ceil(thislen/10) min([nland thislen])])])),:); %HACK hard-coded
		tpix = [];
		tpix(:,1) = thesepix(:,1)*pixel_spacing(1);
		tpix(:,2) = thesepix(:,2)*pixel_spacing(2);
		tpix(:,3) = thesepix(:,3)*pixel_spacing(3);
		minPhys(i) = min(pdist2(cpix,tpix,'Euclidean','Smallest',1));
	end
	%toc;
	
	[closest,i_closest] = sort(minPhys,'ascend');
	i_closest = i_closest(~isnan(closest));
	closest = closest(~isnan(closest));
	thresh = 50;
	i_closest = i_closest(closest<thresh); % HACK hard-coded threshold!!!
	closest = closest(closest<thresh); % HACK hard-coded threshold!!!
end

this = 0;
for i = 1:length(roiIdx)
	this = this+squeeze(stimData.activity(roiIdx(i),:,:))*stimData.volume(roiIdx(i))/sum(stimData.volume(roiIdx));
end


covM = [];
tmp = [];

% for macroscopic comparisons, compute the across-all-repeats-and-stimuli covariance
if ~isempty(i_closest)
	act = cell(1);
	for i = 1:length(i_closest)
		act{i} = squeeze(stimData.activity(i_closest(i),:,:));
		tmp = cov(this,act{i});
		covM(i) = tmp(1,2)/mean([tmp(1,1) tmp(2,2)]);
	end
	act{end+1} = squeeze(mean(stimData.activity(:,:,:),1));
	tmp = cov(this, act{end});
	covM(end+1) = tmp(1,2)/mean([tmp(1,1) tmp(2,2)]);

	[~, covsorti] = sort(covM(1:length(i_closest)),'descend');
	closest = closest(covsorti);
	setappdata(fig,'closest',closest);
	if length(i_closest)>3
		top = i_closest(covsorti(1:3));  %% HACK CURRENTLY HARD-CODED to 3.  COULD CHANGE THIS!!! %%
	else
		top = i_closest(covsorti);
	end

	for i = 1:3
		if i <= length(top)
			setappdata(fig,['ldAx' num2str(i) 'Plot'],[this(:),act{covsorti(i)}(:)]);
			setappdata(fig,['ldAx' num2str(i) 'CorrCoef'], covM(covsorti(i)));
			setappdata(fig,['ldAx' num2str(i) 'Idx'], top(i));
		else
			setappdata(fig,['ldAx' num2str(i) 'Plot'],[]);
			setappdata(fig,['ldAx' num2str(i) 'CorrCoef'], []);
			setappdata(fig,['ldAx' num2str(i) 'Idx'], []);
		end
	end
	if i<=length(top)
		setappdata(fig,'ldAxLastPlot',[this(:),act{end}(:)]);
		setappdata(fig,'ldAxLastCorrCoef',covM(end));
	else
		setappdata(fig,'ldAxLastPlot',[this(:),act{end}(:)]);
		setappdata(fig,'ldAxLastCorrCoef',covM(end));
	end
else
	for i = 1:3
		setappdata(fig,['ldAx' num2str(i) 'Plot'],[]);
		setappdata(fig,['ldAx' num2str(i) 'CorrCoef'], []);
		setappdata(fig,['ldAx' num2str(i) 'Idx'], []);
	end
	act{1} = squeeze(mean(stimData.activity(:,:,:),1));
	tmp = cov(this, act{end});
	covM = tmp(1,2)/mean([tmp(1,1) tmp(2,2)]);
	setappdata(fig,'ldAxLastPlot',[this(:),act{end}(:)]);
	setappdata(fig,'ldAxLastCorrCoef',covM(end));
end

updateClosest(fig);


% --- patchFocus
function patchFocus(hObject,eventdata)
fig = get_parent_fig(hObject);
thisRoi = getappdata(hObject,'thisRoi');
if isempty(thisRoi)
	return;
end
roi = getappdata(fig,'roi');
mergeOrders = getappdata(fig,'mergeOrders');
handles = guihandles(fig);
if isnumeric(roi.roi_defs(thisRoi).label)
	thisstr = num2str(roi.roi_defs(thisRoi).label);
else
	thisstr = roi.roi_defs(thisRoi).label;
end
found = strmatch(thisstr,get(handles.curRoiList,'string'),'exact');
if ~isempty(found)
	set(handles.curRoiList,'value',found);
	curRoiList_Callback(handles.curRoiList,[],[]);
else
	for i = 1:length(mergeOrders)
		found = ~isempty(strmatch(thisstr,mergeOrders(i).from,'exact'));
		if found
			found = mergeOrders(i).to;
			found = strmatch(found,get(handles.curRoiList,'string'),'exact');
			set(handles.curRoiList,'value',found);
			curRoiList_Callback(handles.curRoiList,[],[]);
		end
	end
end


% --- function to update the comparison plots along the bottom ---
function updateClosest(fig)
strng = 'Covariance: ';  % HACK I know, hard-coding is terrible
handles = guihandles(fig);
closest = getappdata(fig,'closest');
roi = getappdata(fig,'roi');
for i = 1:3 % HARD CODED to 3... CAN CHANGE IN FUTURE!
	thisPlot = getappdata(fig,['ldAx' num2str(i) 'Plot']);
	corrCoef = getappdata(fig,['ldAx' num2str(i) 'CorrCoef']);
	thisIdx = getappdata(fig,['ldAx' num2str(i) 'Idx']);
	if ~isempty(thisPlot)
		p = plot(handles.(['ldAx' num2str(i)]), thisPlot(:,1), thisPlot(:,2),'k.');
		set(p,'hittest','off');
	else
		delete(get(handles.(['ldAx' num2str(i)]),'children'));
	end
	set(handles.(['ldAx' num2str(i)]),'tag',['ldAx' num2str(i)'],'box','off','xtick',[],'ytick',[],...
		'ButtonDownFcn', @(hObject,eventdata)roi_merge(['ldAx' num2str(i) '_ButtonDownFcn'],hObject,eventdata,guidata(hObject)));  % FOUND A BUG IN MATLAB PLOT FXN THAT REMOVES TAG

	%	axes(handles.(['ldAx' num2str(i)])); axis equal;
	if ~isempty(thisPlot)
		set(handles.(['nearRoi' num2str(i) '_txt']),'string',['ROI ' num2str(roi.roi_defs(thisIdx).label) ' (' num2str(ceil(closest(i))), 'um)']);
		set(handles.(['nearRoi' num2str(i) '_coef']),'string',[strng sprintf('%0.3f',corrCoef)]);
	else
		set(handles.(['nearRoi' num2str(i) '_txt']),'string',[]);
		set(handles.(['nearRoi' num2str(i) '_coef']),'string',[]);
	end
end
thisPlot = getappdata(fig,'ldAxLastPlot');
plot(handles.ldAxOther, thisPlot(:,1),thisPlot(:,2),'k.');
set(handles.ldAxOther,'tag','ldAxOther','box','off','xtick',[],'ytick',[]);
%axes(handles.ldAxOther); axis equal;
set(handles.nearRoiOther_coef,'string',[strng num2str(getappdata(fig,'ldAxLastCorrCoef'))]);
guidata(fig,handles);

% --- function to calculate the max responses per stimulus, per trial ---
function calcDistances(fig)
roi = getappdata(fig,'roi');
roisz = length(roi.roi_defs);

input = questdlg('Would you like to load saved distances?','Load existing distances?','Yes','No','Yes');
if ~isempty(strmatch(input,'Yes','exact'))
	moveon = 0;
	while ~moveon
		filename = uigetfile('*.roidist','Select file to load distances');
		if ~filename
			break;
		end
		load(filename,'-mat')
		if ~exist('distmatrix','var') || ~exist('cutoff','var')
			input = questdlg('This file lacked distmatrix or cutoff information, try again?','Bad file, Try again?','Yes','No','Yes');
			if isempty(strmatch(input,'Yes','exact'))
				moveon = 1;
			end
		elseif length(distmatrix) ~= roisz
			input = questdlg('This file had distmatrix information of the wrong size, try again?','Wrong size, Try again?','Yes','No','Yes');
			if isempty(strmatch(input,'Yes','exact'))
				moveon = 1;
			end
		else
			setappdata(fig,'distmatrix',distmatrix);
			setappdata(fig,'cutoff',cutoff);
			return;
		end
	end
end

distmatrix = zeros(repmat(roisz,1,2));
pixel_spacing = roi.header.pixel_spacing;
N = length(roi.roi_defs);

% prog  = progress(struct('progress',0,'max',N));
% prog = progress(prog);
max_cutoff = 240;
min_cutoff = 120;
max_to_check = 5e2;
for i = 1:roisz
	centers(i,:) = roi.roi_defs(i).centerInUm;
end
% centers = reshape(centers,[roisz 3]);
center_dists = sqrt(sqrdist(centers',centers'));
base = min(center_dists(:)):1:max(center_dists(:));
nhist = histc(center_dists(:)',base);
fit = fit_gaussian(base,nhist);
cutoff = fit.center - fit.sigma;
cutoff = max([cutoff min_cutoff]);
cutoff = min([cutoff max_cutoff]);
dontcheck = center_dists>cutoff;
distmatrix(dontcheck) = center_dists(dontcheck);
tocheck = ~dontcheck;

progfig = figure('name','distance matrix: 0%');
him = axes('parent',progfig,'xlim',[-0.5 roisz+0.5], 'ylim',[-0.5 roisz+0.5]);
count = 0;
for i = 1:length(roi.roi_defs)
	checknow = find(tocheck(i,:) & ~distmatrix(i,:));
	count = count+1;
	pct = 100*count/roisz;
	set(progfig,'name',['distance matrix: ' sprintf('%0.1f',pct) '%']);
	for j = checknow
		if i == j
			if i == 1 && j == 1
				axes(him);
				imh = imagesc(distmatrix);
			else
				set(imh,'cdata',distmatrix);
				set(gca,'clim',[0 max([1 max(distmatrix(:))])]);
				drawnow;
			end
% 			prog.progress = prog.progress+1;
% 			prog = progress(prog);
			continue;
		end
% 		if distmatrix(i,j) || distmatrix(j,i)
% 			continue;
% 		end
% 		if norm(centers(i,:)-centers(j,:)) > min_center_spacing
% 			distmatrix(i,j) = min_center_spacing;
% 			distmatrix(j,i) = min_center_spacing;
% 		else
			a = []; b = [];
			if size(roi.roi_defs(i).pixels,1) > max_to_check
				samps = randsample(1:size(roi.roi_defs(i).pixels,1),max_to_check);
				a = roi.roi_defs(i).pixels(samps,:);
			else
				a = roi.roi_defs(i).pixels;
			end
			if size(roi.roi_defs(j).pixels,1) > max_to_check
				samps = randsample(1:size(roi.roi_defs(j).pixels,1),max_to_check);
				b = roi.roi_defs(j).pixels(samps,:);
			else
				b = roi.roi_defs(j).pixels;
			end
			distmatrix(i,j) = min(pdist2(a,b,'Euclidean','Smallest',1));
			distmatrix(j,i) = distmatrix(i,j);
% 		end
		set(imh,'cdata',distmatrix);
		set(gca,'clim',[0 max([1 max(distmatrix(:))])]);
		drawnow;
	end
end
setappdata(fig,'distmatrix',distmatrix);
setappdata(fig,'cutoff',cutoff);
input = questdlg('Would you like to save these distances?','Save distances?','Yes','No','Yes');
if ~isempty(strmatch(input,'Yes','exact'));
	filename = uiputfile('*.roidist','Select file to save distances');
	if isempty(strfind(filename,'.roidist'))
		filename = [filename '.roidist'];
	end
	save(filename,'distmatrix','cutoff','-mat');
end

% --- function to calculate the max responses per stimulus, per trial ---
function calcIntervals(fig)
smm = getappdata(fig,'smm');
roi = getappdata(fig,'roi');
inten = getappdata(fig,'inten');
stimData = struct;

[b,a] = cheby2(4,25,[1/1000 1/20]/(0.5/5));
tmp1 = [inten.intensities(:,end:-1:1) inten.intensities];
tmp = filter(b,a,tmp1,[],2);
inten.intensities = tmp(:,size(inten.intensities,2)+1:end);

smmh = smm.header;
stimData.onsets = find(diff(smmh.stim_lookup)>0)+1;
stimData.valvenums = unique(nonzeros(smmh.stim_lookup))';
stimData.ids = smmh.stim_lookup(stimData.onsets);
for i = 1:length(stimData.valvenums)
	stimData.onsetMatrix(i,:) = stimData.onsets(stimData.ids==stimData.valvenums(i));
end
stimData.activity = zeros([size(inten.intensities,1) size(stimData.onsetMatrix)]);

for i = 1:size(stimData.activity,1)
	base = 0;
	for j = -5:-1
		base = base+(1/5)*inten.intensities(i,stimData.onsetMatrix+j);
	end
	base = reshape(base,size(stimData.onsetMatrix));
	act = zeros([1 prod(size(stimData.onsetMatrix))]);
	for j = 0:4
		act = nanmean([act; inten.intensities(i,stimData.onsetMatrix+j)],1);
	end
	act = reshape(act,size(stimData.onsetMatrix));
	stimData.activity(i,:,:) = 100*(act-base)./base;
end
stimData.activity(:,stimData.ids(1),1) = 0; % first trial has the weirdness%
stimData.volume = [roi.roi_defs.roiVolumeInUm3];
stimData.labels = smmh.stim_labels;
setappdata(fig,'stimData',stimData);


% --- change display image dimensions to match 
function img = prepDisplayImg(fig,im2flow)
if ndims(im2flow) > 3
	img = squeeze(max(im2flow,[],4));
end
smm = getappdata(fig,'smm');
newsz = smm.size; newsz = newsz(1:3);
img = array_prolong(img,newsz);


% --- find the center frame for your ROI ---
function findCtr(fig,varargin)
roi = getappdata(fig,'roi');
if isempty(varargin)
	curRoi = getappdata(fig,'curRoi');
else
	curRoi = varargin{1};
end
roiIdx = getIdx(fig,curRoi);
newFrame = roi.roi_defs(roiIdx).centerInPixels(3);
setappdata(fig,'curFrame',newFrame);
updateFrame(fig);


% --- just grabs the index in the roi_defs structure matching a given ROI ---
function idx = getIdx(fig,lbl)
tags = getappdata(fig,'tags');
if ischar(lbl)
	idx = strmatch(lbl,tags,'exact');
elseif iscellstr(lbl)
	for i = 1:length(lbl)
		idx(i) = strmatch(lbl{i},tags,'exact');
	end
end
