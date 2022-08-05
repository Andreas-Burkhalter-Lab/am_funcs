function varargout = stack_roi_rg(varargin)
% STACK_ROI_RG: manually place ROIs, using red/green to code for activity
% Usage:
%   stack_roi_rg
% Prompts you to select the .imagine file via the GUI
%
%   stack_roi_rg(struct('filename','vno_gcamp2.imagine'))
% This version lets you supply a filename via the command line or a program
% 
%   stack_roi_rg(struct('label',roi_label,'display',display_popup_name))
% Assuming the GUI is already open and ROI definitions are loaded, this
% syntax allows you to "navigate" to see a particular ROI.
%
%   stack_roi_rg(struct('roiFilename',roiFilename,'label',roi_label,'display',display_popup_name))
% No need to load ROI definition before hand, but supply the 'roiFilename'
% when "navigate" to a particular ROI.
%
% % keybindings:
%    up/down arrow: change frame #
%    left/right arrow: change roi radius (currently disabled)
%    DEL/d: delete selected ROIs
%    right-drag on circle: change roi radius: 
%    left-drag on circle: move roi position: 
%    middle-click: define a new roi
%    right-click: look at clicked position in frames before and after (in
%                zoom panel)
%    left-click on center: select a roi
%    left-click on circle: deselect a roi
%    pageup/down: change stack #
%

%
% Written by Diwakar Turaga
% Adapted from stack_roi (written by Zhongsheng Guo)
% 2010_06_10
% Small modifications, Timothy E. Holy 2010-09-22
% Small modifications, update show_zoom_images after changing repeats/frames, etc. Pei S. Xu 2012-01-30
% Add subtraction and centralize the handling of pre/post timing, TEH 2012-03-03
% Partial but sizable cleanup of graphics updating and GUI state, TEH 2012-03-04
% small change in initialize_new to fix the bug for showing cluster
%   clickable to call STACK_ROI_RG, change "gui_Singleton" back to 1 for
%   navigating among different datasets without changing windows (may need
%   alternative to keep gui_Singleton == 0
%
% STACK_ROI_RG M-file for stack_roi_rg.fig
%      STACK_ROI_RG, by itself, creates a new STACK_ROI_RG or raises the existing
%      singleton*.
%
%      H = STACK_ROI_RG returns the handle to a new STACK_ROI_RG or the handle to
%      the existing singleton*.
%
%      STACK_ROI_RG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_ROI_RG.M with the given input arguments.
%
%      STACK_ROI_RG('Property','Value',...) creates a new STACK_ROI_RG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_roi_rg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_roi_rg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stack_roi_rg

% Last Modified by GUIDE v2.5 11-Apr-2012 12:59:50

% Begin initialization code - DO NOT EDIT
% if nargin ==2 && strcmp(varargin{2}, 'currentwindow')
% gui_Singleton = 1;
% else
%     gui_Singleton = 0;
% end

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @stack_roi_rg_OpeningFcn, ...
  'gui_OutputFcn',  @stack_roi_rg_OutputFcn, ...
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


% --- Executes just before stack_roi_rg is made visible.
function stack_roi_rg_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_roi_rg (see VARARGIN)

% Choose default command line output for stack_roi_rg
handles.output = hObject;

if isempty(varargin)
    % This is a new call with GUI to select the file
    handles = initialize_new(hObject,handles);
    [fileToOpen,pathName] = uigetfile('*.imagine; *.mat');
    if isequal(fileToOpen,0)
        return
    end
    fileToOpen = [pathName filesep fileToOpen];
    openFile(handles,fileToOpen)
else
    stackRoiParams=varargin{1};
    if ~isstruct(stackRoiParams)
        error('Usage: stack_roi(struct(''filename'', headerFile))');
    end
    if isfield(stackRoiParams,'filename')
        % open a new file
        fileToOpen=stackRoiParams.filename;
        handles = initialize_new(hObject,handles);
        openFile(handles, fileToOpen);
    else
        if isfield(stackRoiParams, 'roiFilename')
            roiFileToLoad = stackRoiParams.roiFilename;
            tt=load(roiFileToLoad, '-mat');
            % delete old rois and related objects
            roi_defs=getappdata(handles.stack_roi_rg, 'roi_defs');
            if(~isempty(roi_defs))
                set([roi_defs.hCircle], 'visible', 'on');
                onDeleteRoi(hObject, []);
            end
            
            roi_defs=tt.roi_defs;
            
            % remove fields
            roi_defs=rmfield(roi_defs, {'posInUm', 'xyradiusInUm'});
            
            idx = findainb(stackRoiParams.label,[roi_defs.label]);
            thisroi = roi_defs(idx);
            addRoi(handles,thisroi);
            
            % This should be a command-line call for an existing GUI; "navigate"
            % to the desired ROI
        else roi_defs = getappdata(handles.stack_roi_rg, 'roi_defs');
            if isempty(roi_defs)
                errordlg('Must first load ROI definitions','','modal');
                return
            end
        end
        idx = findainb(stackRoiParams.label,[roi_defs.label]);
        thisroi = roi_defs(idx);
        setappdata(handles.stack_roi_rg, 'curFrameIdx',thisroi.posInPixel(3));
        if isfield(stackRoiParams,'display')
            strs = get(handles.display_popupmenu,'String');
            strIdx = find(strcmp(stackRoiParams.display,strs));
            if ~isempty(strIdx)
                set(handles.display_popupmenu,'Value',strIdx);
            end
        end
        display_plot(handles.stack_roi_rg);
        show_zoom_images(handles, thisroi.posInPixel(1:2), thisroi.posInPixel(3));
    end
end

function handles = initialize_new(hObject,handles)
  % Set up blank images
  clim = [500 10000];  % default contrast limits
  axes_properties = {'CLim',clim,'Visible','off','YDir','normal'};
  image_properties = {'Visible','on','CDataMapping','scaled'};
  handles.main_image = image(0,'Parent',handles.image_axes,image_properties{:});
  set(handles.image_axes, axes_properties{:});
  % Set up the zoom axes in a more flexible way
  handles.zoom_axes = [handles.zoom1_axes handles.zoom2_axes handles.zoom3_axes handles.zoom4_axes handles.zoom5_axes];
  setappdata(handles.stack_roi_rg,'zoom_zoffsets',-2:2);
  % Set up blank images for zoom axes
  handles.zoom_image = handles.zoom_axes;
  for i = 1:length(handles.zoom_axes)
    handles.zoom_image(i) = image(0,'Parent',handles.zoom_axes(i),image_properties{:});
  end
  set(handles.zoom_axes,axes_properties{:});
  colormap(handles.image_axes,gray(256));

  fig=hObject;
  bind_shortcut(fig, 'uparrow', @onIncFrameIdx);
  bind_shortcut(fig, 'downarrow', @onDecFrameIdx);
  %bind_shortcut(fig, 'rightarrow', @onIncRoiRadius);
  %bind_shortcut(fig, 'leftarrow', @onDecRoiRadius);
  bind_shortcut(fig, 'd', @onDeleteRoi);
  bind_shortcut(fig, 'delete', @onDeleteRoi);
  bind_shortcut(fig, 'pageup', @onIncStackIdx);
  bind_shortcut(fig, 'pagedown', @onDecStackIdx);

  install_mouse_event_handler(handles.image_axes, 'up', @onMouseUpOnAxesImage);

  show_zoom = 1;
  set(handles.zoom_radiobutton,'Value',show_zoom);
%   setappdata(handles.stack_roi_rg,'show_zoom',show_zoom);

  % Set the pre/post timing
  setappdata(handles.stack_roi_rg,'bgIndex',-5:-1);  % background = 5 stacks pre
  setappdata(handles.stack_roi_rg,'respIndex',0:4);  % response   = 5 stacks post
  
  setappdata(handles.stack_roi_rg,'roi_defs',[]);
  
  % Update handles structure
  guidata(hObject, handles);

% if(nargin == 4 && ~isempty(varargin))
%   stackRoiParams=varargin{1};
%   if isfield(stackRoiParams,'filename')
%     fileToOpen=stackRoiParams.filename;
%     openFile(handles, fileToOpen);
%   end
%   if isfield(stackRoiParams,'posInPixel')
%     % Navigate to a particular position
%     frame_num = stackRoiParams.posInPixel(3);
%     set(handles.current_frame, 'string', num2str(frame_num));
%     current_frame_Callback(handles.current_frame, [], handles)
%     show_zoom_images(handles, stackRoiParams.posInPixel(1:2),frame_num);
%   end
% else
%   error('Usage: stack_roi(struct(''filename'', headerFile))');
% end


% UIWAIT makes stack_roi_rg wait for user response (see UIRESUME)
% uiwait(handles.stack_roi_rg);

function pixel=um2pixel(handles, um, dim)
h=getappdata(handles.stack_roi_rg, 'header');
if(dim=='x' || dim=='y')
  if isfield(h, 'pixel_spacing')
    if h.pixel_spacing(1) == h.pixel_spacing(2);
      h.um_per_pixel_xy = h.pixel_spacing(1);
    end
  elseif(h.um_per_pixel_xy==-1)
    h.um_per_pixel_xy=0.71;
  end
  pixel=um/h.um_per_pixel_xy;
else
  h.um_per_pixel_z=abs(diff(h.piezo_start_stop))/h.frames_per_stack;
  if(h.um_per_pixel_z==0)
    h.um_per_pixel_z=1; %
  end % if, piezo didn't move at all
  pixel=um/h.um_per_pixel_z;
end

function radInUm=getRadius(handles)
   radInUm=str2double(get(handles.ROIradius_edit, 'String'));

% plot or update circle/center/label
function hroi=fPlotRois(hAxes, roi)
%axes(hAxes); % NOTE: or another way, set parent property in line()

npts = 100;
th = linspace(0,2*pi,npts);
cth = cos(th);
sth = sin(th);
n_rois = length(roi);
hroi = zeros(n_rois,3);
for idx = 1:n_rois
  x = roi(idx).posInPixel(1) + roi(idx).xyradiusInPixel*cth;
  y = roi(idx).posInPixel(2) + roi(idx).xyradiusInPixel*sth;
  % the circle
  if(isfield(roi, 'hCircle'))
    % update data only
    set(roi(idx).hCircle, 'XData', x, 'YData', y);
    hroi(idx,1)=roi.hCircle;
  else
    hroi(idx,1) = line(x,y, 'parent',hAxes,'linewidth', 1, 'color', 'blue');
  end
  % the circle center:
  if(isfield(roi, 'hCenter'))
    % update data only
    set(roi(idx).hCenter, 'XData',roi(idx).posInPixel(1), 'YData',roi(idx).posInPixel(2));
    hroi(idx,2)=roi(idx).hCenter;
  else
    hroi(idx,2) = line(roi(idx).posInPixel(1), roi(idx).posInPixel(2),  'parent',hAxes,'linewidth', 2, 'marker', '.', 'color', 'green');
  end
end
% the label: (TEH: do these in parallel to speed performance)
pos = cat(1,roi.posInPixel);
r = cat(1,roi.xyradiusInPixel);
if(isfield(roi, 'hLabel'))
  % update data only
  set([roi.hLabel], 'Position', [pos(:,1)+r, pos(:,2)]);
  hroi(:,3)=[roi.hLabel];
else
  label = cat(1,roi.label);
  labelstr = cell(1,n_rois);
  for idx = 1:n_rois
    labelstr{idx} = sprintf('%d',label(idx));
  end
  hroi(:,3) = text(pos(:,1)+r, pos(:,2),labelstr);
  set(hroi(:,3), 'parent',hAxes,'color', 'r', 'FontSize', 14); % text color is red
end


function addRoi(handles, roi)
% PRE:
%    roi: a ROI that has 3 fields:
%                       posInPixel, xyradiusInPixel, label % field 1,2,3
% POST:
%    fig's appdata roi_defs is updated;
%    circle/center/label objects are created
hRois=fPlotRois(handles.image_axes, roi);
for idx = 1:length(roi)
  roi(idx).hCircle=hRois(idx,1);                                % field 4
  roi(idx).hCenter=hRois(idx,2);                                % field 5
  roi(idx).hLabel=hRois(idx,3);                                 % field 6
  roi(idx).showCircle=true;                                    % field 7 of 7
end

roi_defs=getappdata(handles.stack_roi_rg, 'roi_defs');
roi_defs=[roi_defs roi];
setappdata(handles.stack_roi_rg, 'roi_defs', roi_defs);

% drag_line(roi.hCenter, struct('type', 'both'));
for idx = 1:length(roi)
  set(roi(idx).hCenter, 'ButtonDownFcn', @onClickCenter);
  setappdata(roi(idx).hCenter, 'hCircle', roi(idx).hCircle);
  setappdata(roi(idx).hCircle, 'hCenter', roi(idx).hCenter);
  setappdata(roi(idx).hCenter, 'hLabel', roi(idx).hLabel);
  setappdata(roi(idx).hCircle, 'hLabel', roi(idx).hLabel);
  drag_circle(roi(idx).hCircle, struct('onDragStart', @onDragRoiStart, 'onDragDone', @onDragRoiDone));
end

set([roi.hCenter], 'visible', 'off');



function result=onMouseUpOnAxesImage(sender, ~)
  handles=guidata(sender);
  posInPixel=get(sender, 'currentpoint');
  curFrameIdx=getappdata(handles.stack_roi_rg, 'curFrameIdx');
   
   if(is_button_down(sender, 'middle')) % sender: the axes
          
     % add a new roi

   
     roi.posInPixel=[posInPixel(1,1:2) curFrameIdx];      % field 1
     radius=getRadius(handles);
     roi.xyradiusInPixel=um2pixel(handles, radius, 'x');  % field 2
     
     roi_defs=getappdata(handles.stack_roi_rg, 'roi_defs');
     if(isempty(roi_defs))
       existentLabels=[];
     else
       existentLabels=[roi_defs.label];
     end
     roi.label=find_first_slot(existentLabels);           % field 3
     
     addRoi(handles, roi);
     
     if get(handles.zoom_radiobutton,'Value') %getappdata(handles.stack_roi_rg,'show_zoom')
       show_zoom_images(handles, posInPixel, curFrameIdx);
     end
     
     result=0;
   elseif(is_button_down(sender,'right'))
     % make zoom feature work
     if get(handles.zoom_radiobutton,'Value') %getappdata(handles.stack_roi_rg,'show_zoom')
       show_zoom_images(handles, posInPixel, curFrameIdx);
     end
     
     
     result =0;
   else
     
     % show_popup_menu(sender);
     result=0;
   end

function updateShowCircleField(hCircle)
   fig=get_parent_fig(hCircle);
   roi_defs=getappdata(fig, 'roi_defs');
   idx=[roi_defs.hCircle]==hCircle;
   isVisible=strcmp(get(hCircle, 'visible'), 'on');
   roi_defs(idx).showCircle=isVisible;
   setappdata(fig, 'roi_defs', roi_defs);
   
function onClickCenter(sender, ~) % sender: the center
   if(is_button_down(sender, 'left'))
      hCircle=getappdata(sender, 'hCircle');
      hLabel=getappdata(sender, 'hLabel');
      if(strcmp(get(hCircle, 'visible'), 'on'))
	 set([hCircle, hLabel], 'visible', 'off');
      else
	 set([hCircle, hLabel], 'visible', 'on');
	 set(sender, 'visible', 'off');
      end
      updateShowCircleField(hCircle);
   end
   
   
function isContinue=onDragRoiStart(sender, ~)
   if(is_button_down(sender, 'middle')) % sender: the line obj
      isContinue=0;
   else
      isContinue=1;
   end

   
function result=isequal_roi(a,b)   
   result=isequal(round(a),round(b));

   
function onDragRoiDone(sender, event_args)
   fig=get_parent_fig(sender); % sender: the circle line
   isInDefTformMode=getappdata(fig, 'isInDefTformMode');
   
   newPosAndRad=event_args.def;
   roi_defs=getappdata(fig, 'roi_defs');
   idx=find([roi_defs.hCircle]==sender);
   if(isequal_roi([roi_defs(idx).posInPixel(1:2) roi_defs(idx).xyradiusInPixel], newPosAndRad))
      % if, just mouse down/up
      if(isInDefTformMode)
	 % select/deselect as control ROI
	 linecolor=get(sender, 'color');
	 if(isequal(linecolor, [0 0 1]))
	    linecolor='blue';
	 elseif(isequal(linecolor, [0 1 0]))
	    linecolor='green';
	 elseif(isequal(linecolor, [1 0 0]))
	    linecolor='red';
	 else
	    error('this line color should not be here');
	 end
	    
	 if(strcmp(linecolor, 'blue'))
	    set(sender, 'color', 'green'); % become ctrl roi
	 elseif(strcmp(linecolor, 'green'))
	    % TODO: may need to make sure only one ctrl roi is selected.
	    set(sender, 'color', 'red'); % become selected ctrl roi
	 else
	    set(sender, 'color', 'blue'); % become non-ctrl roi
	 end
      else
	 hLabel=getappdata(sender, 'hLabel');
	 set([sender,hLabel], 'visible', 'off');
	 set(getappdata(sender, 'hCenter'), 'visible', 'on');
	 updateShowCircleField(sender);
      end
   else
      % roi def is changed:
      roi_defs(idx).posInPixel(1:2)=newPosAndRad(1:2);
      roi_defs(idx).xyradiusInPixel=newPosAndRad(3);
      setappdata(fig, 'roi_defs', roi_defs);
      
      % update the position of center object and label object
      fPlotRois(get(sender, 'parent'), roi_defs(idx)); % NOTE: the line obj's parent is axes.
      
   end


% --- Outputs from this function are returned to the command line.
function varargout = stack_roi_rg_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ROIradius_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ROIradius_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ROIradius_edit as text
%        str2double(get(hObject,'String')) returns contents of ROIradius_edit as a double


% --- Executes during object creation, after setting all properties.
function ROIradius_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROIradius_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function current_frame_Callback(hObject, eventdata, handles)
% hObject    handle to current_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_frame as text
%        str2double(get(hObject,'String')) returns contents of current_frame as a double

fig=handles.stack_roi_rg;
curFrameIdx=getappdata(fig, 'curFrameIdx');
newFrameIdx=str2double(get(hObject,'String'));
stackSize=getappdata(fig, 'stackSize');
if(newFrameIdx<1 || newFrameIdx>stackSize(3))
  newFrameIdx = curFrameIdx;
end
setappdata(fig,'curFrameIdx', newFrameIdx);
update_images(handles);


% --- Executes during object creation, after setting all properties.
function current_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_pushbutton.
function load_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fig=handles.stack_roi_rg;

   
   if(isunix)
      file_filter='*.roidef';
   else
      file_filter='*.roidef';
   end
   [filename, pathname] = uigetfile(file_filter, 'Pick a ROI definition file');
   if(filename==0)
      return;
   end
   filename  =fullfile(pathname, filename);
   
   tt=load(filename, '-mat');
   % tt.roi_defs
   % TODO: delete old circles/centers objects, create new and install event handles like what add does

   % delete old rois and related objects
   roi_defs=getappdata(fig, 'roi_defs');
   if(~isempty(roi_defs))
      set([roi_defs.hCircle], 'visible', 'on');
      onDeleteRoi(hObject, []);
   end
   
   roi_defs=tt.roi_defs;
   

   % remove fields
   roi_defs=rmfield(roi_defs, {'posInUm', 'xyradiusInUm'});
   
%   for idx=1:length(roi_defs)
%      addRoi(handles, roi_defs(idx));
%   end
   addRoi(handles,roi_defs);


% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
headerFile=getappdata(handles.stack_roi_rg, 'headerFile');
   header=getappdata(handles.stack_roi_rg, 'header');
   [pathstr,name,ext] = fileparts(headerFile);
   [thefilename,pathname] = uiputfile([pathstr name '.roidef'],'Save data to file...');
   roi_defs=getappdata(handles.stack_roi_rg, 'roi_defs');
   
   roi_defs=rmfield(roi_defs, {'hCircle','hCenter','hLabel', 'showCircle'});
   
   % convert pixel to um
   pixelPerUm(1)=um2pixel(handles, 1, 'x');
   pixelPerUm(2)=pixelPerUm(1);
   pixelPerUm(3)=um2pixel(handles, 1, 'z');
   for idxRoi=1:length(roi_defs)
      roi_defs(idxRoi).posInUm=roi_defs(idxRoi).posInPixel./pixelPerUm; % NOTE: this is relative pos to piezo start pos
      roi_defs(idxRoi).xyradiusInUm=roi_defs(idxRoi).xyradiusInPixel/pixelPerUm(1);
   end
   
   % when save, save two new fields (posInUm, xyradiusInUm) in roi_defs, 
   %    and also save header and pixelPerUm
   save([pathname thefilename],'header','roi_defs', 'pixelPerUm', '-mat');


% --- Executes on button press in done_pushbutton.
function done_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to done_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;


% --- Executes on button press in contrast_pushbutton.
function contrast_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
newclim = imrangegui(get(handles.main_image, 'CData'), get(handles.image_axes, 'CLim'), 0);
if ~isempty(newclim)
  set([handles.image_axes handles.zoom_axes], 'CLim', newclim);
end


% --- Executes on button press in deselect_pushbutton.
function deselect_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deselect_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.stack_roi_rg;
   roi_defs=getappdata(fig, 'roi_defs');
   for idx=1:length(roi_defs)
      roi_defs(idx).showCircle=0;
   end
   setappdata(fig, 'roi_defs', roi_defs);
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   display_plot(fig);
   
   % set([roi_defs.hCircle roi_defs.hLabel], 'visible', 'off');
   % set([roi_defs.hCenter], 'visible', 'on');
   


function current_stack_Callback(hObject, eventdata, handles)
% hObject    handle to current_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_stack as text
%        str2double(get(hObject,'String')) returns contents of
%        current_stack as a double
fig=handles.stack_roi_rg;
curStackIdx=getappdata(fig, 'curStackIdx');
newStackIdx=str2double(get(hObject,'String'));
stack=getappdata(fig, 'stack');
stackSize=stack.size;
if(newStackIdx<1 || newStackIdx>stackSize(4))
  newStackIdx = curStackIdx;
end
setappdata(fig,'curStackIdx', newStackIdx);
update_images(handles);



% --- Executes during object creation, after setting all properties.
function current_stack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dfof_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dfof_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dfof_min_edit as text
%        str2double(get(hObject,'String')) returns contents of dfof_min_edit as a double
update_images(handles)


% --- Executes during object creation, after setting all properties.
function dfof_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfof_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dfof_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to dfof_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dfof_max_edit as text
%        str2double(get(hObject,'String')) returns contents of
%        dfof_max_edit as a double
update_images(handles)


% --- Executes during object creation, after setting all properties.
function dfof_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dfof_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function repeat_edit_Callback(hObject, eventdata, handles)
% hObject    handle to repeat_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repeat_edit as text
%        str2double(get(hObject,'String')) returns contents of repeat_edit as a double
update_images(handles);


% --- Executes during object creation, after setting all properties.
function repeat_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repeat_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function openFile(handles, headerFile)
   smm = stackmm(headerFile); % stack memory map
   stackSize = smm.size;
   setappdata(handles.stack_roi_rg,'stack',smm);
   setappdata(handles.stack_roi_rg,'stackSize',stackSize);
   curFrameIdx=1;
   curStackIdx=1;
   setappdata(handles.stack_roi_rg, 'headerFile', headerFile);
   h=imreadheader(headerFile);
   setappdata(handles.stack_roi_rg, 'header', h);
   setappdata(handles.stack_roi_rg, 'curFrameIdx', curFrameIdx);
   setappdata(handles.stack_roi_rg, 'curStackIdx', curStackIdx);
   
   stim_lookup = h.stim_lookup;
   toffset = [getappdata(handles.stack_roi_rg,'bgIndex'), getappdata(handles.stack_roi_rg,'respIndex')];
   [onset_list,ustim] = find_stimulus_start(stim_lookup,[],[max(1,1-min(toffset)) min(h.nstacks,h.nstacks-max(toffset))]);
   setappdata(handles.stack_roi_rg, 'onset_list', onset_list);
   setappdata(handles.stack_roi_rg, 'ustim', ustim);
   
   stim_labels = h.stim_labels;
   setappdata(handles.stack_roi_rg, 'stim_labels', stim_labels);

   popup_labels{1} = 'grey';
   popup_labels(2:length(stim_labels)+1) = stim_labels;
   set(handles.display_popupmenu,'String',popup_labels)
      
   set(handles.image_axes,'XLim',[0 stackSize(1)]+0.5,'YLim',[0 stackSize(2)]+0.5);
   
   fig = handles.stack_roi_rg;
   display_plot(fig);
   
   center = round(stackSize(1:2)/2);
%    center(2,:) = center;
   show_zoom_images(handles, center, curFrameIdx)
   
   
   %--------------------------
   
function update_images(handles)
  display_plot(handles);
  if get(handles.zoom_radiobutton,'Value')
    zc = getappdata(handles.stack_roi_rg,'zoom_center');
    curFrameIdx=getappdata(handles.stack_roi_rg, 'curFrameIdx');
    show_zoom_images(handles, zc, curFrameIdx)
  end

function display_plot(handles)
  if ~isstruct(handles)
    fig = handles;
    handles = guidata(fig);
  end
    
    display_mode = get(handles.display_popupmenu,'Value');
    
    if display_mode == 1 % grey mode
      
      set(handles.dfof_max_edit,'Visible','off');
      set(handles.dfof_min_edit,'Visible','off');
      set(handles.repeat_edit,'Visible','off');
      set(handles.repeat_text,'Visible','off');
      set(handles.dfof_text,'Visible','off');
      set(handles.checkboxWhiteBg,'Visible','off');
      set(handles.checkboxSubtract,'Visible','off');
      
      set(handles.current_stack,'Visible','on');
      set(handles.stack_text,'Visible','on');
      
      display_grey(handles);
    else % rg
      
      set(handles.dfof_max_edit,'Visible','on');
      set(handles.dfof_min_edit,'Visible','on');
      set(handles.repeat_edit,'Visible','on');
      set(handles.repeat_text,'Visible','on');
      set(handles.dfof_text,'Visible','on');
      set(handles.checkboxWhiteBg,'Visible','on');
      set(handles.checkboxSubtract,'Visible','on');
      
      set(handles.current_stack,'Visible','off');
      set(handles.stack_text,'Visible','off');
      
      display_rg(handles)  
    end

function imrgb = calculate_rg(handles,rng)
  smm = getappdata(handles.stack_roi_rg,'stack');
  onset_list = getappdata(handles.stack_roi_rg, 'onset_list');
  current_stim = get(handles.display_popupmenu,'Value') - 1;
  current_repeat = str2double(get(handles.repeat_edit,'String'));
  
  num_repeats = length(onset_list{1});
  if current_repeat < 0 || current_repeat > num_repeats % taking care of bad user input
    if current_repeat < 0
      current_repeat = 1;
    else
      current_repeat = num_repeats;
    end
    set(handles.repeat_edit,'String',num2str(current_repeat));
  end
  indx = onset_list{current_stim}(current_repeat);
  setappdata(handles.stack_roi_rg,'curStackIdx', indx);
  
  clim_raw = get(handles.image_axes,'CLim');
  clim_dfof(1) = str2double(get(handles.dfof_min_edit,'String'));
  clim_dfof(2) = str2double(get(handles.dfof_max_edit,'String'));
  if clim_dfof(1)>=clim_dfof(2) % taking care of bad user input
    clim_dfof(2) = clim_dfof(1)+0.1;
    set(handles.dfof_max_edit,'String',num2str(clim_dfof(2)));
  end
  
  im1 = mean(single(smm(rng{:}, indx+getappdata(handles.stack_roi_rg,'bgIndex'))),4); %bg
  im2 = mean(single(smm(rng{:}, indx+getappdata(handles.stack_roi_rg,'respIndex'))),4); %stim
  % this flip is to make it compatible with Zhongsheng's stack_roi code
  im1 = permute(im1,[2 1 3]);
  im2 = permute(im2,[2 1 3]);
   
  % Calculate dfof
  dfof = im2./im1-1;
  dfof(im1==0) = 0;
  if (get(handles.checkboxSubtract,'Value') == 1)
    % Subtract the control response
    ctrl = getappdata(handles.stack_roi_rg,'ctrl');
    ctrl = permute(ctrl(rng{:}),[2 1 3]);
    r = dfof./ctrl;
    alpha = nanmedian(r(:));
    dfof = dfof - alpha*ctrl;
  end
  
  if get(handles.checkboxWhiteBg,'Value')
    % Display on a white background
    raw = ones(size(im1),'single');
  else
    % Display on the grayscale background
    % Scale the image
    raw = im1 - clim_raw(1);
    raw(raw<0) = 0;
    raw = raw/diff(clim_raw);
    raw(raw>1) = 1;
  end
  
  % Scale the dfof
  if (clim_dfof(2) > 0)
    mask = dfof > 0;
    dfof(mask) = dfof(mask)/clim_dfof(2);
  end
  if (clim_dfof(1) < 0)
    mask = dfof < 0;
    dfof(mask) = dfof(mask)/abs(clim_dfof(1));
  end
  dfof(dfof>1) = 1;
  dfof(dfof<-1) = -1;
  
  % Make the colorized image
  imrgb = colorize_pixels(raw,dfof);  
  
function display_rg(handles)
  frame_num = getappdata(handles.stack_roi_rg, 'curFrameIdx');
  set(handles.current_frame, 'string', num2str(frame_num));
  
  rng = {':',':',frame_num};
  imrgb = calculate_rg(handles,rng);
  
  set(handles.main_image, 'CData', imrgb);
  
  show_frame_rois(handles);

    
function display_grey(handles)
  
   idxFrame = getappdata(handles.stack_roi_rg, 'curFrameIdx');
   idxStack = getappdata(handles.stack_roi_rg, 'curStackIdx');
  
   stack = getappdata(handles.stack_roi_rg,'stack');
   stackSize=stack.size;
   
   imageData = stack(:,:,idxFrame,idxStack); % when get data in matlab, index is 1-based
   imageData=squeeze(imageData);
   imageData=imageData';

   imageW=size(imageData, 2);
   imageH=size(imageData, 1); 
   
   axsize = [imageW imageH];

   
   oldPosAxes=get(handles.image_axes, 'position');
   %set(handles.image_axes, 'position', [oldPosAxes(1:2) axsize]);
   set(handles.image_axes, 'position', oldPosAxes);
   oldPosFig=get(handles.stack_roi_rg, 'position');
   
   set(handles.main_image,'CData',imageData);
   
   set(handles.current_frame, 'string', num2str(idxFrame));
   set(handles.current_stack, 'string', num2str(idxStack));
   
   show_frame_rois(handles);

   
   
  function show_frame_rois(handles)
    
    % show rois if needed
    roi_defs=getappdata(handles.stack_roi_rg, 'roi_defs');
    idxFrame = getappdata(handles.stack_roi_rg, 'curFrameIdx');
    if(~isempty(roi_defs))
        % tt=findobj([roi_defs.hCircle roi_defs.hCenter roi_defs.hLabel], 'visible', 'on');
        % set(tt, 'visible', 'on');
        tt=cat(1,roi_defs.posInPixel);
        idxRoisToShow=find(tt(:,3)==idxFrame);
        
        for idx=make_vector(idxRoisToShow, 'row')
            % apply the tform
            roi=roi_defs(idx);
            fPlotRois(handles.image_axes, roi); % move the position
        end
        % Set the visibility flags appropriately
        showCircle = logical([roi_defs(idxRoisToShow).showCircle]);
        circIdx = idxRoisToShow(showCircle);
        centerIdx = idxRoisToShow(~showCircle);
        if ~isempty(circIdx)
            set([roi_defs(circIdx).hCircle roi_defs(circIdx).hLabel], 'visible', 'on');
            set([roi_defs(circIdx).hCenter], 'visible', 'off');
        end
        if ~isempty(centerIdx)
            set([roi_defs(centerIdx).hCircle roi_defs(centerIdx).hLabel], 'visible', 'off');
            set([roi_defs(centerIdx).hCenter], 'visible', 'on');
        end
        idxRoisToHide=setdiff(1:length(roi_defs), idxRoisToShow);
        if(~isempty(idxRoisToHide))
            set([roi_defs(idxRoisToHide).hCircle roi_defs(idxRoisToHide).hLabel ...
                roi_defs(idxRoisToHide).hCenter], 'visible', 'off');
        end
    end
    


% --- Executes on selection change in display_popupmenu.
function display_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to display_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns display_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from display_popupmenu
update_images(handles)

% --- Executes during object creation, after setting all properties.
function display_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function onIncFrameIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx')+1;
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curFrameIdx==stackSize(3)+1)
      curFrameIdx=curFrameIdx-1;
   end
   setappdata(fig,'curFrameIdx', curFrameIdx);
   update_images(guidata(fig));

function onDecFrameIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx')-1;
   if(curFrameIdx==0)
      curFrameIdx=1;
   end
   setappdata(fig,'curFrameIdx', curFrameIdx);
   update_images(guidata(fig));

function onIncStackIdx(sender, event)
   fig=get_parent_fig(sender);
   curStackIdx=getappdata(fig, 'curStackIdx')+1;
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curStackIdx==stackSize(4)+1)
      %curStackIdx=curStackIdx-1;      % stay at the last stack
      curStackIdx=1;                  % wrap around to the beginning
   end
   setappdata(fig,'curStackIdx', curStackIdx);
   update_images(guidata(fig));

   
function onDecStackIdx(sender, event)
   fig=get_parent_fig(sender);
   curStackIdx=getappdata(fig, 'curStackIdx')-1;
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curStackIdx==0)
      %curStackIdx=0;                % stack with the first stack
      curStackIdx = stackSize(4);   % wrap around to the end
   end
   setappdata(fig,'curStackIdx', curStackIdx);
   update_images(guidata(fig));

  
   

function onDeleteRoi(sender, event)
   fig=get_parent_fig(sender);
   handles=guidata(fig);
   roi_defs=getappdata(fig, 'roi_defs');
   visibleRois=findobj([roi_defs.hCircle], 'visible', 'on');
   indicesToDel=[];
   for hroi=visibleRois' % for's loop set must be row vector
      idx=find([roi_defs.hCircle]==hroi);
      indicesToDel(end+1)=idx;
      delete(roi_defs(idx).hCircle);
      delete(roi_defs(idx).hCenter);
      delete(roi_defs(idx).hLabel);
   end
   roi_defs(indicesToDel)=[];
   
   setappdata(fig, 'roi_defs', roi_defs);

   
   % --- Executes on button press in zoom_radiobutton.
  function zoom_radiobutton_Callback(hObject, eventdata, handles)
    % hObject    handle to zoom_radiobutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Hint: get(hObject,'Value') returns toggle state of zoom_radiobutton
    show_zoom = get(hObject,'Value');
    if ~show_zoom
      set(handles.zoom_image,'Visible','off');
      set(handles.zoom_text,'Visible','off');
      set(handles.zoom_size_edit,'Visible','off');
    else
      set(handles.zoom_image,'Visible','on');
      set(handles.zoom_text,'Visible','on');
      set(handles.zoom_size_edit,'Visible','on');
    end



function zoom_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_size_edit as text
%        str2double(get(hObject,'String')) returns contents of zoom_size_edit as a double
  zoom_size = str2double(get(hObject,'String'));
  % Check for valid input
  zoom_size_default = 100;
  zoom_size_min = 3;
  stackSize = getappdata(handles.stack_roi_rg,'stackSize');
  zoom_size_max = floor(min(stackSize(1:2))/2-1);
  if isnan(zoom_size)
    zoom_size = zoom_size_default;
  end
  zoom_size = max(zoom_size_min,min(zoom_size,zoom_size_max));
  set(handles.zoom_size_edit,'String',num2str(zoom_size))
  
  zoomimage = get(handles.zoom_image(1),'CData');
  szz = size(zoomimage);
  szz = szz(1);
  if (szz >= 2*zoom_size+1)
    % We can just adjust the display limits
    extra = (szz-(2*zoom_size+1))/2;
    lim = [extra+1 szz-extra];
    set(handles.zoom_axes,'XLim',lim,'YLim',lim);
  else
    % Need to recalculate
    zc = getappdata(handles.stack_roi_rg,'zoom_center');
    curFrameIdx=getappdata(handles.stack_roi_rg, 'curFrameIdx');
    show_zoom_images(handles, zc, curFrameIdx)
  end

  

% --- Executes during object creation, after setting all properties.
function zoom_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

  function show_zoom_images(handles, posInPixel, curFrameIdx)
    delete(findobj(handles.stack_roi_rg,'Tag','overlay'))
    stackSize = getappdata(handles.stack_roi_rg,'stackSize');
    
    posInPixel = round(posInPixel(1,:));
    setappdata(handles.stack_roi_rg,'zoom_center',posInPixel(1:2));
    
    zoom_size = str2double(get(handles.zoom_size_edit, 'String'));
    
    axsize = [zoom_size*2+1 zoom_size*2+1];
    set(handles.zoom_axes,'XLim',[0 axsize(1)]+0.5,'YLim',[0 axsize(2)]+0.5);
    
    display_mode = get(handles.display_popupmenu,'Value');
    
    xsnip = posInPixel(1)-zoom_size:posInPixel(1)+zoom_size;
    keepx = xsnip >=1 & xsnip <= stackSize(1);
    ysnip = posInPixel(2)-zoom_size:posInPixel(2)+zoom_size;
    keepy = ysnip >=1 & ysnip <= stackSize(2);
    show_frames = curFrameIdx + getappdata(handles.stack_roi_rg,'zoom_zoffsets');
    keepframes = show_frames >=1 & show_frames <= stackSize(3);

    if display_mode == 1 % show grey scale
      smm = getappdata(handles.stack_roi_rg,'stack');
      idxStack = getappdata(handles.stack_roi_rg, 'curStackIdx'); 
      im = smm(xsnip(keepx),ysnip(keepy),show_frames(keepframes),idxStack);
      imageData = nan([axsize length(show_frames)],'single');
      imageData(keepx,keepy,keepframes) = im;
      imageData = permute(imageData,[2 1 3]);
      colormap(gray(256));
      
    else % display mode in stim colorize dfof
      rng = {xsnip(keepx),ysnip(keepy),show_frames(keepframes)};
      imrgb = calculate_rg(handles,rng);
      if ndims(imrgb) == 3
        imrgb = reshape(imrgb,[sum(keepx) sum(keepy) 1 3]);
      end
      imageData = nan([axsize length(show_frames) 3],'single');
      imageData(keepx,keepy,keepframes,:) = imrgb;
    end
    
    % Display the data
    for i = 1:length(handles.zoom_image)
      set(handles.zoom_image(i),'CData',squeeze(imageData(:,:,i,:)));
    end
    
    % show the ROI defs
    roi_defs = getappdata(handles.stack_roi_rg, 'roi_defs');
    if isempty(roi_defs)
      return
    end
    pos = cat(1,roi_defs.posInPixel);
    xsnip = posInPixel(1) + [-1 1]*zoom_size;
    ysnip = posInPixel(2) + [-1 1]*zoom_size;
    isInRectFlag = is_in_rect(pos(:,1:2),[xsnip(1) ysnip(1) xsnip(2) ysnip(2)]);
    hax = [handles.zoom1_axes handles.zoom2_axes handles.zoom3_axes handles.zoom4_axes handles.zoom5_axes];
    for idx = 1:5
      keepFlag = pos(:,3) == show_frames(idx);
      haxg = axes('Position',get(hax(idx),'Position'),'Visible','off','Tag','overlay');
      hcirc = copyobj([roi_defs(keepFlag).hCircle],haxg);
      set(hcirc,'Visible','on')
      hlbl = copyobj([roi_defs(keepFlag & isInRectFlag).hLabel],haxg);
      set(hlbl,'Visible','on')
      set(haxg,'XLim',xsnip([1 end]),'YLim',ysnip([1 end]))
    end
    
    
    
    


% --- Executes on button press in select_pushbutton.
function select_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to select_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

   fig=handles.stack_roi_rg;
   roi_defs=getappdata(fig, 'roi_defs');
   for idx=1:length(roi_defs)
      roi_defs(idx).showCircle=1;
   end
   setappdata(fig, 'roi_defs', roi_defs);
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   display_plot(fig);


% --- Executes on button press in checkboxSubtract.
function checkboxSubtract_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxSubtract (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxSubtract
  ctrl = getappdata(handles.stack_roi_rg,'ctrl');
  if isempty(ctrl)
    % Prompt the user to specify the identity of the control stimulus
    stim_labels = getappdata(handles.stack_roi_rg, 'stim_labels');
    choice = listdlg('ListString',stim_labels,'SelectionMode','single','PromptString','Choose the control stimulus');
    if isempty(choice)
      set(hObject,'Value',0);
      return
    else
      % Calculate the response to the control stimulus
      onset_list = getappdata(handles.stack_roi_rg, 'onset_list');
      onset = onset_list{choice};
      smm = getappdata(handles.stack_roi_rg,'stack');
      bgIndex = getappdata(handles.stack_roi_rg,'bgIndex');
      respIndex = getappdata(handles.stack_roi_rg,'respIndex');
      progress_bar(struct('max',length(onset),'progress',0,'what','Calculating response to control...'));
      for i = 1:length(onset)
        im1 = mean(single(smm(:,:,:,onset(i)+bgIndex)),4);
        im2 = mean(single(smm(:,:,:,onset(i)+respIndex)),4);
        dfof = im2./im1 - 1;
        dfof(im1==0) = 0;
        if isempty(ctrl)
          ctrl = dfof;
        else
          ctrl = ctrl + dfof;
        end
        progress_bar(struct('max',length(onset),'progress',i,'what','Calculating response to control...'));
      end
      ctrl = ctrl/length(onset);
      setappdata(handles.stack_roi_rg,'ctrl',ctrl);
    end
  end
  update_images(handles);


% --- Executes on button press in checkboxWhiteBg.
function checkboxWhiteBg_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxWhiteBg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxWhiteBg
  update_images(handles);
