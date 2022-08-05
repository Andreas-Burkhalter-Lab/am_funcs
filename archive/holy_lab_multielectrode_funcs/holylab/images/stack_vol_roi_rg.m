function varargout = stack_vol_roi_rg(varargin)
%
% % keybinding:
%    up/down arrow: change frame #
%    left/right arrow: change roi radius: NOT FUNCTIONAL
%    DEL/d: delete selected ROIs
%    pageup/down: change stack #
%    middle-click: begin/end definition of a new volmetric ROI
%    ---- while defining a volumetric ROI----
%           left-click: add a vertex to the ROI for this frame
%           right-click: remove a vertex from the ROI for this frame
%    ---- while not defining a volumetric ROI----
%           right-click: look at clicked position in frames 
%                        before and after (in zoom panel)
%           left-click on center: select a roi
%    control+e: enter/exit "edit mode" to edit existing ROI
%         *to edit an existing ROI, after pressing control+e,
%          click on the ROI center point to select the ROI,
%          then use right/left clicks just as described for
%          "while defining a volumetric ROI" above
%         *to finish editing, first hit "control+e" again to exit
%          editing mode, then middle-click to close the ROI.
%       
%    note/todo: stack_vol_roi_rg was written in a way that is not
%               back-compatible with circle-rois.  In the future, the 
%               added features in this version may be fused with 
%               stack_roi_rg (with non-volumetric circular ROIs) 
%               to allow both types to be drawn.
%
% usage: stack_roi(struct('filename','aob_gcamp2.imagine'))
%
% Written by Julian P. Meeks,
% Adapted from stack_roi (written by Zhongsheng Guo)
% and stack_roi_rg (written by Diwakar Turaga)
% 2010_12_09
%
% STACK_VOL_ROI_RG M-file for STACK_VOL_ROI_RG.fig
%      STACK_VOL_ROI_RG, by itself, creates a new STACK_VOL_ROI_RG or raises the existing
%      singleton*.
%
%      H = STACK_VOL_ROI_RG returns the handle to a new STACK_VOL_ROI_RG or the handle to
%      the existing singleton*.
%
%      STACK_VOL_ROI_RG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_VOL_ROI_RG.M with the given input arguments.
%
%      STACK_VOL_ROI_RG('Property','Value',...) creates a new STACK_VOL_ROI_RG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STACK_VOL_ROI_RG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STACK_VOL_ROI_RG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STACK_VOL_ROI_RG

% Last Modified by GUIDE v2.5 09-Dec-2010 20:24:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @STACK_VOL_ROI_RG_OpeningFcn, ...
  'gui_OutputFcn',  @STACK_VOL_ROI_RG_OutputFcn, ...
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


% --- Executes just before STACK_VOL_ROI_RG is made visible.
function STACK_VOL_ROI_RG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STACK_VOL_ROI_RG (see VARARGIN)

% Choose default command line output for STACK_VOL_ROI_RG
handles.output = hObject;

% Set up blank images
set(handles.image_axes, ...
  'CLim',[500 10000]);

fig=hObject;
bind_shortcut(fig, 'uparrow', @onIncFrameIdx);
bind_shortcut(fig, 'downarrow', @onDecFrameIdx);
%bind_shortcut(fig, 'rightarrow', @onIncRoiRadius);
%bind_shortcut(fig, 'leftarrow', @onDecRoiRadius);
bind_shortcut(fig, 'd', @onDeleteRoi);
bind_shortcut(fig, 'delete', @onDeleteRoi);
bind_shortcut(fig, 'pageup', @onIncStackIdx);
bind_shortcut(fig, 'pagedown', @onDecStackIdx);
bind_shortcut(fig, 'control+e', @onEditToggle);

install_mouse_event_handler(handles.image_axes, 'up', @onMouseUpOnAxesImage);

show_zoom = 1;
set(handles.zoom_radiobutton,'Value',show_zoom);
setappdata(handles.stack_vol_roi_rg,'show_zoom',show_zoom);
setappdata(handles.stack_vol_roi_rg,'roiInProgress',0);
setappdata(handles.stack_vol_roi_rg,'rg_mode','raw');
setappdata(handles.stack_vol_roi_rg,'edit_status',0);

% Update handles structure
guidata(hObject, handles);

if(nargin == 4 && ~isempty(varargin))
  stackRoiParams=varargin{1};
  if ischar(stackRoiParams)
      temp = struct('filename',stackRoiParams);
      stackRoiParams = temp;
      clear temp;
  end
  fileToOpen=stackRoiParams.filename;
  openFile(handles, fileToOpen);
else
  error('Usage: stack_roi(struct(''filename'', headerFile))');
end


% UIWAIT makes STACK_VOL_ROI_RG wait for user response (see UIRESUME)
% uiwait(handles.STACK_VOL_ROI_RG);

function pixel=um2pixel(handles, um, dim)
h=getappdata(handles.stack_vol_roi_rg, 'header');
if(dim=='x' || dim=='y')
  if(h.um_per_pixel_xy==-1)
    h.um_per_pixel_xy=0.71;
  end
  pixel=um/h.um_per_pixel_xy;
else
  h.um_per_pixel_z=abs(diff(h.piezo_start_stop))/(h.frames_per_stack-1);
  if(h.um_per_pixel_z==0)
    h.um_per_pixel_z=1; %
  end % if, piezo didn't move at all
  pixel=um/h.um_per_pixel_z;
end

function radInUm=getRadius(handles)
   radInUm=str2num(get(handles.ROIradius_edit, 'String'));

% plot or update circle/center/label
function hroi=fPlotRois(hAxes, roi)
%axes(hAxes); % NOTE: or another way, set parent property in line()
curAx = getappdata(get(hAxes,'parent'),'curFrameIdx');
centerColor = [0.7 0 0.7];
tmp = [];
% for i = 1:length(roi.vtxInPixels{curAx})
    tmp = roi.vtxInPixels{curAx};
% end
if isempty(tmp)
    if isfield(roi,'hOutline')
        hroi{1} = roi.hOutline;
        hroi{1}(curAx) = 0;
    end
    % the circle center:
    if(isfield(roi, 'hCenter'))
        % update data only
        set(roi.hCenter, 'XData',roi.centerInPixels(1), 'YData',roi.centerInPixels(2) );
        hroi{2}=roi.hCenter;
    else
        hroi{2} = line(roi.centerInPixels(1), roi.centerInPixels(2),  'parent',hAxes,'linewidth', 2, 'marker', '.', 'color', centerColor );
    end
    % the label:
    if(isfield(roi, 'hLabel'))
        % update data only
        set(roi.hLabel, 'Position', [roi.centerInPixels(1)+5, roi.centerInPixels(2)] );
        hroi{3}=roi.hLabel;
    else
        hroi{3} = text(roi.centerInPixels(1)+5, roi.centerInPixels(2), num2str(roi.label));
        set(hroi{3}, 'color', centerColor, 'FontSize', 14 ); % text color is red
    end
else
    if ~isempty(roi.vtxInPixels{curAx})
        x = roi.vtxInPixels{curAx}(:,1);
        y = roi.vtxInPixels{curAx}(:,2);
        % the ROI Outline for this frame
        if(isfield(roi, 'hOutline'))
            if ~isempty(roi.hOutline)
                if ~isempty(intersect(curAx, find(roi.hOutline)))
                    % update data only
                    set(roi.hOutline(curAx), 'XData', x, 'YData', y );
                    hroi{1}=roi.hOutline;
                else
                    hroi{1}(curAx) = patch(x,y,'b','parent',hAxes,'linewidth', 1, 'edgecolor', 'blue','facecolor','none' );
                end
            else
                hroi{1}(curAx) = patch(x,y,'b','parent',hAxes,'linewidth', 1, 'edgecolor', 'blue','facecolor','none' );
            end
        else
            hroi{1}(curAx) = patch(x,y,'b','parent',hAxes,'linewidth', 1, 'edgecolor', 'blue','facecolor','none' );
        end
        % the circle center:
        if isfield(roi, 'hCenter')
            % update data only
%             set(roi.hCenter, 'XData',roi.centerInPixels(1), 'YData',roi.centerInPixels(2));
            hroi{2}=roi.hCenter;
        else
            hroi{2} = line(roi.centerInPixels(1), roi.centerInPixels(2),  'parent',hAxes,'linewidth', 2, 'marker', '.', 'color', centerColor );
        end
        % the label:
        if isfield(roi, 'hLabel')
            % update data only
%             set(roi.hLabel, 'Position', [roi.centerInPixels(1)+5, roi.centerInPixels(2)]);
            hroi{3}=roi.hLabel;
        else
            hroi{3} = text(roi.centerInPixels(1)+5, roi.centerInPixels(2), num2str(roi.label));
            set(hroi{3}, 'color', centerColor, 'FontSize', 14 ); % text color is red
        end
    end
end


function updateRoi(handles, roi)
% PRE:
%    roi: a ROI that has 3 fields:
%                       posInPixel, xyradiusInPixel, label % field 1,2,3
% POST:
%    fig's appdata roi_defs is updated;
%    circle/center/label objects are created
curFrame = getappdata(handles.stack_vol_roi_rg,'curFrameIdx');
hRois=fPlotRois(handles.image_axes, roi);
if ~isempty(hRois{1}) && length(hRois{1}) >= curFrame
    roi.hOutline(curFrame)=hRois{1}(curFrame);     % field 4
else
    roi.hOutline = [];
end
roi.hCenter=hRois{2};                                % field 5
roi.hLabel=hRois{3};                                 % field 6
roi.showOutline=1;                                    % field 7 of 7

roi_defs=getappdata(handles.stack_vol_roi_rg, 'roi_defs');
if isempty(roi_defs)
    roi_defs = fillRoi(struct);
end
idx = find([roi_defs.label]==roi.label,1);
if ~isempty(idx)
    roi_defs(idx)=fillRoi(roi);
elseif isempty(roi_defs(end).label)
    roi_defs(1) = fillRoi(roi);
else
    roi_defs(end+1)=fillRoi(roi);
end
setappdata(handles.stack_vol_roi_rg, 'roi_defs', roi_defs);

if isempty(get(roi.hCenter,'ButtonDownFcn'))
    set(roi.hCenter, 'ButtonDownFcn', @onClickCenter);
end
setappdata(roi.hCenter,'hCenter',roi.hCenter);
setappdata(roi.hCenter,'hLabel',roi.hLabel);
setappdata(roi.hCenter,'hOutline',roi.hOutline);
if ~isempty(roi.hOutline)
    set(roi.hOutline(curFrame),'hittest','off');
end
set(roi.hCenter, 'visible', 'on');


function roi_out = fillRoi(roi_in)
% to ensure standardization of roi_defs structure array
roi_out = struct('label',[],'centerInPixels',[],'vtxInPixels',[],'hOutline',[],'hCenter',[],'hLabel',[],'showOutline',[]);
if isfield(roi_in,'label'); roi_out.label = roi_in.label; end;
if isfield(roi_in,'centerInPixels'); roi_out.centerInPixels = roi_in.centerInPixels; end;
if isfield(roi_in,'vtxInPixels'); roi_out.vtxInPixels = roi_in.vtxInPixels; end;
if isfield(roi_in,'hOutline'); roi_out.hOutline = roi_in.hOutline; end;
if isfield(roi_in,'hCenter'); roi_out.hCenter = roi_in.hCenter; end;
if isfield(roi_in,'hLabel'); roi_out.hLabel = roi_in.hLabel; end;
if isfield(roi_in,'showOutline'); roi_out.showOutline = roi_in.showOutline; end;

function roi_out = findRoi(handles,roi)
%
roi_defs = getappdata(handles.stack_vol_roi_rg,'roi_defs');
existentLabels = [roi_defs.label];
idx = find(existentLabels == roi,1);
roi_out = roi_defs(idx);


function result=onMouseUpOnAxesImage(sender, event_args)
  handles=guidata(sender);
  posInPixel=get(sender, 'currentpoint');
  edit_status = getappdata(handles.stack_vol_roi_rg, 'edit_status');
  if isempty(edit_status)
      edit_status = 0;
      setappdata(handles.stack_vol_roi_rg,'edit_status',edit_status);
  end
  curFrameIdx=getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
  curRepeatIdx=getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
  roiInProgress = getappdata(handles.stack_vol_roi_rg,'roiInProgress');
  roi_defs=getappdata(handles.stack_vol_roi_rg, 'roi_defs');
  header = getappdata(handles.stack_vol_roi_rg, 'header');
  if(isempty(roi_defs))
      existentLabels=[];
  else
      existentLabels=[roi_defs.label];
  end
  
  if(is_button_down(sender, 'middle')) % sender: the axes
      if roiInProgress && edit_status == 0 % if an ROI is being defined, CLOSE editing of current ROI
          roi = findRoi(handles,roiInProgress);
          if length(roi) == 0
              setappdata(handles.stack_vol_roi_rg,'roiInProgress',0);
              result = 0;
              return;
          end
          tmp = [];
          for i = 1:length(roi.vtxInPixels)
              if ~isempty(roi.vtxInPixels{i})
                  roi.vtxInPixels{i} = [roi.vtxInPixels{i}; roi.vtxInPixels{i}(1,:)];
              end
          end
          if ~isempty(tmp)
              roi.vtxInPixels{curFrameIdx} = [roi.vtxInPixels{curFrameIdx};roi.vtxInPixels{curFrameIdx}(1,:)];
          end
          roiInProgress = 0;
          setappdata(handles.stack_vol_roi_rg,'roiInProgress',roiInProgress);
          status_message = ['finished defining ROI #' num2str(roi.label)];
          setappdata(handles.stack_vol_roi_rg,'status_message',status_message);
          update_status_message(handles.stack_vol_roi_rg);
      elseif edit_status == 1% if an ROI needs to be edited, don't start a new ROI, but allow other fxns to set editing ROI
          roi = [];
      else % if no ROI is being defined, start a new ROI
          roi.label = find_first_slot(existentLabels);
          roi.centerInPixels = [posInPixel(1,1:2) curFrameIdx];
          roi.vtxInPixels = cell([1 header.frames_per_stack]);
          roiInProgress = roi.label;
          setappdata(handles.stack_vol_roi_rg,'roiInProgress',roiInProgress);
          status_message = ['defining new ROI #' num2str(roiInProgress)];
          setappdata(handles.stack_vol_roi_rg,'status_message',status_message);
          update_status_message(handles.stack_vol_roi_rg);
      end
      if getappdata(handles.stack_vol_roi_rg,'show_zoom')
          show_zoom_images(handles, posInPixel, curFrameIdx);
      end
      
      if ~isempty(roi)
          updateRoi(handles,roi);
      end
      result = 0;
      
  elseif(is_button_down(sender,'left'))
        handles=guidata(sender);
        roiInProgress = getappdata(handles.stack_vol_roi_rg,'roiInProgress');
        if roiInProgress
            roi = findRoi(handles,roiInProgress);
            newvtx=get(sender, 'currentpoint');
            curFrameIdx=getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
            roi.vtxInPixels{curFrameIdx} = [roi.vtxInPixels{curFrameIdx}; newvtx(1,1:2)];
            updateRoi(handles,roi);
            if getappdata(handles.stack_vol_roi_rg,'show_zoom')
                show_zoom_images(handles, roi.centerInPixels, curFrameIdx);
            end
            setappdata(handles.stack_vol_roi_rg,'lastROI',roi.label);
            setappdata(handles.stack_vol_roi_rg,'hLastROI',roi.hCenter);
        end        
        result = 0; 
        
  elseif(is_button_down(sender,'right'))
      curpt = get(sender,'currentpoint');
      handles=guidata(sender);
      roiInProgress = getappdata(handles.stack_vol_roi_rg,'roiInProgress');
      curFrameIdx=getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
      if roiInProgress
          roi = findRoi(handles,roiInProgress);
          roi.vtxInPixels{curFrameIdx} = roi.vtxInPixels{curFrameIdx}(1:end-1,:);
          updateRoi(handles,roi);
      else
          if isappdata(handles.stack_vol_roi_rg,'lastROI')
              if ishandle(getappdata(handles.stack_vol_roi_rg,'hLastROI'))
                  setappdata(handles.stack_vol_roi_rg,'lastROI',[]);
                  setappdata(handles.stack_vol_roi_rg,'hLastROI',[]);
              end
          end
          result = 0;
      end
      if getappdata(handles.stack_vol_roi_rg,'show_zoom')
          show_zoom_images(handles, curpt(1,1:2), curFrameIdx);
      end
      result = 0;
  else
      % show_popup_menu(sender);
      result = 0;
  end

% function updateShowCircleField(hCircle)
%    fig=get_parent_fig(hCenter);
%    roi_defs=getappdata(fig, 'roi_defs');
%    idx=find([roi_defs.hCircle]==hCircle);
%    isVisible=strcmp(get(hCircle, 'visible'), 'on');
%    roi_defs(idx).showCircle=isVisible;
%    setappdata(fig, 'roi_defs', roi_defs);
   
function onEditToggle(sender, eventdata)
    fig = get_parent_fig(sender);
    edit_status = getappdata(fig,'edit_status');
    if isempty(edit_status)
        edit_status = 1;
        status_message = 'ROI edit mode selected';
    elseif edit_status == 0
        edit_status = 1;
        status_message = 'ROI edit mode selected';
    else
        edit_status = 0;
        status_message = 'ROI edit mode exited';
    end
    setappdata(fig,'edit_status',edit_status);
    setappdata(fig,'status_message',status_message);
    
    update_status_message(fig);
    
function update_status_message(fig)
    handles = getappdata(fig,'UsedByGUIData_m');
    set(handles.status_txt_message,'string',getappdata(fig,'status_message'));
   
   
function onClickCenter(sender, eventdata) 
    % sender: the center
    fig = get_parent_fig(sender);
    roi_defs = getappdata(fig,'roi_defs');
    idxROI = roi_defs(find([roi_defs.hCenter]==sender,1)).label;
    edit_status = getappdata(fig,'edit_status');
    handles = getappdata(fig,'UsedByGUIData_m');
    if is_button_down(sender,'left')
        if ~isempty(roi_defs)
            setappdata(fig,'hLastROI',sender);
            setappdata(fig,'lastROI',roi_defs(idxROI).label);
            status_message = ['selected ROI #' num2str(roi_defs(idxROI).label)];
            setappdata(fig ,'status_message',status_message);
            update_status_message(fig);
            set(handles.selectRoi,'string',num2str(idxROI));
        end
    elseif is_button_down(sender,'middle')
        if ~isempty(idxROI) && edit_status == 1
            setappdata(fig, 'roiInProgress', idxROI);
        end
    end
   
function isContinue=onDragRoiStart(sender, event_args)
   if(is_button_down(sender, 'middle')) % sender: the line obj
      isContinue=0;
   else
      isContinue=1;
   end

   
function result=isequal_roi(a,b)   
   result=isequal(round(a),round(b));

   
% function onDragRoiDone(sender, event_args)
%    fig=get_parent_fig(sender); % sender: the circle line
%    isInDefTformMode=getappdata(fig, 'isInDefTformMode');
%    
%    newPosAndRad=event_args.def;
%    roi_defs=getappdata(fig, 'roi_defs');
%    idx=find([roi_defs.hCircle]==sender);
%    if(isequal_roi([roi_defs(idx).posInPixel(1:2) roi_defs(idx).xyradiusInPixel], newPosAndRad))
%       % if, just mouse down/up
%       if(isInDefTformMode)
% 	 % select/deselect as control ROI
% 	 linecolor=get(sender, 'color');
% 	 if(isequal(linecolor, [0 0 1]))
% 	    linecolor='blue';
% 	 elseif(isequal(linecolor, [0 1 0]))
% 	    linecolor='green';
% 	 elseif(isequal(linecolor, [1 0 0]))
% 	    linecolor='red';
% 	 else
% 	    error('this line color should not be here');
% 	 end
% 	    
% 	 if(strcmp(linecolor, 'blue'))
% 	    set(sender, 'color', 'green'); % become ctrl roi
% 	 elseif(strcmp(linecolor, 'green'))
% 	    % TODO: may need to make sure only one ctrl roi is selected.
% 	    set(sender, 'color', 'red'); % become selected ctrl roi
% 	 else
% 	    set(sender, 'color', 'blue'); % become non-ctrl roi
% 	 end
%       else
% 	 hLabel=getappdata(sender, 'hLabel');
% 	 set([sender,hLabel], 'visible', 'off');
% 	 set(getappdata(sender, 'hCenter'), 'visible', 'on');
% 	 updateShowCircleField(sender);
%       end
%    else
%       % roi def is changed:
%       roi_defs(idx).posInPixel(1:2)=newPosAndRad(1:2);
%       roi_defs(idx).xyradiusInPixel=newPosAndRad(3);
%       setappdata(fig, 'roi_defs', roi_defs);
%       
%       % update the position of center object and label object
%       fPlotRois(get(sender, 'parent'), roi_defs(idx)); % NOTE: the line obj's parent is axes.
%       
%    end


% --- Outputs from this function are returned to the command line.
function varargout = STACK_VOL_ROI_RG_OutputFcn(hObject, eventdata, handles)
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

fig=handles.stack_vol_roi_rg;
curFrameIdx=getappdata(fig, 'curFrameIdx');
newFrameIdx=str2double(get(handles.current_frame,'String'));
stack=getappdata(fig, 'stack');
stackSize=stack.size;
if(newFrameIdx<1 || newFrameIdx>stackSize(3))
  newFrameIdx = curFrameIdx;
end
setappdata(fig,'curFrameIdx', newFrameIdx);
display_plot(fig);


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
fig=handles.stack_vol_roi_rg;

   
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
      set([roi_defs.hCenter], 'visible', 'on');
      setappdata(fig,'lastROI',[roi_defs.label]);
      setappdata(fig,'hLastROI',[roi_defs.hCenter]);
      onDeleteRoi(hObject, []);
   end
   
   roi_defs_in=tt.roi_defs; clear tt;
   

   % remove fields
   roi_defs_in=rmfield(roi_defs_in, {'roiVolumeInUm3', 'vtxInUm','centerInUm'});
   orig_frame = getappdata(fig,'curFrameIdx');
   for i=1:length(roi_defs_in)
       for j = 1:length(roi_defs_in(i).vtxInPixels)
           setappdata(fig,'curFrameIdx',j);
           if j == 1
               updateRoi(handles, roi_defs_in(i));
           else
               roi_defs = getappdata(fig,'roi_defs');
               updateRoi(handles, roi_defs(i));
               %show_frame_rois(handles);
           end
       end
   end
   
   for i = 1:length(roi_defs_in(end).vtxInPixels)
       setappdata(fig,'curFrameIdx',i);
       show_frame_rois(handles);
   end
   
   setappdata(fig,'curFrameIdx',orig_frame);
   
   current_frame_Callback(hObject, eventdata, handles)


% --- Executes on button press in save_pushbutton.
function save_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to save_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
headerFile=getappdata(handles.stack_vol_roi_rg, 'headerFile');
   header=getappdata(handles.stack_vol_roi_rg, 'header');
   [pathstr,name,ext] = fileparts(headerFile);
   [thefilename,pathname] = uiputfile([pathstr name '.roidef'],'Save data to file...');
   roi_defs=getappdata(handles.stack_vol_roi_rg, 'roi_defs');
   
   roi_defs=rmfield(roi_defs, {'hOutline','hCenter','hLabel', 'showOutline'});
   
   % convert pixel to um
   pixelPerUm(1)=um2pixel(handles, 1, 'x');
   pixelPerUm(2)=pixelPerUm(1);
   pixelPerUm(3)=um2pixel(handles, 1, 'z');
   volPerPixel = pixelVol(pixelPerUm); % volume in um^3

   for idxRoi=1:length(roi_defs)
       temp = [];
       roi_defs(idxRoi).centerInUm = roi_defs(idxRoi).centerInPixels./pixelPerUm; % NOTE: this is relative pos to piezo start pos
       for idxVtx=1:length(roi_defs(idxRoi).vtxInPixels)
           if ~isempty(roi_defs(idxRoi).vtxInPixels{idxVtx})
               roi_defs(idxRoi).vtxInUm{idxVtx}(:,1) = roi_defs(idxRoi).vtxInPixels{idxVtx}(:,1)/pixelPerUm(1);
               roi_defs(idxRoi).vtxInUm{idxVtx}(:,2) = roi_defs(idxRoi).vtxInPixels{idxVtx}(:,2)/pixelPerUm(2);           
               temp = [temp polyarea(roi_defs(idxRoi).vtxInPixels{idxVtx}(:,1),roi_defs(idxRoi).vtxInPixels{idxVtx}(:,2))];
           else
               roi_defs(idxRoi).vtxInUm{idxVtx} = [];
           end
       end
       roi_defs(idxRoi).roiVolumeInUm3 = sum(temp)*volPerPixel;
   end
   
   % when save, save two new fields (posInUm, xyradiusInUm) in roi_defs, 
   %    and also save header and pixelPerUm
   save([pathname thefilename],'header','roi_defs', 'pixelPerUm', '-mat');

function volPerPixel = pixelVol(pixelPerUm)
    temp = 1./pixelPerUm;
    volPerPixel = prod(temp);
    
    
   
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
hAxes=handles.image_axes;
hImg=getappdata(handles.stack_vol_roi_rg, 'hImage');
[newclim, newcs] = imrangegui(get(hImg, 'CData'), get(hAxes, 'clim'), 0);
set(hAxes, 'clim', newclim);


% --- Executes on button press in deselect_pushbutton.
function deselect_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to deselect_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.stack_vol_roi_rg;
   roi_defs=getappdata(fig, 'roi_defs');
   setappdata(fig,'hLastROI',0);
   setappdata(fig,'lastROI',0);
   set(handles.selectRoi,'string','0');
   
   % set([roi_defs.hCircle roi_defs.hLabel], 'visible', 'off');
   % set([roi_defs.hCenter], 'visible', 'on');
   


function current_stack_Callback(hObject, eventdata, handles)
% hObject    handle to current_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_stack as text
%        str2double(get(hObject,'String')) returns contents of
%        current_stack as a double
fig=handles.stack_vol_roi_rg;
curStackIdx=getappdata(fig, 'curStackIdx');
newStackIdx=str2double(get(hObject,'String'));
stack=getappdata(fig, 'stack');
stackSize=stack.size;
if(newStackIdx<1 || newStackIdx>stackSize(4))
  newStackIdx = curStackIdx;
end
setappdata(fig,'curStackIdx', newStackIdx);
display_plot(fig);



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



function clim_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to clim_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clim_min_edit as text
%        str2double(get(hObject,'String')) returns contents of
%        clim_min_edit as a double
display_plot(handles.stack_vol_roi_rg)

% --- Executes during object creation, after setting all properties.
function clim_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clim_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function clim_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to clim_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clim_max_edit as text
%        str2double(get(hObject,'String')) returns contents of
%        clim_max_edit as a double
display_plot(handles.stack_vol_roi_rg)



% --- Executes during object creation, after setting all properties.
function clim_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clim_max_edit (see GCBO)
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
display_plot(handles.stack_vol_roi_rg)


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
display_plot(handles.stack_vol_roi_rg)


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
display_plot(handles.stack_vol_roi_rg)


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
   imsz = size(smm(:,:,1,1));
   setappdata(handles.stack_vol_roi_rg,'stack',smm);
   curFrameIdx=1;
   curStackIdx=1;
   setappdata(handles.stack_vol_roi_rg, 'headerFile', headerFile);
   h=imreadheader(headerFile);
   setappdata(handles.stack_vol_roi_rg, 'header', h);
   setappdata(handles.stack_vol_roi_rg, 'curFrameIdx', curFrameIdx);
   setappdata(handles.stack_vol_roi_rg, 'curStackIdx', curStackIdx);
   
   stim_lookup = h.stim_lookup;
   [onset_list,ustim] = find_stimulus_start(stim_lookup);
   smm_sz = smm.size;
   for i = 1:length(onset_list)
       onset_list{i}(onset_list{i}+4>smm_sz(4)) = [];
   end
   setappdata(handles.stack_vol_roi_rg, 'onset_list', onset_list);
   setappdata(handles.stack_vol_roi_rg, 'ustim', ustim);
   
   stim_labels = h.stim_labels;
   setappdata(handles.stack_vol_roi_rg, 'stim_labels', stim_labels);

   popup_labels{1} = 'grey';
   popup_labels(2:length(stim_labels)+1) = stim_labels;
   set(handles.display_popupmenu,'String',popup_labels)
   
   fig = handles.stack_vol_roi_rg;
   display_plot(fig);
   
   center = round(imsz/2);
   center(2,:) = center;
   show_zoom_images(handles, center, curFrameIdx)
   
   
   %--------------------------
  function display_plot(fig)
    
    handles = guidata(fig);
    
    display_mode = get(handles.display_popupmenu,'Value');
    
    if display_mode == 1 % grey mode
      
      set(handles.dfof_max_edit,'Visible','off');
      set(handles.dfof_min_edit,'Visible','off');
      set(handles.clim_max_edit,'Visible','off');
      set(handles.clim_min_edit,'Visible','off');
      set(handles.repeat_edit,'Visible','off');
      set(handles.repeat_text,'Visible','off');
      set(handles.repeatinc,'Visible','off');
      set(handles.repeatdec,'Visible','off');
      set(handles.clim_text,'Visible','off');
      set(handles.dfof_text,'Visible','off');
      set(handles.rg_mode_text,'Visible','off');
      set(handles.rg_display_popup,'Visible','off');
      
      set(handles.current_stack,'Visible','on');
      set(handles.stack_text,'Visible','on');
      set(handles.contrast_pushbutton,'Visible','on');
      
      display_grey(handles);
    else % rg
      
      set(handles.dfof_max_edit,'Visible','on');
      set(handles.dfof_min_edit,'Visible','on');
      set(handles.clim_max_edit,'Visible','on');
      set(handles.clim_min_edit,'Visible','on');
      set(handles.repeat_edit,'Visible','on');
      set(handles.repeat_text,'Visible','on');
      set(handles.repeatinc,'Visible','on');
      set(handles.repeatdec,'Visible','on');
      set(handles.clim_text,'Visible','on');
      set(handles.dfof_text,'Visible','on');
      set(handles.rg_mode_text,'Visible','on');
      set(handles.rg_display_popup,'Visible','on');
      
      set(handles.current_stack,'Visible','off');
      set(handles.stack_text,'Visible','off');
      set(handles.contrast_pushbutton,'Visible','off');
      
      rg_mode = getappdata(fig,'rg_mode');
      if ~isempty(strmatch(rg_mode,'raw'))
          set(handles.repeat_edit,'Visible','on');
          set(handles.repeat_text,'Visible','on');
          set(handles.repeatinc,'Visible','on');
          set(handles.repeatdec,'Visible','on');
      else
          set(handles.repeat_edit,'Visible','off');
          set(handles.repeat_text,'Visible','off');
          set(handles.repeatinc,'Visible','off');
          set(handles.repeatdec,'Visible','off');
      end
      
      display_rg(handles)
    end
    
function display_rg(handles)
  
  smm = getappdata(handles.stack_vol_roi_rg,'stack');
  frame_num = getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
  onset_list = getappdata(handles.stack_vol_roi_rg, 'onset_list');
  ustim = getappdata(handles.stack_vol_roi_rg, 'ustim');
  rg_mode = getappdata(handles.stack_vol_roi_rg,'rg_mode');
  
  set(handles.current_frame, 'string', num2str(frame_num));
  
  current_stim = get(handles.display_popupmenu,'Value') - 1;
  if ~isempty(strmatch(rg_mode,'raw','exact'))
      current_repeat = str2num(get(handles.repeat_edit,'String'));
  else
      current_repeat = 1:length(onset_list{current_stim});
  end
  
  clim(1) = str2num(get(handles.clim_min_edit,'String'));
  clim(2) = str2num(get(handles.clim_max_edit,'String'));
  
  if clim(1)>=clim(2) % taking care of bad user input
    clim(2) = clim(1) + 1;
    set(handles.clim_max_edit,'String',num2str(clim(2)));
  end
  
  dfof(1) = str2num(get(handles.dfof_min_edit,'String'));
  dfof(2) = str2num(get(handles.dfof_max_edit,'String'));
  
  if dfof(1)>=dfof(2) % taking care of bad user input
    dfof(2) = dfof(1)+0.1;
    set(handles.dfof_max_edit,'String',num2str(dfof(2)));
  end
  
  num_repeats = length(onset_list{1});
  if length(current_repeat) == 1
      if current_repeat < 0 || current_repeat > num_repeats % taking care of bad user input
          current_repeat = 1;
          set(handles.repeat_edit,'String',num2str(current_repeat));
      end
  end
  indx = onset_list{current_stim}(current_repeat)';
  setappdata(handles.stack_vol_roi_rg,'curStackIdx', indx);
  
  if ~isempty(strmatch(rg_mode,'zscore'))
      smm_sz = smm.size;
      imdiff = NaNs([smm_sz(1:2),length(onset_list{current_stim})]);
      imstd = NaNs(smm_sz(1:2));
      immean = imstd;
      imrgb = imstd;
      imbg = ones(smm_sz(1:2));
      for i = 1:length(onset_list{current_stim})
          imdiff(:,:,i) = mean(smm(:,:,frame_num,indx(i):indx(i)+4),4)...
                    - mean(smm(:,:,frame_num,indx(i)-5:indx(i)-1),4);
      end 
      imstd = std(imdiff,[],3);
      immean = mean(imdiff,3);
      imrgb = immean./(imstd+0.1*median(imdiff,3));
      
      options.clim_dfof = [-1.65 1.65];
      options.clim_raw = clim;
      imrgb = colorize_dfof({cat(3,imbg'+1,imrgb'+1)},options);
      imrgb = imrgb{1};
  else
      indices{1} = []; indices{2} = [];
      if length(current_repeat)==1
          indices{1} = indx-5:indx-1; % bg
          indices{2} = indx:indx+4; % stim
      else
          for i = 1:length(current_repeat)
              indices{1} = [indices{1} indx(i)-5:indx(i)-1];
              indices{2} = [indices{2} indx(i):indx(i)+4];
          end
      end
      im1 = mean(smm(:,:,frame_num, indices{1}), 4); %bg
      im2 = mean(smm(:,:,frame_num, indices{2}), 4); %stim
  
      im(:,:,1) = im1'; % this flip is to make it compatible with Zhongsheng's stack_roi code
      im(:,:,2) = im2';
  
      frames{1} = im;
      options.clim_dfof = dfof;
      %options.sigma = 2;
      options.clim_raw = clim;
  
      imrgb_tmp = colorize_dfof(frames, options);
      imrgb = imrgb_tmp{1};
  end
  
  axes(handles.image_axes);
  hImage=getappdata(handles.stack_vol_roi_rg, 'hImage');
  set(hImage, 'CData', imrgb);
  set(gca,'layer','bottom');
  
  show_frame_rois(handles);

    
function display_grey(handles)
  
   idxFrame = getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
   idxStack = getappdata(handles.stack_vol_roi_rg, 'curStackIdx');
  
   stack = getappdata(handles.stack_vol_roi_rg,'stack');
   stackSize = stack.size;
   
   imageData = stack(:,:,idxFrame,idxStack); % when get data in matlab, index is 1-based
   imageData  = nanmean(squeeze(imageData),3);
   imageData = imageData';

   imageW=size(imageData, 2);
   imageH=size(imageData, 1); 
   
   axsize = [imageW imageH];

   
   oldPosAxes=get(handles.image_axes, 'position');
   %set(handles.image_axes, 'position', [oldPosAxes(1:2) axsize]);
   set(handles.image_axes, 'position', oldPosAxes);
   oldPosFig=get(handles.stack_vol_roi_rg, 'position');
   
   colormap(gray(256));

   axes(handles.image_axes);
   hImage=getappdata(handles.stack_vol_roi_rg, 'hImage');
   if(~isempty(hImage) && ishandle(hImage))
      %set(fig, 'position', [oldPosFig(1:2) figsize]); % change fig size only
      set(hImage, 'CData', imageData);
   else
      %set(fig, 'position', [12 55 figsize]); % change both fig position and size
      hImage = imagesc(imageData,get(handles.image_axes, 'clim'));
      % hide the ticks and make sure the y dir is reverse
      set(handles.image_axes,'Visible','off', 'ydir', 'normal');
      
      setappdata(handles.stack_vol_roi_rg, 'hImage', hImage);
   end
   set(gca,'layer','bottom');
   
   set(handles.current_frame, 'string', num2str(idxFrame));
   set(handles.current_stack, 'string', num2str(idxStack));
   
   show_frame_rois(handles);

   
   
  function show_frame_rois(handles)
    
    % show rois if needed
   roi_defs=getappdata(handles.stack_vol_roi_rg, 'roi_defs');
      idxFrame = getappdata(handles.stack_vol_roi_rg, 'curFrameIdx');
   if(~isempty(roi_defs))
     % tt=findobj([roi_defs.hCircle roi_defs.hCenter roi_defs.hLabel], 'visible', 'on');
     % set(tt, 'visible', 'on');
     outlinesToShow = [];
     outlineCenter = [];
     idxCentersToShow = [];
     highlightCenter = [];
     highlightLabel = [];
     for i = 1:length(roi_defs)
         idxs = idxFrame-1:idxFrame+1; idxs(idxs<1)=[];
         tmp = intersect(find(roi_defs(i).hOutline), idxs);
         tmpc = intersect(find(roi_defs(i).hOutline), idxFrame);
         tmphi = intersect(roi_defs(i).centerInPixels(3), idxFrame);
         if ~isempty(tmp)
             outlinesToShow = unique([outlinesToShow roi_defs(i).hOutline(tmp)]);
             idxCentersToShow = [idxCentersToShow i];
         end
         if ~isempty(tmpc)
             outlineCenter = [outlineCenter roi_defs(i).hOutline(tmpc)];
         end
         if ~isempty(tmphi)
             highlightCenter = [highlightCenter roi_defs(i).hCenter];
             highlightLabel = [highlightLabel roi_defs(i).hLabel];
         end
     end
     centersToShow = [roi_defs(idxCentersToShow).hCenter];
     labelsToShow = [roi_defs(idxCentersToShow).hLabel];
     
     for hi=make_vector([outlinesToShow centersToShow labelsToShow], 'row')
         set(hi,'visible','on');
     end
     for hi=make_vector([centersToShow labelsToShow], 'row')
         set(hi,'color', [0.9 0.5 0.9]);
     end
     for hi=make_vector([highlightCenter highlightLabel], 'row')
         set(hi,'color', [0.7 0 0.7]);
     end
     set(outlinesToShow,'edgecolor',[0.6 0.6 1],'edgealpha',0.3);
     set(outlineCenter,'edgecolor','b','edgealpha',1);
     
     for idx=make_vector(idxCentersToShow, 'row')
       % apply the tform
       roi=roi_defs(idx);
       fPlotRois(handles.image_axes, roi); % move the position
       
       if(roi_defs(idx).showOutline)
         set([roi_defs(idx).hCenter roi_defs(idx).hLabel roi_defs(idx).hOutline(find(roi_defs(idx).hOutline))], 'visible', 'on' );
       else
         set([roi_defs(idx).hCenter roi_defs(idx).hLabel roi_defs(idx).hOutline(find(roi_defs(idx).hOutline))], 'visible', 'off' );
       end
     end
     idxRoisToHide=setdiff(1:length(roi_defs), idxCentersToShow);
     if(~isempty(idxRoisToHide))
         for i = idxRoisToHide
             set([roi_defs(i).hCenter roi_defs(i).hLabel roi_defs(i).hOutline(find(roi_defs(i).hOutline))], 'visible', 'off')
         end;
     end
   end
    


% --- Executes on selection change in display_popupmenu.
function display_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to display_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns display_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from display_popupmenu

display_plot(handles.stack_vol_roi_rg)

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
   display_plot(fig);

function onDecFrameIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx')-1;
   if(curFrameIdx==0)
      curFrameIdx=1;
   end
   setappdata(fig,'curFrameIdx', curFrameIdx);
   display_plot(fig);

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
   display_plot(fig);

   
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
   display_plot(fig);

  
   

function onDeleteRoi(sender, event)
   fig=get_parent_fig(sender);
   handles=guidata(fig);
   roi_defs=getappdata(fig, 'roi_defs');
   visibleRois=findobj([roi_defs.hOutline], 'visible', 'on');
   indicesToDel=[];
   if isappdata(fig,'hLastROI')
       hLastROI = getappdata(fig,'hLastROI');
   else
       hLastROI = [];
   end
   if ~isempty(hLastROI)
       if ishandle(hLastROI)
           hroi = intersect(hLastROI,visibleRois);
           if ~isempty(hroi)
               %            for hroi=visibleRois' % for's loop set must be row vector
               %                idx=find([roi_defs.hCenter]==hroi);
               %                indicesToDel(end+1)=idx;
               %                delete(roi_defs(idx).hOutline);
               %                delete(roi_defs(idx).hCenter);
               %                delete(roi_defs(idx).hLabel);
               %            end
               %        else
               idx=find([roi_defs.hCenter]==hroi);
               indicesToDel=idx;
               for i = 1:length(idx)
                   tmp = getappdata(roi_defs(idx(i)).hCenter);
                   delete(tmp.hOutline(find([tmp.hOutline])));
                   delete(tmp.hLabel);
                   delete(tmp.hCenter);
                   status_message = ['deleted ROI #' num2str(roi_defs(idx(i)).label)];
                   setappdata(fig,'status_message',status_message);
                   update_status_message(fig);
               end
           end
       end
   end
   roi_defs(indicesToDel)=[];
   setappdata(fig, 'hLastROI',[]);
   setappdata(fig, 'lastROI',0);
   setappdata(fig, 'roi_defs', roi_defs);

   
   % --- Executes on button press in zoom_radiobutton.
  function zoom_radiobutton_Callback(hObject, eventdata, handles)
    % hObject    handle to zoom_radiobutton (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Hint: get(hObject,'Value') returns toggle state of zoom_radiobutton
    show_zoom = get(hObject,'Value');
    setappdata(handles.stack_vol_roi_rg,'show_zoom',show_zoom);
    zoom1handle=getappdata(handles.stack_vol_roi_rg, 'zoom1handle');
    zoom2handle=getappdata(handles.stack_vol_roi_rg, 'zoom2handle');
    zoom3handle=getappdata(handles.stack_vol_roi_rg, 'zoom3handle');
    zoom4handle=getappdata(handles.stack_vol_roi_rg, 'zoom4handle');
    zoom5handle=getappdata(handles.stack_vol_roi_rg, 'zoom5handle');
    if ~show_zoom
      set(zoom1handle,'Visible','off');
      set(zoom2handle,'Visible','off');
      set(zoom3handle,'Visible','off');
      set(zoom4handle,'Visible','off');
      set(zoom5handle,'Visible','off');
      set(handles.zoom_text,'Visible','off');
      set(handles.zoom_size_edit,'Visible','off');
    else
      set(zoom1handle,'Visible','on');
      set(zoom2handle,'Visible','on');
      set(zoom3handle,'Visible','on');
      set(zoom4handle,'Visible','on');
      set(zoom5handle,'Visible','on');
      set(handles.zoom_text,'Visible','on');
      set(handles.zoom_size_edit,'Visible','on');
    end



function zoom_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoom_size_edit as text
%        str2double(get(hObject,'String')) returns contents of zoom_size_edit as a double
rmappdata(handles.stack_vol_roi_rg, 'zoom1handle');

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
    smm = getappdata(handles.stack_vol_roi_rg,'stack');
    smmsz=smm.size;
    
    posInPixel = round(posInPixel);
    zoom_size = str2double(get(handles.zoom_size_edit, 'String'));
    
    if zoom_size < 1 || zoom_size>min((smmsz(1:2)/2)-1) % taking care of bad user input
      zoom_size = 50; % default value
      set(handles.zoom_size_edit,'String',num2str(zoom_size))
    end
    
    axsize = [zoom_size*2+1 zoom_size*2+1];
    
    display_mode = get(handles.display_popupmenu,'Value');
    
    
    smmsz(1) = min(smmsz(1:2));
    smmsz(2) = min(smmsz(1:2));
    
    xsnip = posInPixel(1)-zoom_size:posInPixel(1)+zoom_size;
    xsnip(xsnip<1) = 1;
    xsnip(xsnip>smmsz(1)) = smmsz(1);
    
    ysnip = posInPixel(2)-zoom_size:posInPixel(2)+zoom_size;
    ysnip(ysnip<1) = 1;
    ysnip(ysnip>smmsz(2)) = smmsz(2);
    
    show_frames = curFrameIdx-2:curFrameIdx+2;
    show_frames(show_frames<1) = 1;
    show_frames(show_frames>smmsz(3)) = smmsz(3);
    
    
    if display_mode == 1 % show grey scale
      idxStack = getappdata(handles.stack_vol_roi_rg, 'curStackIdx');  
      for indx = 1:5
        imageData{indx} = smm(xsnip,ysnip,show_frames(indx),idxStack)'; % when get data in matlab, index is 1-based  
      end
      colormap(gray(256));
      
    else % display mode in stim colorize dfof
      onset_list = getappdata(handles.stack_vol_roi_rg, 'onset_list');
      ustim = getappdata(handles.stack_vol_roi_rg, 'ustim');
      current_stim = get(handles.display_popupmenu,'Value') - 1;
      current_repeat = str2num(get(handles.repeat_edit,'String'));
      
      clim(1) = str2num(get(handles.clim_min_edit,'String'));
      clim(2) = str2num(get(handles.clim_max_edit,'String'));
      
      dfof(1) = str2num(get(handles.dfof_min_edit,'String'));
      dfof(2) = str2num(get(handles.dfof_max_edit,'String'));
      
      indx = onset_list{current_stim}(current_repeat);
      setappdata(handles.stack_vol_roi_rg,'curStackIdx', indx);
      
      im1 = mean(smm(xsnip,ysnip,show_frames, indx-5:indx-1), 4); %bg
      im2 = mean(smm(xsnip,ysnip,show_frames, indx:indx+4), 4); %stim
      
      for indx = 1:5
        im(:,:,1) = im1(:,:,indx)'; 
        im(:,:,2) = im2(:,:,indx)';
        frames{indx} = im;
      end

      options.clim_dfof = dfof;
      %options.sigma = 2;
      options.clim_raw = clim;
      
      imageData = colorize_dfof(frames, options);
      
    end
    
    
    zoom1handle=getappdata(handles.stack_vol_roi_rg, 'zoom1handle');
    zoom2handle=getappdata(handles.stack_vol_roi_rg, 'zoom2handle');
    zoom3handle=getappdata(handles.stack_vol_roi_rg, 'zoom3handle');
    zoom4handle=getappdata(handles.stack_vol_roi_rg, 'zoom4handle');
    zoom5handle=getappdata(handles.stack_vol_roi_rg, 'zoom5handle');
    
    if(~isempty(zoom1handle) && ishandle(zoom1handle))
      set(zoom1handle, 'CData', imageData{1});
      set(zoom2handle, 'CData', imageData{2});
      set(zoom3handle, 'CData', imageData{3});
      set(zoom4handle, 'CData', imageData{4});
      set(zoom5handle, 'CData', imageData{5});
      
    else
      axes(handles.zoom1_axes);
      zoom1handle = imagesc(imageData{1},get(handles.image_axes, 'clim'));
      
      axes(handles.zoom2_axes);
      zoom2handle = imagesc(imageData{2},get(handles.image_axes, 'clim'));
      
      axes(handles.zoom3_axes);
      zoom3handle = imagesc(imageData{3},get(handles.image_axes, 'clim'));
      
      axes(handles.zoom4_axes);
      zoom4handle = imagesc(imageData{4},get(handles.image_axes, 'clim'));
      
      axes(handles.zoom5_axes);
      zoom5handle = imagesc(imageData{5},get(handles.image_axes, 'clim'));
      
      
      % hide the ticks and make sure the y dir is reverse
      set(handles.zoom1_axes,'Visible','off', 'ydir', 'normal');
      set(handles.zoom2_axes,'Visible','off', 'ydir', 'normal');
      set(handles.zoom3_axes,'Visible','off', 'ydir', 'normal');
      set(handles.zoom4_axes,'Visible','off', 'ydir', 'normal');
      set(handles.zoom5_axes,'Visible','off', 'ydir', 'normal');
      
      setappdata(handles.stack_vol_roi_rg, 'zoom1handle', zoom1handle);
      setappdata(handles.stack_vol_roi_rg, 'zoom2handle', zoom2handle);
      setappdata(handles.stack_vol_roi_rg, 'zoom3handle', zoom3handle);
      setappdata(handles.stack_vol_roi_rg, 'zoom4handle', zoom4handle);
      setappdata(handles.stack_vol_roi_rg, 'zoom5handle', zoom5handle);
    end
    
    
    



function merge1_Callback(hObject, eventdata, handles)
% hObject    handle to merge1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of merge1 as text
%        str2double(get(hObject,'String')) returns contents of merge1 as a double
tmp = guidata(hObject);
fig = tmp.stack_vol_roi_rg; 
thisval = str2num(get(tmp.merge1,'string')); clear tmp;
if isappdata(fig,'tomerge')
    tomerge = getappdata(fig,'tomerge');
else
    tomerge = zeros([1 2]);
end
tomerge(1) = thisval;
setappdata(fig,'tomerge',tomerge);

% --- Executes during object creation, after setting all properties.
function merge1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to merge1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function merge2_Callback(hObject, eventdata, handles)
% hObject    handle to merge2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of merge2 as text
%        str2double(get(hObject,'String')) returns contents of merge2 as a double
tmp = guidata(hObject);
fig = tmp.stack_vol_roi_rg; 
thisval = str2num(get(tmp.merge2,'string')); clear tmp;
if isappdata(fig,'tomerge')
    tomerge = getappdata(fig,'tomerge');
else
    tomerge = zeros([1 2]);
end
tomerge(2) = thisval;
setappdata(fig,'tomerge',tomerge);

% --- Executes during object creation, after setting all properties.
function merge2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to merge2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in merge_push.
function merge_push_Callback(hObject, eventdata, handles)
% hObject    handle to merge_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = guidata(hObject);
fig = tmp.stack_vol_roi_rg;
appdata = getappdata(fig);
roi = appdata.roi_defs;
mergeto = appdata.tomerge(1);
mergefrom = appdata.tomerge(2);
for i = 1:length(roi)
    match(i) = roi(i).label == mergefrom;
end
if sum(match) > 0
    [roi(match).label] = deal(mergeto);
    for i = 1:length(match)
        if match(i)
            set(roi(i).hLabel,'string',num2str(mergeto));
        end
    end
end
setappdata(fig,'roi_defs',roi);
show_frame_rois(handles);


% --- Executes on button press in frameadvance.
function frameadvance_Callback(hObject, eventdata, handles)
% hObject    handle to frameadvance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
onIncFrameIdx(hObject, eventdata);


% --- Executes on button press in framedecrease.
function framedecrease_Callback(hObject, eventdata, handles)
% hObject    handle to framedecrease (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
onDecFrameIdx(hObject, eventdata);


% --- Executes on button press in repeatinc.
function repeatinc_Callback(hObject, eventdata, handles)
% hObject    handle to repeatinc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = guidata(hObject);
fig = tmp.stack_vol_roi_rg;
rg_mode = getappdata(handles.stack_vol_roi_rg,'rg_mode');
if ~isempty(strmatch(rg_mode,'raw','exact'))
    onset_list = getappdata(handles.stack_vol_roi_rg, 'onset_list');
    num_repeats = length(onset_list{1});
    current_repeat = str2num(get(handles.repeat_edit,'String'))+1;
    if current_repeat > num_repeats % taking care of bad user input
        current_repeat = num_repeats;
    end
    set(handles.repeat_edit,'String',num2str(current_repeat));
    display_plot(fig);
else
    set(handles.repeat_edit,'String',str2num(get(handles.repeat_edit,'String')));
end

% --- Executes on button press in repeatdec.
function repeatdec_Callback(hObject, eventdata, handles)
% hObject    handle to repeatdec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = guidata(hObject);
fig = tmp.stack_vol_roi_rg;
rg_mode = getappdata(handles.stack_vol_roi_rg,'rg_mode');
if ~isempty(strmatch(rg_mode,'raw','exact'))
    onset_list = getappdata(handles.stack_vol_roi_rg, 'onset_list');
    num_repeats = length(onset_list{1});
    current_repeat = str2num(get(handles.repeat_edit,'String'))-1;
    if current_repeat < 1 % taking care of bad user input
        current_repeat = 1;
    end
    set(handles.repeat_edit,'String',num2str(current_repeat));
    display_plot(fig);
else
    set(handles.repeat_edit,'String',str2num(get(handles.repeat_edit,'String')));
end


% --- Executes during object creation, after setting all properties.
function image_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate image_axes


% --- Executes on selection change in rg_display_popup.
function rg_display_popup_Callback(hObject, eventdata, handles)
% hObject    handle to rg_display_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choices = get(hObject,'string');
setappdata(handles.stack_vol_roi_rg,'rg_mode',choices{get(hObject,'value')});
display_plot(handles.stack_vol_roi_rg);

% --- Executes during object creation, after setting all properties.
function rg_display_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rg_display_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'string',{'raw','mean','zscore'});



function selectRoi_Callback(hObject, eventdata, handles)
% hObject    handle to selectRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig = get_parent_fig(hObject);
roi_defs = getappdata(fig,'roi_defs');
roiselect = str2num(get(hObject,'string'));
idxselect = find([roi_defs.label]==roiselect);
hselect = roi_defs(idxselect).hCenter;
if ~isempty(roiselect)
    setappdata(fig, 'lastROI',idxselect);
    setappdata(fig, 'hLastROI',hselect);
end


% --- Executes during object creation, after setting all properties.
function selectRoi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
