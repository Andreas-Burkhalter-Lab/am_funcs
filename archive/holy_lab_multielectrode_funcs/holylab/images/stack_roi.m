function varargout = stack_roi(varargin)
% STACK_ROI M-file for stack_roi.figDefStackRoi
%      STACK_ROI, by itself, creates a new STACK_ROI or raises the existing
%      singleton*.
%
%      H = STACK_ROI returns the handle to a new STACK_ROI or the handle to
%      the existing singleton*.
%
%      STACK_ROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STACK_ROI.M with the given input arguments.
%
%      STACK_ROI('Property','Value',...) creates a new STACK_ROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stack_roi_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stack_roi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% keybinding:
%    up/down arrow: change frame # 
%    left/right arrow: change roi radius
%    DEL/d: delete selected ROIs
%    right-drag on circle: change roi radius
%    left-drag on circle: move roi position
%    middle-click: define a new roi
%    left-click on center: select a roi
%    left-click on circle: deselect a roi
%    pageup/down: change stack #

% Edit the above text to modify the response to help stack_roi

% Last Modified by GUIDE v2.5 03-Oct-2006 16:31:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stack_roi_OpeningFcn, ...
                   'gui_OutputFcn',  @stack_roi_OutputFcn, ...
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


% --- Executes just before stack_roi is made visible.
function stack_roi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stack_roi (see VARARGIN)

% Choose default command line output for stack_roi
handles.output = hObject;


% % Set up blank images
% handles.hImage = image('parent',handles.axesImage);
% set(handles.hImage, ...
%    'CDataMapping','scaled'); %,'EraseMode','none'); % TODO: erase mode?
set(handles.axesImage, ...
    'CLim',[500 3000]);
% colormap(gray(256));

fig=hObject;
bind_shortcut(fig, 'uparrow', @onIncFrameIdx);
bind_shortcut(fig, 'downarrow', @onDecFrameIdx);
bind_shortcut(handles.btnOpenFile, 'uparrow', @onIncFrameIdx);
bind_shortcut(handles.btnOpenFile, 'downarrow', @onDecFrameIdx);
bind_shortcut(fig, 'rightarrow', @onIncRoiRadius);
bind_shortcut(fig, 'leftarrow', @onDecRoiRadius);
bind_shortcut(fig, 'd', @onDeleteRoi);
bind_shortcut(fig, 'delete', @onDeleteRoi);
bind_shortcut(fig, 'pageup', @onIncStackIdx);
bind_shortcut(fig, 'pagedown', @onDecStackIdx);


install_mouse_event_handler(handles.axesImage, 'up', @onMouseUpOnAxesImage);

% Update handles structure
guidata(hObject, handles);

set(handles.editRoiRadius, 'string', 5);
set([handles.btnOpenFile], 'visible', 'off'); % TODO: temp disabled

% TODO: handle the callback (i.e. when user press ENTER)
set([handles.editCurZpos], 'enable', 'inactive');

if(nargin == 4 && ~isempty(varargin))
   stackRoiParams=varargin{1};
   fileToOpen=stackRoiParams.filename;
   openFile(handles, fileToOpen);
else
   error('Usage: stack_roi(struct(''filename'', headerFile))');
end


% UIWAIT makes stack_roi wait for user response (see UIRESUME)
% uiwait(handles.figDefStackRoi);

function pixel=um2pixel(handles, um, dim)
   h=getappdata(handles.figDefStackRoi, 'header');
   if(dim=='x' || dim=='y')
      if(h.um_per_pixel_xy==-1)
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
   radInUm=str2num(get(handles.editRoiRadius, 'String'));
   
% plot or update circle/center/label    
function hroi=fPlotRois(hAxes, roi)
   %axes(hAxes); % NOTE: or another way, set parent property in line()
   
   npts = 100;
   th = linspace(0,2*pi,npts);
   x = roi.posInPixel(1) + roi.xyradiusInPixel*cos(th);
   y = roi.posInPixel(2) + roi.xyradiusInPixel*sin(th);
   % the circle
   if(isfield(roi, 'hCircle'))
      % update data only
      set(roi.hCircle, 'XData', x, 'YData', y);
      hroi(1)=roi.hCircle;
   else
      hroi(1) = line(x,y, 'parent',hAxes,'linewidth', 1, 'color', 'blue');
   end
   % the circle center:
   if(isfield(roi, 'hCenter'))
      % update data only
      set(roi.hCenter, 'XData',roi.posInPixel(1), 'YData',roi.posInPixel(2));  
      hroi(2)=roi.hCenter;
   else
      hroi(2) = line(roi.posInPixel(1), roi.posInPixel(2),  'parent',hAxes,'linewidth', 2, 'marker', '.', 'color', 'red');
   end
   % the label:
   if(isfield(roi, 'hLabel'))
      % update data only
      set(roi.hLabel, 'Position', [roi.posInPixel(1)+roi.xyradiusInPixel, roi.posInPixel(2)]);
      hroi(3)=roi.hLabel;
   else
      hroi(3) = text(roi.posInPixel(1)+roi.xyradiusInPixel, roi.posInPixel(2), num2str(roi.label));
      set(hroi(3), 'color', 'r', 'FontSize', 14); % text color is red
   end

   
function addRoi(handles, roi)
% PRE: 
%    roi: a ROI that has 3 fields: 
%                       posInPixel, xyradiusInPixel, label % field 1,2,3
% POST:
%    fig's appdata roi_defs is updated;
%    circle/center/label objects are created
      hRois=fPlotRois(handles.axesImage, roi);
      roi.hCircle=hRois(1);                                % field 4
      roi.hCenter=hRois(2);                                % field 5 
      roi.hLabel=hRois(3);                                 % field 6
      roi.showCircle=1;                                    % field 7 of 7
      
      roi_defs=getappdata(handles.figDefStackRoi, 'roi_defs');
      roi_defs=[roi_defs roi];
      setappdata(handles.figDefStackRoi, 'roi_defs', roi_defs);
      
      % drag_line(roi.hCenter, struct('type', 'both'));
      set(roi.hCenter, 'ButtonDownFcn', @onClickCenter);
      setappdata(roi.hCenter, 'hCircle', roi.hCircle);
      setappdata(roi.hCircle, 'hCenter', roi.hCenter);
      setappdata(roi.hCenter, 'hLabel', roi.hLabel);
      setappdata(roi.hCircle, 'hLabel', roi.hLabel);
      drag_circle(roi.hCircle, struct('onDragStart', @onDragRoiStart, 'onDragDone', @onDragRoiDone));
      
      set(roi.hCenter, 'visible', 'off');
   
   
function result=onMouseUpOnAxesImage(sender, event_args)
   handles=guidata(sender);
   fig=handles.figDefStackRoi;
   isInDefTformMode=getappdata(fig, 'isInDefTformMode');

   if(is_button_down(sender, 'middle')) % sender: the axes
      if(isInDefTformMode)
	 return;
      end

      % add a new roi
      posInPixel=get(sender, 'currentpoint');
      curFrameIdx=getappdata(handles.figDefStackRoi, 'curFrameIdx');
      roi.posInPixel=[posInPixel(1,1:2) curFrameIdx];      % field 1
      radius=getRadius(handles); 
      roi.xyradiusInPixel=um2pixel(handles, radius, 'x');  % field 2

      roi_defs=getappdata(handles.figDefStackRoi, 'roi_defs');
      if(isempty(roi_defs))
	 existentLabels=[];
      else
	 existentLabels=[roi_defs.label];
      end
      roi.label=find_first_slot(existentLabels);           % field 3

      addRoi(handles, roi);
      
      result=0;
   else
      % show_popup_menu(sender);
      
      result=0;
   end

function updateShowCircleField(hCircle)
   fig=get_parent_fig(hCircle);
   roi_defs=getappdata(fig, 'roi_defs');
   idx=find([roi_defs.hCircle]==hCircle);
   isVisible=strcmp(get(hCircle, 'visible'), 'on');
   roi_defs(idx).showCircle=isVisible;
   setappdata(fig, 'roi_defs', roi_defs);
   
function onClickCenter(sender, eventdata) % sender: the center
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
   
   
function isContinue=onDragRoiStart(sender, event_args)
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
function varargout = stack_roi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   isInDefTformMode=getappdata(fig, 'isInDefTformMode');
   if(isInDefTformMode)
      return;
   end
   
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
   fig=handles.figDefStackRoi;
   
   % delete old rois and related objects
   roi_defs=getappdata(fig, 'roi_defs');
   if(~isempty(roi_defs))
      set([roi_defs.hCircle], 'visible', 'on');
      onDeleteRoi(hObject, []);
   end
   
   % tform_info=getappdata(fig, 'tform_info');
   
   
   % TODO: like calc_roi_inten(), here should check if header is match and may re-calc positions
   
   roi_defs=tt.roi_defs;
   
   if(isfield(tt, 'tform_info'))
      tform_info=tt.tform_info;
   else
      % this fields will be saved 
      tform_info.roi_def_time=0; % ROIs defined at stack 0
      tform_info.tform_def_times=[];
      tform_info.tforms=[];
   end
   tform_info.cur_tform_def_time=[]; % this cur_ field are for gui only, won't be saved in file
   % tform_info.cur_tform=[];
   setappdata(fig, 'tform_info', tform_info);
   
   isInDefTformMode=0;
   setappdata(fig, 'isInDefTformMode', isInDefTformMode);
   
   % remove fields
   roi_defs=rmfield(roi_defs, {'posInUm', 'xyradiusInUm'});
   
   for idx=1:length(roi_defs)
      addRoi(handles, roi_defs(idx));
   end


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   headerFile=getappdata(handles.figDefStackRoi, 'headerFile');
   header=getappdata(handles.figDefStackRoi, 'header');
   [pathstr,name,ext] = fileparts(headerFile);
   [thefilename,pathname] = uiputfile([pathstr name '.roidef'],'Save data to file...');
   roi_defs=getappdata(handles.figDefStackRoi, 'roi_defs');
   
   roi_defs=rmfield(roi_defs, {'hCircle','hCenter','hLabel', 'showCircle'});
   
   % convert pixel to um
   pixelPerUm(1)=um2pixel(handles, 1, 'x');
   pixelPerUm(2)=pixelPerUm(1);
   pixelPerUm(3)=um2pixel(handles, 1, 'z');
   for idxRoi=1:length(roi_defs)
      roi_defs(idxRoi).posInUm=roi_defs(idxRoi).posInPixel./pixelPerUm; % NOTE: this is relative pos to piezo start pos
      roi_defs(idxRoi).xyradiusInUm=roi_defs(idxRoi).xyradiusInPixel/pixelPerUm(1);
   end

   % tform infos
   tform_info=getappdata(handles.figDefStackRoi, 'tform_info');
   
   % when save, save two new fields (posInUm, xyradiusInUm) in roi_defs, 
   %    and also save header and pixelPerUm
   save([pathname thefilename],'header','roi_defs', 'pixelPerUm', 'tform_info', '-mat');


% --- Executes on button press in btnDone.
function btnDone_Callback(hObject, eventdata, handles)
% hObject    handle to btnDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   close;
   

function display_ref_image(handles)
   stack = getappdata(handles.figDefStackRoi,'stack');
   stackSize=stack.size;
   tform_info=getappdata(handles.figDefStackRoi, 'tform_info');
   if(isempty(tform_info.tforms))
      idxStack=tform_info.roi_def_time;
   else
      sorted_tform_def_times=sort(tform_info.tform_def_times);
      tt=find(sorted_tform_def_times<tform_info.cur_tform_def_time);
      if(isempty(tt))
	 idxStack=tform_info.roi_def_time;
      else
	 idxStack=sorted_tform_def_times(tt(end));
      end
   end
   
   tform=interpolate3dTform(tform_info, idxStack);
   
   % now find out the selected ctrl ROI
   roi_defs=getappdata(handles.figDefStackRoi, 'roi_defs');
   selectedCtrlRoi=findobj([roi_defs.hCircle], 'color', 'red');
   if(length(selectedCtrlRoi)>1)
      msgbox('Please select only one contrl ROI');
      return;
   end
   
   if(isempty(selectedCtrlRoi))
      msgbox('Please select one contrl ROI');
      return;
   end
   
   selectedCtrlRoiIdx=find([roi_defs.hCircle]==selectedCtrlRoi);
   roi=roi_defs(selectedCtrlRoiIdx);
   roi.posInPixel=tformfwd(tform, roi.posInPixel); 
   roi.posInPixel(3)=max(min(round(roi.posInPixel(3)), stackSize(4)-1), 0);
   idxFrame=roi.posInPixel(3);
   
   % load data and draw image
   imageData = stack(:,:,idxFrame+1,idxStack+1); % when get data in matlab, index is 1-based
   imageData=squeeze(imageData);
   imageData=imageData';
   imageData=flipud(imageData);
   
   refFig=figure('name', 'reference image');
   colormap(gray(256));
   hImage = imagesc(imageData,get(handles.axesImage, 'clim'));
   % hide the ticks and make sure the y dir is reverse
   set(gca,'Visible','off', 'ydir', 'reverse');
   
   % now draw the ROI too
   roi=rmfield(roi, {'hCircle','hCenter','hLabel'});
   hroi=fPlotRois(gca, roi);
   set(hroi(1), 'lineStyle', ':');
   
   
% NOTE: indices here is 0-based
function display_frame(fig, idxStack, idxFrame)
   handles=guidata(fig);
   stack = getappdata(fig,'stack');
   stackSize=stack.size;
   tform_info=getappdata(fig, 'tform_info');
   
   imageData = stack(:,:,idxFrame+1,idxStack+1); % when get data in matlab, index is 1-based
   imageData=squeeze(imageData);
   imageData=imageData';
%    imageData=flipud(imageData);
   imageW=size(imageData, 2);
   imageH=size(imageData, 1); 
   
   axsize = [imageW imageH];
   tPosPanelCtrl=get(handles.uipanelCtrl, 'position');

   tMaxVisibleAxSize=get(0, 'screensize');
   tMaxVisibleAxSize=tMaxVisibleAxSize(3:4)-tPosPanelCtrl([3 2]) -[100 100] ;
   tRatioHtoW=imageH/imageW;
   if(tMaxVisibleAxSize(1)*tRatioHtoW<=tMaxVisibleAxSize(2))
      tMaxVisibleAxSize(2)=tMaxVisibleAxSize(1)*tRatioHtoW;
   else
      tMaxVisibleAxSize(1)=tMaxVisibleAxSize(2)/tRatioHtoW;
   end
   axsize=min([axsize; tMaxVisibleAxSize]);
   figsize = axsize + tPosPanelCtrl([3 2])+[20,10];
   if(figsize(2)<800)
      figsize(2)=800;
   end
   
   oldPosAxes=get(handles.axesImage, 'position');
   set(handles.axesImage, 'position', [oldPosAxes(1:2) axsize]);
   oldPosFig=get(fig, 'position');
   
   colormap(gray(256));

   axes(handles.axesImage);
   hImage=getappdata(fig, 'hImage');
   if(~isempty(hImage) && ishandle(hImage))
      set(fig, 'position', [oldPosFig(1:2) figsize]); % change fig size only
      set(hImage, 'CData', imageData);
   else
      set(fig, 'position', [12 55 figsize]); % change both fig position and size
      hImage = imagesc(imageData,get(handles.axesImage, 'clim'));
      % hide the ticks and make sure the y dir is reverse
      set(handles.axesImage,'Visible','off', 'ydir', 'normal');
      
      setappdata(fig, 'hImage', hImage);
   end
   
   set(handles.editCurZpos, 'string', num2str(idxFrame));
   set(handles.editCurStack, 'string', num2str(idxStack));
   
   % set(handles.hImage,'CData', );
   
   % show rois if needed
   isInDefTformMode=getappdata(fig, 'isInDefTformMode');
   roi_defs=getappdata(fig, 'roi_defs');
   if(~isempty(roi_defs))
      % tt=findobj([roi_defs.hCircle roi_defs.hCenter roi_defs.hLabel], 'visible', 'on');
      % set(tt, 'visible', 'on');
      tt=cat(1,roi_defs.posInPixel);
      idxRoisToShow=find(tt(:,3)==idxFrame);
      if(isInDefTformMode)
	 ctrlROIs=setdiff([roi_defs.hCircle]', findobj([roi_defs.hCircle], 'color', 'blue'));
	 [tt, idxCtrlRois, tt2]=intersect([roi_defs.hCircle], ctrlROIs);
	 idxRoisToShow=union(idxRoisToShow, idxCtrlRois);
      end
      for idx=make_vector(idxRoisToShow, 'row')
	 % apply the tform
	 roi=roi_defs(idx);
	 tform=interpolate3dTform(tform_info, idxStack);
	 roi.posInPixel=tformfwd(tform, roi.posInPixel); 
	 roi.posInPixel(3)=max(min(round(roi.posInPixel(3)), stackSize(4)-1), 0);
	 fPlotRois(handles.axesImage, roi); % move the position
	 
	 if(roi_defs(idx).showCircle)
	    set([roi_defs(idx).hCircle roi_defs(idx).hLabel], 'visible', 'on');
	    set(roi_defs(idx).hCenter, 'visible', 'off');
	 else
	    set([roi_defs(idx).hCircle roi_defs(idx).hLabel], 'visible', 'off');
	    set(roi_defs(idx).hCenter, 'visible', 'on');
	 end
      end
      idxRoisToHide=setdiff(1:length(roi_defs), idxRoisToShow);
      if(~isempty(idxRoisToHide))
	 set([roi_defs(idxRoisToHide).hCircle roi_defs(idxRoisToHide).hLabel ...
	      roi_defs(idxRoisToHide).hCenter], 'visible', 'off');
      end
   end

function openFile(handles, headerFile)
   smm = stackmm(headerFile); % stack memory map
   setappdata(handles.figDefStackRoi,'stack',smm);
   curFrameIdx=0;
   curStackIdx=0;
   isInDefTformMode=0;
   setappdata(handles.figDefStackRoi, 'headerFile', headerFile);
   h=imreadheader(headerFile);
   setappdata(handles.figDefStackRoi, 'header', h);
   setappdata(handles.figDefStackRoi, 'curFrameIdx', curFrameIdx);
   setappdata(handles.figDefStackRoi, 'curStackIdx', curStackIdx);
   setappdata(handles.figDefStackRoi, 'isInDefTformMode', isInDefTformMode);
   
   tform_info.roi_def_time=0;
   tform_info.tform_def_times=[];
   tform_info.tforms=[];
   tform_info.cur_tform_def_time=[];
   setappdata(handles.figDefStackRoi, 'tform_info', tform_info);
   
   display_frame(handles.figDefStackRoi, curStackIdx, curFrameIdx);
   
   
% --- Executes on button press in btnOpenFile.
function btnOpenFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   [basename,pathname] = uigetfile('*','Select header file');
   if(basename==0)
      return;
   end
   headerFile=[pathname filesep basename];
   % setappdata(handles.figDefStackRoi, 'file', headerFile);
   openFile(handles, headerFile);


function editRoiRadius_Callback(hObject, eventdata, handles)
% hObject    handle to editRoiRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRoiRadius as text
%        str2double(get(hObject,'String')) returns contents of editRoiRadius as a double


% --- Executes during object creation, after setting all properties.
function editRoiRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRoiRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCurZpos_Callback(hObject, eventdata, handles)
% hObject    handle to editCurZpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCurZpos as text
%        str2double(get(hObject,'String')) returns contents of editCurZpos as a double


% --- Executes during object creation, after setting all properties.
function editCurZpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCurZpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function onIncFrameIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx')+1;
   curStackIdx=getappdata(fig, 'curStackIdx');
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curFrameIdx==stackSize(3))
      curFrameIdx=curFrameIdx-1;
   end
   setappdata(fig,'curFrameIdx', curFrameIdx);
   display_frame(fig, curStackIdx, curFrameIdx);

function onDecFrameIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx')-1;
   curStackIdx=getappdata(fig, 'curStackIdx');
   stack=getappdata(fig, 'stack');
   if(curFrameIdx==-1)
      curFrameIdx=0;
   end
   setappdata(fig,'curFrameIdx', curFrameIdx);
   display_frame(fig, curStackIdx, curFrameIdx);

function onIncStackIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   curStackIdx=getappdata(fig, 'curStackIdx')+1;
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curStackIdx==stackSize(4))
      %curStackIdx=curStackIdx-1;      % stay at the last stack
      curStackIdx=0;                  % wrap around to the beginning
   end
   setappdata(fig,'curStackIdx', curStackIdx);
   display_frame(fig, curStackIdx, curFrameIdx);

   
function onDecStackIdx(sender, event)
   fig=get_parent_fig(sender);
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   curStackIdx=getappdata(fig, 'curStackIdx')-1;
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curStackIdx==-1)
      %curStackIdx=0;                % stack with the first stack
      curStackIdx = stackSize(4)-1;   % wrap around to the end
   end
   setappdata(fig,'curStackIdx', curStackIdx);
   display_frame(fig, curStackIdx, curFrameIdx);
   
   
   
function step=getRadiusStep(handles)
   step=str2num(get(handles.editRadiusStep, 'string'));   

   
function minRadius=getMinRadius(handles)
   % TODO:
   minRadius=2;
   
   
function changeRoiRadius(sender, isInc)
   fig=get_parent_fig(sender);
   handles=guidata(fig);
   roi_defs=getappdata(fig, 'roi_defs');
   visibleRois=findobj([roi_defs.hCircle], 'visible', 'on');
   radiusStep=getRadiusStep(handles);
   radiusStepInPixel=um2pixel(handles, radiusStep, 'x');
   minRadius=getMinRadius(handles);
   minRadiusInPixel=um2pixel(handles, minRadius, 'x');
   if(~isInc) 
      radiusStepInPixel=-radiusStepInPixel;
   end
   for hroi=make_vector(visibleRois,'row') % for's loop set must be row vector
      idx=find([roi_defs.hCircle]==hroi);
      roi_defs(idx).xyradiusInPixel=roi_defs(idx).xyradiusInPixel+radiusStepInPixel;
      roi_defs(idx).xyradiusInPixel=max(roi_defs(idx).xyradiusInPixel, minRadiusInPixel);
      fPlotRois(handles.axesImage, roi_defs(idx));
   end
   
   setappdata(fig, 'roi_defs', roi_defs);

   
function onIncRoiRadius(sender, event)
   changeRoiRadius(sender, 1);

   
function onDecRoiRadius(sender, event)
   changeRoiRadius(sender, 0);
   

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
   


% --- Executes on button press in btnChangeContrast.
function btnChangeContrast_Callback(hObject, eventdata, handles)
% hObject    handle to btnChangeContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hAxes=handles.axesImage;
   hImg=getappdata(handles.figDefStackRoi, 'hImage');
   [newclim, newcs] = imrangegui(get(hImg, 'CData'), get(hAxes, 'clim'), 0);
   set(hAxes, 'clim', newclim);

 


% --- Executes on button press in btnDeselectAll.
function btnDeselectAll_Callback(hObject, eventdata, handles)
% hObject    handle to btnDeselectAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   roi_defs=getappdata(fig, 'roi_defs');
   for idx=1:length(roi_defs)
      roi_defs(idx).showCircle=0;
   end
   setappdata(fig, 'roi_defs', roi_defs);
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   display_frame(fig, curStackIdx, curFrameIdx);
   
   % set([roi_defs.hCircle roi_defs.hLabel], 'visible', 'off');
   % set([roi_defs.hCenter], 'visible', 'on');
   


% --- Executes on button press in btnSetAsRoiDefTime.
function btnSetAsRoiDefTime_Callback(hObject, eventdata, handles)
% hObject    handle to btnSetAsRoiDefTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   curStackIdx=getappdata(fig, 'curStackIdx');
   tform_info=getappdata(fig, 'tform_info');
   tform_info.roi_def_time=curStackIdx;
   setappdata(fig, 'tform_info', tform_info);
   

% --- Executes on button press in btnSetAsTformTime.
function btnSetAsTformTime_Callback(hObject, eventdata, handles)
% hObject    handle to btnSetAsTformTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   curStackIdx=getappdata(fig, 'curStackIdx');
   tform_info=getappdata(fig, 'tform_info');
   tform_info.cur_tform_def_time=curStackIdx;
   setappdata(fig, 'tform_info', tform_info);
   

% --- Executes on button press in btnBeginDefTform.
function btnBeginDefTform_Callback(hObject, eventdata, handles)
% hObject    handle to btnBeginDefTform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   
   % ask user about current tform time
   tform_info=getappdata(fig, 'tform_info');
   headerFile=getappdata(fig, 'headerFile');
   t=stackviewer(headerFile, struct('calledByStackRoi', 1, 'existingTformTimes', ...
				    [tform_info.tform_def_times tform_info.roi_def_time]+1));
   if(~isempty(t))
      curStackIdx=t-1;
      tform_info.cur_tform_def_time=curStackIdx;
      setappdata(fig, 'tform_info', tform_info);
      setappdata(fig, 'curStackIdx', curStackIdx);
      curFrameIdx=getappdata(fig, 'curFrameIdx');
      display_frame(fig, curStackIdx, curFrameIdx);
   end
   
   isInDefTformMode=1;
   setappdata(fig, 'isInDefTformMode', isInDefTformMode);
   
   % backup the roi definitions b/c when drag roi, roi def changes.
   roi_defs=getappdata(fig, 'roi_defs');
   backup_roi_defs=roi_defs;
   setappdata(fig, 'backup_roi_defs', backup_roi_defs); 
   
   set(fig, 'name', 'Define Stack Tforms');
   

% --- Executes on button press in btnEndDefTform.
function btnEndDefTform_Callback(hObject, eventdata, handles)
% hObject    handle to btnEndDefTform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   roi_defs=getappdata(fig, 'roi_defs');
   selectedCtrlROIs=findobj([roi_defs.hCircle], 'color', 'red');
   deselectedCtrlROIs=findobj([roi_defs.hCircle], 'color', 'green');
   ctrlROIs=[selectedCtrlROIs;deselectedCtrlROIs]; % ctrlROIs is a col vector
   [tt, indices, tt2]=intersect([roi_defs.hCircle], ctrlROIs);
   roiNewPos=cat(1, roi_defs(indices).posInPixel);
   
   backup_roi_defs=getappdata(fig, 'backup_roi_defs');
   roiOldPos=cat(1, backup_roi_defs(indices).posInPixel);
   
   % calc tform (from roi position at roi def time to what at roi tform time)
   tform=cp2tform_3d_affine(roiOldPos, roiNewPos);
   if(isempty(tform))
      msgbox('failed to calc the tform from control ROIs');
      return;
   end
   
   % save the tform info
   tform_info=getappdata(fig, 'tform_info');
   tform_info.tform_def_times(end+1)=tform_info.cur_tform_def_time;
   tform_info.tforms=[tform_info.tforms  tform];
   setappdata(fig, 'tform_info', tform_info);
   
   switch2DefRoiMod(fig);
   
   
   
function switch2DefRoiMod(fig)
   % restore roi_defs
   backup_roi_defs=getappdata(fig, 'backup_roi_defs');
   setappdata(fig, 'roi_defs', backup_roi_defs);
   
   % change mode:
   isInDefTformMode=0;
   setappdata(fig, 'isInDefTformMode', isInDefTformMode);
   
   % change control rois' color
   roi_defs=getappdata(fig, 'roi_defs');
   selectedCtrlROIs=findobj([roi_defs.hCircle], 'color', 'red');
   deselectedCtrlROIs=findobj([roi_defs.hCircle], 'color', 'green');
   ctrlROIs=[selectedCtrlROIs;deselectedCtrlROIs]; % ctrlROIs is a col vector
   set(ctrlROIs, 'color', 'blue');
   
   % redraw the screen
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   curStackIdx=getappdata(fig, 'curStackIdx');
   display_frame(fig, curStackIdx, curFrameIdx);
   
   set(fig, 'name', 'Define Stack ROIs');
   
   
   
% --- Executes on button press in btnCancelDefTform.
function btnCancelDefTform_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancelDefTform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   switch2DefRoiMod(handles.figDefStackRoi);
   
   

% --- Executes on button press in btnShowRefImage.
function btnShowRefImage_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowRefImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   isInDefTformMode=getappdata(handles.figDefStackRoi, 'isInDefTformMode');
   if(~isInDefTformMode) 
      return;
   end
   
   display_ref_image(handles);
   
   

% --- Executes on button press in btnChgCurCtrlRoiZpos.
function btnChgCurCtrlRoiZpos_Callback(hObject, eventdata, handles)
% hObject    handle to btnChgCurCtrlRoiZpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figDefStackRoi;
   roi_defs=getappdata(fig, 'roi_defs');
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   selectedCtrlROIs=findobj([roi_defs.hCircle], 'color', 'red');
   for idxRoi=1:length(roi_defs)
      if(any(roi_defs(idxRoi).hCircle==selectedCtrlROIs))
	 roi_defs(idxRoi).posInPixel(3)=curFrameIdx;
      end
   end
   setappdata(fig, 'roi_defs', roi_defs);
   



function editCurStack_Callback(hObject, eventdata, handles)
% hObject    handle to editCurStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCurStack as text
%        str2double(get(hObject,'String')) returns contents of editCurStack as a double
   fig=get_parent_fig(hObject);
   curFrameIdx=getappdata(fig, 'curFrameIdx');
   %curStackIdx=round(getappdata(fig, 'curStackIdx'));
   curStackIdx = str2num(get(hObject,'String'));
   stack=getappdata(fig, 'stack');
   stackSize=stack.size;
   if(curStackIdx>=stackSize(4))
      curStackIdx=curStackIdx-1;      % stay at the last stack
   elseif (curStackIdx < 0)
      curStackIdx=0;                  % stay at the beginning
   end
   setappdata(fig,'curStackIdx', curStackIdx);
   display_frame(fig, curStackIdx, curFrameIdx);


