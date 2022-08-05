function varargout = manipulate_rois(varargin)
% MANIPULATE_ROIS M-file for manipulate_rois.fig
%      MANIPULATE_ROIS, by itself, creates a new MANIPULATE_ROIS or raises the existing
%      singleton*.
%
%      H = MANIPULATE_ROIS returns the handle to a new MANIPULATE_ROIS or the handle to
%      the existing singleton*.
%
%      MANIPULATE_ROIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANIPULATE_ROIS.M with the given input arguments.
%
%      MANIPULATE_ROIS('Property','Value',...) creates a new MANIPULATE_ROIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manipulate_rois_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manipulate_rois_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manipulate_rois

% Last Modified by GUIDE v2.5 11-Aug-2005 12:17:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manipulate_rois_OpeningFcn, ...
                   'gui_OutputFcn',  @manipulate_rois_OutputFcn, ...
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


% --- Executes just before manipulate_rois is made visible.
function manipulate_rois_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manipulate_rois (see VARARGIN)
   args=varargin{1};
   setappdata(hObject, 'args', args);
   
   set(handles.editNumOfRoisPerCol, 'string', '1');
   
   init(hObject);

   % Choose default command line output for manipulate_rois
   handles.output = hObject;

   % Update handles structure
   guidata(hObject, handles);

   % UIWAIT makes manipulate_rois wait for user response (see UIRESUME)
   % uiwait(handles.figureRois);


% --- Outputs from this function are returned to the command line.
function varargout = manipulate_rois_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


  
% todo: pick next line's color by: c = uisetcolor([1 0 0], 'DialogTitle');

% note: in this func, fig's guidata is readonly
function init(fig)
   handles=guidata(fig);
   args=getappdata(fig, 'args');
   index=getappdata(args.axes, 'index');
   figSummary=getappdata(get_parent_fig(args.axes), 'figSummary');
   ip=getappdata(figSummary, 'ip');
   imageData=imphysfetch(ip.dfdraw(index));
   imageTime=ip.dfdraw(index).stacknum;
   
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
   
   oldPosAxes=get(handles.axesImage, 'position');
   set(handles.axesImage, 'position', [oldPosAxes(1:2) axsize]);
   oldPosFig=get(fig, 'position');
   set(fig, 'position', [12 55 figsize]);
   
   if(args.options_imstimviewdf.invert)
      colormap(1-gray(256))
   else
      colormap(gray(256))
   end

   axes(handles.axesImage);
   hImage = imagesc(imageData,get(args.axes, 'clim'));
   set(handles.axesImage,'Visible','off')

   setappdata(fig, 'hImage', hImage);
   setappdata(fig, 'imageTime', imageTime);
   
   set(fig, 'name', ['Manage ROIs (frame=' num2str(imageTime) ': valve=' ...
      num2str(getappdata(get_parent_fig(args.axes), 'valve')) ', trial=' ...
      num2str(getappdata(args.axes, 'trial')) ')' ]);
   
   oldRois=get_roi(figSummary, imageTime);
   if(isempty(oldRois))
      nOldRois=0;
   else
      nOldRois=length(oldRois.type);
   end
   
   % todo: maybe should move to get_roi()?
   for idx=1:nOldRois
      oldRois.action(idx)={'none'}; % stands for no change. Other values: 'modified', 'deleted', 'added'.
   end
   
   hRois=fPlotRois(fig, oldRois);
   for idx=1:length(hRois)
      drag_circle(hRois(idx), struct('onDragDone', @onDragRoiDone, 'onDragStart', @onDragRoiStart));
   end
   if(~isempty(hRois))
      oldRois.handle=hRois;
   end
   setappdata(fig, 'roi_defs', oldRois);
   
% --- Executes on button press in btnAddRoi.
function btnAddRoi_Callback(hObject, eventdata, handles)
% hObject    handle to btnAddRoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=get_parent_fig(hObject);
   roi_defs=getappdata(fig, 'roi_defs');
   if(isempty(roi_defs) || isempty(roi_defs.label) ) 
      nextLabel=1;
   else
      nextLabel=max(roi_defs.label)+1;
   end
   defaultRoi=[100 100 str2num(get(handles.editDefaultRoiSize, 'string'))]; % the values are before translation. todo: should use the value after translation
   tRoi.type='c';
   tRoi.label=nextLabel;
   tRoi.x=defaultRoi(1);
   tRoi.y=defaultRoi(2);
   tRoi.xyradius=defaultRoi(3);
   tRoi.action={'added'};
   tRoi.handle=[];

   % prepare for merging fields:
   if(~isempty(roi_defs) && ~isfield(roi_defs, 'action'))
      roi_defs.action={}; 
   end
   if(~isempty(roi_defs) && ~isfield(roi_defs, 'handle'))
      roi_defs.handle=[]; 
   end
   
   roi_defs=merge_roi_defs(roi_defs, tRoi);
   if(~isfield(roi_defs, 'tform'))
      roi_defs.tform=unit_tform;
   end
   tRoi.tform=roi_defs.tform;
   hRoi=fPlotRois(fig, tRoi);
   drag_circle(hRoi, struct('onDragDone', @onDragRoiDone, 'onDragStart', @onDragRoiStart));
   roi_defs.handle=[roi_defs.handle hRoi];
   setappdata(fig, 'roi_defs', roi_defs);

% --- Executes on button press in btnDelSelected.
function btnDelSelected_Callback(hObject, eventdata, handles)
% hObject    handle to btnDelSelected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   circlesToDel=findobj(handles.figureRois, 'type', 'line', 'marker', '*');
   deleteRois(handles.figureRois, circlesToDel);

% --- Executes on button press in btnDelAll.
function btnDelAll_Callback(hObject, eventdata, handles)
% hObject    handle to btnDelAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   circlesToDel=findobj(handles.figureRois, 'type', 'line');
   deleteRois(handles.figureRois, circlesToDel);

function result=unit_tform
   result=maketform('affine',[1 0 0; 0 1 0; 0 0 1]);

function hrois=fPlotRois(fig, rois)
   handles=guidata(fig);
   hrois = roiplot(handles.axesImage,rois);
   
function isContinue=onDragRoiStart(sender, event_args)
   isContinue=1;
   if(strcmp(event_args.event_type, 'resized'))
      resizing=1;
   else
      resizing=0;
   end
   if(is_modifier_down(sender, 'shift')) % sender: the line object
      isChangeAll=1;
   else
      isChangeAll=0;
   end
   % todo: should forbid resizing all?
   setappdata(sender, 'isChangeAll', isChangeAll); % note: drag_circle() isn't aware of "change all" so we use appdata here
   if(isChangeAll)
      fig=get_parent_fig(sender);
      circles=findobj(fig, 'type', 'line');
      set(circles, 'color', 'red');
   end

function onDragRoiDone(sender, event_args)
   fig=get_parent_fig(sender); % sender: the line
   isChangeAll=getappdata(sender, 'isChangeAll');
   
   % change circles' color back if necessary:
   if(isChangeAll)
      circles=findobj(fig, 'type', 'line');
      set(circles, 'color', 'blue');
   end
   
   roiLabel=getappdata(sender, 'label'); 
   rois=getappdata(fig, 'roi_defs');
   idx=find(rois.label==roiLabel);
   event_args.def(1:2)=tforminv(rois.tform, event_args.def(1:2));
   % todo: should round() new position here?
   if(isequal_roi([rois.x(idx) rois.y(idx) rois.xyradius(idx)], event_args.def))
      % if just mouse down/up
      oldMarker=get(sender, 'marker');
      if(strcmp(oldMarker, '*'))
         newMarker='none';
      else
         newMarker='*';
      end
      set(sender, 'marker', newMarker);
   else % else, roi def is changed
      if(strcmp(event_args.event_type, 'resized'))
         resized=1;
      else
         resized=0;
      end
      
      % here need distinguish "change tform" from "change a single ROI's def"
      if(isChangeAll) % if, change all ROIs (either tform or radius)
         if(resized)
            for idxEachRoi=1:length(rois.x)
               rois.xyraduis(idxEachRoi)=event_args.def(3);
               if(strcmp(rois.action{idxEachRoi}, 'none'))
                  rois.action{idxEachRoi}='modified';
               end
            end
         else % else, translation, so only need change tform
            delta=event_args.def(1:2)-[rois.x(idx) rois.y(idx)];
            % todo: should round() delta?
            rois.tform=translate_tform(rois.tform, delta);
            rois.tform_changed=1;
         end
      else % else, change a single ROI's def
         if(resized) 
            rois.xyradius(idx)=event_args.def(3);
            if(strcmp(rois.action{idx}, 'none'))
               rois.action{idx}='modified';
            end
         else % else, moved
            rois.x(idx)=event_args.def(1);
            rois.y(idx)=event_args.def(2);
            if(strcmp(rois.action{idx}, 'none'))
               rois.action{idx}='modified';
            end
         end
      end
      
      setappdata(fig, 'roi_defs', rois);
      
      % update ploting
      tt=fPlotRois(fig, rois);
   end 

   
% --- Executes on button press in btnApply.
function btnApply_Callback(hObject, eventdata, handles)
% hObject    handle to btnApply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=get_parent_fig(hObject);
   roi_defs=getappdata(fig, 'roi_defs');
   if(isempty(roi_defs) || isempty(roi_defs.label) ) 
      uiwait(msgbox('No change to apply.','Manage ROIs','modal'));
      return;
   end
   roi_defs_changed=0;
   for idx=1:length(roi_defs.label)
      if(~strcmp(roi_defs.action{idx}, 'none'))
         roi_defs_changed=1;
         break;
      end
   end
   if(roi_defs_changed || (isfield(roi_defs, 'tform_changed') && roi_defs.tform_changed) )
      % do nothing
   else
      uiwait(msgbox('No change to apply.','Manage ROIs','modal'));
      return;
   end
   
   args=getappdata(fig, 'args');
   figSummary=getappdata(get_parent_fig(args.axes), 'figSummary');
   imageTime=getappdata(fig, 'imageTime');
   roi_defs=set_roi(figSummary, imageTime, roi_defs); % will delete to-be-deleted rois
   setappdata(fig, 'roi_defs', roi_defs);
   
   
function result=isequal_roi(a,b)   
   result=isequal(round(a),round(b));
   
   
function deleteRois(fig, circlesToDel)
   if(isempty(circlesToDel)) return; end
   
   roi_defs=getappdata(fig, 'roi_defs');
   for idxRoiToDel=1:length(circlesToDel)
      hRoi=circlesToDel(idxRoiToDel);
      label=getappdata(hRoi, 'label');
      labelObject=getappdata(hRoi, 'labelObject');
      delete(labelObject); 
      delete(hRoi);
      idx=find(roi_defs.label==label);
      if(strcmp(roi_defs.action{idx},'added')) % if, delete newly added
         roi_defs=delete_roi_defs(roi_defs, label);
      else % else, delete orignal (then defer the action until "apply") 
         roi_defs.action{idx}='deleted';
      end
   end
   setappdata(fig, 'roi_defs', roi_defs); % save back the roi definitions
   
   
function newTform=translate_tform(oldTform, translation)   
   translationOld=oldTform.tdata.T(end, 1:2);
   translationNew=translationOld+translation;
   affine_matrix=[eye(2); translationNew];
   newTform=maketform('affine', affine_matrix);

% pre: fig: figureRois
% note: plot the intensities of current ROIs (not what saved in figSummary)
function plotRoiIntensities(fig)   
   roi_defs=getappdata(fig, 'roi_defs');

   args=getappdata(fig, 'args');
   % index=getappdata(args.axes, 'index');
   figSummary=getappdata(get_parent_fig(args.axes), 'figSummary');
   ip=getappdata(figSummary, 'ip');
   % imageData=imphysfetch(ip.dfdraw(index));
   % imageTime=ip.dfdraw(index).stacknum;
   
   valveFigures=getappdata(figSummary, 'valveFigures');
   validVlvFig=[];
   for idxValveFig=1:length(valveFigures)
      if(ishandle(valveFigures(idxValveFig)))
         validVlvFig=[validVlvFig valveFigures(idxValveFig)];
      end
   end
   valveFigures=validVlvFig;
   
   for idxValveFig=1:length(valveFigures)
      curVlvFig=valveFigures(idxValveFig);
      curValve=getappdata(curVlvFig, 'valve');
      % allAxes=findobj(curVlvFig, 'type', 'axes');
      allAxes=getappdata(curVlvFig, 'allAxes');
      % for curAxes=allAxes
      for idxAxes=1:length(allAxes)
         curAxes=allAxes(idxAxes);
         index=getappdata(curAxes, 'index');
         trials(idxValveFig).index(idxAxes)=index;
         trials(idxValveFig).trialNum(idxAxes)=getappdata(curAxes, 'trial');
         trials(idxValveFig).valve(idxAxes)=curValve;
         % im=imphysfetch(ip.dfcalc(index));
         % trials(idxValveFig).intensities(idxAxes)={roimeasure(im, roi_defs)};
      end
   end % for, each valve
   
   cellarrayRoi_defs=split_roi_defs(roi_defs); % 
   handles=guidata(fig);
   nRoiPerRow=str2num(get(handles.editNumOfRoisPerRow, 'string'));
   nRoiPerCol=str2num(get(handles.editNumOfRoisPerCol, 'string'));
   needNewFigure=1;
   allRoiPlotAxes=[];
   for idxRoi=1:length(cellarrayRoi_defs)
      curRoiDef=cellarrayRoi_defs{idxRoi};
      progress_bar(struct('progress', idxRoi-1, 'max', length(cellarrayRoi_defs), 'what', ['calculating roi ' num2str(curRoiDef.label) '...']));
      images=imphysfetch(ip.dfcalc([trials.index]));
      intensities=roimeasure(images, curRoiDef);
      clear images;
      if(needNewFigure)
         figRoiIntensities=figure('NumberTitle','off','Name','Roi Intensities');
         needNewFigure=0;
      end
      nRoisPerFig=nRoiPerRow*nRoiPerCol;
      position=mod(idxRoi, nRoisPerFig);
      if(position==0) position=nRoisPerFig; end
      figure(figRoiIntensities);
      tAxes=subplot(nRoiPerCol, nRoiPerRow, position);
      allRoiPlotAxes(idxRoi)=tAxes;
      plot(tAxes, [trials.valve], intensities, 'ro');
      set(tAxes, 'xtick', unique([trials.valve]));
      ylabel(tAxes, ip.dfcalc([trials(1).index(1)]).computation); % note: use first valve's first trial's computation name
      xlabel(tAxes, 'Valve');
      title(tAxes, ['ROI ' num2str(curRoiDef.label)]);
      setappdata(tAxes, 'roi_label', curRoiDef.label);
      setappdata(tAxes, 'figSummary', figSummary);
      setappdata(tAxes, 'figureRois', fig);
      install_mouse_event_handler(tAxes, 'up', @onClickOnAxes);
      if(position==nRoisPerFig)
         needNewFigure=1;
      end
      progress_bar(struct('progress', idxRoi, 'max', length(cellarrayRoi_defs), 'what', ['done with calculating roi ' num2str(curRoiDef.label) '...']));
   end % for, each roi
   ylims=get(allRoiPlotAxes, 'ylim');
   if(iscell(ylims))
      ylims=cell2mat(ylims);
   end
   ylim=[min(ylims(:,1)) max(ylims(:,2))];
   set(allRoiPlotAxes, 'ylim', ylim);
   

% --- Executes on button press in btnPlotRoisIntensities.
function btnPlotRoisIntensities_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlotRoisIntensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   plotRoiIntensities(handles.figureRois);
   

function result=onClickOnAxes(sender, event_args)
   if(is_button_down(sender, 'left')) % sender: the axes
      plot_valve_transition(sender);
   else
      % show_popup_menu(sender);
   end
   
   result=0;

   
% --- Executes on button press in btnResetTforms.
function btnResetTforms_Callback(hObject, eventdata, handles)
% hObject    handle to btnResetTforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fig=handles.figureRois;
   roi_defs=getappdata(fig, 'roi_defs');

   if(~isempty(roi_defs.label))
      uiwait(msgbox('There are ROIs defined. You must remove all ROIs before reset tforms','Manage ROIs','modal'));
      return;
   end

   % get the roi def from summary fig
   args=getappdata(fig, 'args');
   % index=getappdata(args.axes, 'index');
   figSummary=getappdata(get_parent_fig(args.axes), 'figSummary');
   rois=getappdata(figSummary, 'rois');
   if(~isempty(rois.defs_orig.label))
      uiwait(msgbox('There are ROIs defined. You must remove all ROIs before reset tforms','Manage ROIs','modal'));
      return;
   end
   
   for idx=1:length(rois.tform)
      rois.tform(idx)=unit_tform;
   end
   
   setappdata(figSummary, 'rois', rois);
   
   % now change the working copy also:
   roi_defs.tform=unit_tform;
   setappdata(fig, 'roi_defs', roi_defs);
   

