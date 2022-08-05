function varargout = implayer(varargin)
% IMPLAYER M-file for implayer.fig
%      IMPLAYER, by itself, creates a new IMPLAYER or raises the existing
%      singleton*.
%
%      H = IMPLAYER returns the handle to a new IMPLAYER or the handle to
%      the existing singleton*.
%
%      IMPLAYER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPLAYER.M with the given input arguments.
%
%      IMPLAYER('Property','Value',...) creates a new IMPLAYER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before implayer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to implayer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% 20041118: 
%   usage:
%      implayer();
%      implayer(struct('mayShink', 0)); % default mayShink=1
% 
% See also: GUIDE, GUIDATA, GUIHANDLES
   
% Edit the above text to modify the response to help implayer
   
% Last Modified by GUIDE v2.5 17-Mar-2005 14:10:43

% @history:
%    2/27/2004: 1. when user sets 0 as fps, play movie at full speed and
%                  update editFps with the measured speed.
%      2. change newly added controls' unit from "character" to pixel, and
%         adjust their position when an image is loaded.
%      3. 
%    -2/26/2004: 1. user can set movie play speed by fps. 
%      2. when user clicks slider bar during movie playback, stop playing;
%      3. add title/name to every displayed figure/window;
%      4. show the image filename on title bar of the main window;
%      5. ...

% @note: 
%    1. selectcircle() is used to handle mouse_move and mouse_up event,
%       of which mouse_move has to be the figure's event instead of axes's.
%       Matlab doesn't support control's mouse move event.

% Begin initialization code - DO NOT EDIT
   gui_Singleton = 1;
   gui_State = struct('gui_Name',       mfilename, ...
                      'gui_Singleton',  gui_Singleton, ...
                      'gui_OpeningFcn', @implayer_OpeningFcn, ...
                      'gui_OutputFcn',  @implayer_OutputFcn, ...
                      'gui_LayoutFcn',  [] , ...
                      'gui_Callback',   []);
   if nargin & isstr(varargin{1})
      gui_State.gui_Callback = str2func(varargin{1});
   end
   
   if nargout
      [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
   else
      gui_mainfcn(gui_State, varargin{:});
   end
% End initialization code - DO NOT EDIT


% --- Executes just before implayer is made visible.
function implayer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to implayer (see VARARGIN)
   
% Load in a new image stack
%[handles.imstruct,handles.image] = imanalyzeload(handles);
   handles.ip=varargin{1};
   handles.ip_raw=varargin{2};
   if(nargin == 6)
      handles.imanalyzeParams=varargin{3};
   else
      handles.imanalyzeParams=[];
   end
   if(~isfield(handles.imanalyzeParams, 'mayShink'))
      handles.imanalyzeParams.mayShink=1;
   end
   handles.image = [];
   handles.imstruct = struct;
   handles.isPlaying = 0;
   handles.contrastindex = [];
   handles.rawcontrastrange = [];
   handles.subcontrastrange = [];
   handles.colorsaturated = 0;
   handles.saturated = [0 4095];
   handles.regions = zeros(3,0);  % Rows are x,y,r
   handles.regionshcircles = [];
   handles.currentframe = 1;
   handles.isPlayForward=1; % 

   % Choose default command line output for implayer
   handles.output = hObject;
   
   % Update handles structure
   guidata(hObject, handles);
   
   loadbutton_Callback(hObject, eventdata, handles);

   
% UIWAIT makes implayer wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = implayer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% Get default command line output from handles structure
   varargout{1} = handles.output;
   
% Slider
% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
   usewhitebg = 1;
   if usewhitebg
      set(hObject,'BackgroundColor',[.9 .9 .9]);
   else
      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   end
   
% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
   set(handles.playbutton, 'value', 0);
   handles.currentframe = round(get(hObject,'Value'));
   guidata(hObject,handles);
   implayer_showframe(handles);
   implayer_guistate(handles);
   
% --- Executes on button press in playbutton.
function playbutton_Callback(hObject, eventdata, handles)
% Determine whether we are playing or stopping
   handles.isPlaying = get(hObject,'Value');
   if ~handles.isPlaying
      set(handles.playbutton,'String','Play');
      implayer_guistate(handles);
   else
      set(handles.playbutton,'String','Stop');
      implayer_guistate(handles);
      handles.currentframe = implayer_playmovie(handles);
      % We may have reached the end of the movie without being interrupted
      % if (length(handles.imstruct) == handles.currentframe)
         handles.isPlaying = 0;
         set(handles.playbutton,'Value',0,'String','Play')
      % end
      guidata(hObject,handles);
      implayer_guistate(handles);
   end
   
% Framenumber
% --- Executes during object creation, after setting all properties.
function framenumber_CreateFcn(hObject, eventdata, handles)
   if ispc
      set(hObject,'BackgroundColor','white');
   else
      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   end
   
function framenumber_Callback(hObject, eventdata, handles)
   currentframe = str2num(get(hObject,'String'));
   if ~isempty(currentframe)
      handles.currentframe = currentframe;
      guidata(hObject,handles);
      implayer_showframe(handles);
      implayer_guistate(handles);
   end
   
% --- Executes on button press in adjustcontrast.
function adjustcontrast_Callback(hObject, eventdata, handles, is_ask_user)
% First, obtain the values in the image (either raw or
% background-subtracted)
   if(nargin==3) is_ask_user=1; end
   
   imgtmp = double(handles.imstruct(handles.currentframe).data);
   rng = handles.rawcontrastrange;
   subtracting = 0;
   if(is_ask_user)
      % Now call the range-setting GUI
      [tmprange,tmpcs] = imrangegui(imgtmp,rng,handles.colorsaturated);
   else % use the param from struct handles
      tmprange=handles.rawcontrastrange;
      tmpcs =handles.colorsaturated;
   end
   if isempty(tmprange)
      return; % User hit cancel, skip the rest
   end
   handles.colorsaturated = tmpcs;
   if subtracting
      handles.subcontrastrange = tmprange;
   else
      tmprange = round(tmprange);
      handles.rawcontrastrange = tmprange;
      % handles.contrastindex = imrange_to_colormap(tmprange,handles.saturated(2)+1);
      if (handles.colorsaturated)
         % Color saturated pixels blue and red
         handles.contrastindex(1,:) = [0 0 1]; % Blue
         handles.contrastindex(end,:) = [1 0 0]; % Red
         set(handles.figure1, 'colormap', [ [0 0 1]; gray(254); [1 0 0] ]);
      else
         set(handles.figure1, 'colormap', [  gray(256) ]);
      end
      % colormap(handles.contrastindex);
      set(gca, 'clim', tmprange);
      
      % jason: a hack to solve matlab colormap problem on windows 20040922
      if(isa(handles.imstruct(1).data(1), 'uint8'))
         colormap(handles.contrastindex(1:256,:));   
      end
      
      % jason: another hack to solve matlab color map problem on windows
      if(~isunix && isa(handles.imstruct(1).data(1), 'uint16'))
         if(handles.colorsaturated)
            tColorBegin=[0 0 1]; % Blue
            tColorEnd  =[1 0 0]; % Red;
         else
            tColorBegin=[0 0 0]; % black
            tColorEnd  =[1 1 1]; % white
         end
         
         % for direct mapping, a bug in matlab on windows prevent following
         % code from working:
         % tColors=[repmat(tColorBegin, tmprange(1), 1); gray(tmprange(2)-tmprange(1)); repmat(tColorEnd, 1024-tmprange(2), 1) ];
         % colormap(tColors);  
         
         % for scaled mapping:
         set(handles.image, 'CDataMapping', 'scaled');
         tColors=[tColorBegin; gray(254); tColorEnd];
         colormap(tColors);
         set(gca, 'clim', tmprange); % gca = handles.imaxes
      end

   end
   guidata(hObject,handles);
   % Put GUI in consistent state
   implayer_guistate(handles);
   implayer_showframe(handles);
  
   
% --- Executes on button press in showregions.
function showregions_Callback(hObject, eventdata, handles)
   if get(hObject,'Value')
      state = 'on';
   else
      state = 'off';
   end
   set(handles.regionshcircles,'Visible',state);
   % Because erase mode is set to 'none', we have to redraw the screen if the
   % user is turning the regions off
   if ~get(hObject,'Value')
      set(handles.image,'CData',[]); % This seems needed to force a redraw
      implayer_showframe(handles);
   end
   
function result=weighted_mean_old1(preFrame, curFrame)
% preFrame: the region's raw data in previous frame
% curFrame: the region's raw data in current frame
% result: the weighted mean value of the region of current frame
   % result=mean(curFrame);
   curFrame=double(curFrame);
   preFrame=double(preFrame);
   % prefer "increasing edge", may use exp later:
   tWeight=1+2*(curFrame-preFrame)./(curFrame+preFrame);
   tIndices=find(tWeight<1);
   tWeight(tIndices)=1;
   tIndices=find(isinf(tWeight));
   tWeight(tIndices)=1;
   result=mean(curFrame.*tWeight);
   
function result=weighted_mean(preFrame, curFrame, min_mean, max_mean)
% preFrame: the region's raw data in previous frame
% curFrame: the region's raw data in current frame
% result: the weighted mean value of the region of current frame
   result=mean(curFrame);
   result=result*exp(3*min_mean/(max_mean-min_mean)*(result-min_mean)/min_mean);

% --- Executes on button press in showintensities.
function showintensities_Callback(hObject, eventdata, handles)
   ttHandles=handles;
   ttHandles.imstruct=ttHandles.imstruct_raw;
   handles.intensity=showintensities(ttHandles, 'ROI Intensity - raw'); % this is a temp workaround. should modify showintensities

   guidata(hObject,handles);

   
function intensity=showintensities(handles, caption)
% Calculate all the intensities
%
% First, get the list of pixels inside each circle
   nregions = size(handles.regions,2);
   nframes = length(handles.imstruct);
   for i = 1:nregions
      % Form the square region, then reject points
      r = handles.regions(3,i);
      xrange = round(handles.regions(1,i)) + (-ceil(r):ceil(r));
      yrange = round(handles.regions(2,i)) + (-ceil(r):ceil(r));
      xlist = repmat(xrange,1,length(xrange));
      ylist = repmat(yrange,length(yrange),1);
      xylist = [xlist(:) ylist(:)]; % Here's the set of points within square of width r
      dxylist = xylist - repmat(handles.regions(1:2,i)',size(xylist,1),1);
      distxylist = sqrt(sum(dxylist.^2,2));
      indxkeep = find(distxylist < r);
      % Convert x,y coords to an index
      sz = size(get(handles.image,'CData'));
      imindx{i} = sub2ind(sz,xylist(indxkeep,2),xylist(indxkeep,1));
   end

   % Loop over frames and find the min of mean
   mag = zeros(nframes,nregions);
   for i = 1:nframes
      img = handles.imstruct(i).data;
      for j = 1:nregions
         mag(i,j) = mean(img(imindx{j}));
      end
   end
   min_mean=min(mag, [], 1);
   max_mean=max(mag, [], 1);
   
   % Loop over frames and compute mean intensities
   mag = zeros(nframes,nregions);
   for i = 1:nframes
      img = handles.imstruct(i).data;
      preImg=handles.imstruct(max(i-1,1)).data;
      for j = 1:nregions
         if(get(handles.cbWeightedMean, 'value'))
            mag(i,j) = weighted_mean(preImg(imindx{j}), img(imindx{j}), min_mean(j), max_mean(j));
         else
            mag(i,j) = mean(img(imindx{j}));
         end
      end
   end
   % Plot the data
   figure
   tFilterScope=str2num(get(handles.editFilterScope, 'string'));
   if(tFilterScope>0)
      if(tFilterScope<1)
         tN=ceil(size(mag, 1)*tFilterScope);
      else
         tN=ceil(tFilterScope);
      end
      mag=medfilt1(mag, tN);
   end
   hlines = plot(mag);
   % Make sure colors/markers are the same and install mouse-click event to
   % toggle markers
   for i = 1:nregions
      set(hlines(i),'Color',get(handles.regionshcircles(i),'Color'));
      set(hlines(i), 'marker', get(handles.regionshcircles(i), 'marker') );
      set(hlines(i), 'userdata', handles.regionshcircles(i) );
      set(hlines(i), 'buttonDownFcn', @intensity_line_mousedown_callback);
   end
   xlabel('Frame');
   ylabel('Mean Pixel Intensity');
   %title('ROI Intensity');
   set(gcf, 'name', caption);
   set(gca, 'xlim', [1 nframes]); % ignore the last frame
   if(~isempty(handles.stimSeqInFrame))
      % set(gca, 'xtick', handles.stimSeqInFrame(:,2)-handles.nFramesToDiscard);
      % Xs=handles.stimSeqInFrame(:,2)-handles.nFramesToDiscard;
      
      % set(gca, 'xgrid', 'on');
      % set(gca, 'xminortick', 'on');
      Ys=get(gca, 'ylim');
      yOffset=(Ys(2)-Ys(1))/16;
      ttY=mean(Ys)-yOffset;
      for idx=1:size(handles.ip,2)
	 ttX=idx;
         h=line([ttX, ttX], [Ys(1), Ys(2)]);
         if(handles.ip(idx).stimulus==0)
            set(h, 'color', 'blue');
         else
            set(h, 'color', 'red');
         end
         set(h, 'linestyle', ':');
	 text(ttX, ttY+mod(idx,2)*2*yOffset, num2str(handles.ip(idx).stimulus));
      end
   end
   
   % set(gca, 'xtickLabel', [handles.ip.stacknum]);
   set(gca, 'xtickLabel', [handles.ip(get(gca, 'xtick')).stacknum]);
   
   % Save the data
   intensity = mag;

   
function intensity_line_mousedown_callback(sender, eventdata)
   if(get(sender, 'marker')=='*')
      set(sender, 'marker', 'none');
   else
      set(sender, 'marker', '*');
   end
   h=get(sender, 'userdata');
   if(ishandle(h))
      set(h, 'marker', get(sender, 'marker') );
      refresh_image(h, []);
   end
   
% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
   if(~isfield(handles, 'intensity'))
      errordlg(['click button "show intensities" to calculate and' 10 ...
		'verify the intensities then save them again'], ...
	       'No intensities exist');
      return;
   end
   
   [pathstr,name,ext] = fileparts(handles.filename);
   [thefilename,pathname] = uiputfile([pathstr name '.mat'],'Save data to file...')
   tRegions = handles.regions;
   tRegions = relRegion2absRegion(tRegions, handles.ip);
   tIntensity = handles.intensity;
   tFilename = handles.filename;
   save([pathname thefilename],'tFilename','tRegions','tIntensity');
% Still to be written


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
   handles=guidata(hObject);
   [handles.imstruct,handles.image,handles.filename, stimSeqInFrame, handles.imstruct_raw] = imanalyzeload(handles);
   handles.stimSeqInFrame=stimSeqInFrame;
   % Now set a reasonable default contrasts--stretch the min/max of the image
   % into the available range of greyscales
   immin = min(min(handles.imstruct(1).data));
   immax = max(max(handles.imstruct(1).data));
   handles.rawcontrastrange = double([immin immax]);
   % handles.contrastindex = imrange_to_colormap([immin immax],handles.saturated(2)+1);
   if (handles.colorsaturated)
      % Color saturated pixels blue and red
      handles.contrastindex(1,:) = [0 0 1]; % Blue
      handles.contrastindex(end,:) = [1 0 0]; % Red
      set(handles.figure1, 'colormap', [ [0 0 1]; gray(254); [1 0 0] ]);
   else
      set(handles.figure1, 'colormap', [  gray(256) ]);
   end
   set(0, 'CurrentFigure', handles.figure1);
   axes(handles.imaxes);
   % colormap(handles.contrastindex);
   % set(handles.figure1, 'Colormap', handles.contrastindex);
   set(handles.imaxes, 'clim', handles.rawcontrastrange);
   
   % jason: a hack to solve matlab colormap problem on windows 20040922
   if(isa(handles.imstruct(1).data(1), 'uint8'))
      % colormap(handles.contrastindex(1:256,:)); 
      set(handles.figure1, 'Colormap', handles.contrastindex(1:256,:) );
   end
   
   % jason: another hack to solve matlab color map problem on windows
   if(~isunix && isa(handles.imstruct(1).data(1), 'uint16'))
         if(handles.colorsaturated)
            tColorBegin=[0 0 1]; % Blue
            tColorEnd  =[1 0 0]; % Red;
         else
            tColorBegin=[0 0 0]; % black
            tColorEnd  =[1 1 1]; % white
         end
         
         % for direct mapping, a bug in matlab on windows prevent following
         % code from working:
         % tColors=[repmat(tColorBegin, tmprange(1), 1); gray(tmprange(2)-tmprange(1)); repmat(tColorEnd, 1024-tmprange(2), 1) ];
         % colormap(tColors);  
         
         % for scaled mapping:
         set(handles.image, 'CDataMapping', 'scaled');
         tColors=[tColorBegin; gray(254); tColorEnd];
         % colormap(tColors);
         set(handles.figure1, 'Colormap', tColors);
         set(gca, 'clim', handles.rawcontrastrange); % gca = handles.imaxes
   end

   
   % clear old region info:
   % delete(handles.regionshcircles);
   % handles.regionshcircles=[];
   % handles.regions=zeros(3,0);  % Rows are x,y,r
   
   guidata(hObject,handles);
   % Put GUI in consistent state
   implayer_guistate(handles);
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refresh_image(hObject, eventdata)
   handles = guidata(hObject);
   set(handles.image,'CData',[]); % This seems needed to force a redraw
   implayer_showframe(handles);

% Function to handle region selection
function regions_callback(hObject, eventdata)
   if(get(hObject, 'marker')=='*')
      set(hObject, 'marker', 'none');
   else
      set(hObject, 'marker', '*');
   end
   refresh_image(hObject, eventdata);
   
function implayer_selectstart(hObject, eventdata)
   handles = guidata(hObject);
   if ~get(handles.showregions,'Value')
      set(handles.showregions,'Value',1); % Show all existing regions
      set(handles.regionshcircles,'Visible','on');
   end
   % Figure out the appropriate color for the new circle and create it
   co = get(handles.imaxes,'ColorOrder');
   colindex = size(handles.regions,2);
   col = co(mod(colindex,size(co,1))+1,:);
   options.color = col;
   
   % Create the selection circle
   [cp,hcircle] = selectcircle(hObject, eventdata, options);
   tMinRegionRadius=str2num(get(handles.editMinRegionRadius, 'string'));
   if(cp(3)<tMinRegionRadius)
      delete(hcircle);
      return;
   end
   
   % Add this new information
   handles.regions = [handles.regions, cp'];
   handles.regionshcircles(end+1) = hcircle;
   set(hcircle,'EraseMode','none');
   set(hcircle, 'ButtonDownFcn', @regions_callback);
   guidata(hObject,handles);
   implayer_guistate(handles);
   implayer_showframe(handles);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_ip(hObject,eventData,handles, varargin)
   set(0, 'currentFigure', hObject);
   newip=varargin{1};
   newip_raw=varargin{2};
   gd=guidata(hObject);
   gd.ip=newip;
   gd.ip_raw=newip_raw;
   guidata(hObject, gd);
   loadbutton_Callback(hObject, eventData, gd);

% A function to handle loading of a new image file, and updating the GUI
% accordingly.
function [imstruct,himg,imfilename, stimSeqInFrame, imstruct_raw] = imanalyzeload(handles)
   ip=handles.ip;
   ip_raw=handles.ip_raw;
   imfilename = ip(1).imfile; % todo:
   [pathstr,basename,extname] = fileparts(imfilename);

   nimages = length(ip);
   for i = 1:nimages
      imstruct(i).data = imphysfetch(ip(i));
      imstruct(i).width =  size(imstruct(i).data,2);
      imstruct(i).height = size(imstruct(i).data,1);
      imstruct_raw(i).data = imphysfetch(ip_raw(i));
      imstruct_raw(i).width =  size(imstruct_raw(i).data,2);
      imstruct_raw(i).height = size(imstruct_raw(i).data,1);
      if(mod(i,5)==1 || i==nimages)
         progress_bar(struct('progress', i, 'max', nimages, 'what', ['loading ' [basename extname] ' ...']));
      end
   end

   descFileName=ip(1).headerfile; % todo:
   if(isempty(descFileName))
      stimSeqInFrame=[];
   else
      % open the stim file
      [s, txt]=load_text_file(descFileName);
      if(s~=0)
	 error(['error to open file: ' descFileName]);
      else
	 strStim=key2value(txt, 'stimulus sequence, time in frames');
	 stimSeqInFrame=eval(['[' strStim ']' ]);
      end
   end
   
   % Resize the figure to adjust to the size of the image
   axsize = [imstruct(1).width imstruct(1).height];
   tPosPanelIntensity=get(handles.uipanelIntensity, 'position');

   if(handles.imanalyzeParams.mayShink)
      tMaxVisibleAxSize=get(0, 'screensize');
      tMaxVisibleAxSize=tMaxVisibleAxSize(3:4)-tPosPanelIntensity([3 2]) -[100 100] ;
      tRatioHtoW=imstruct(1).height/imstruct(1).width;
      if(tMaxVisibleAxSize(1)*tRatioHtoW<=tMaxVisibleAxSize(2))
	 tMaxVisibleAxSize(2)=tMaxVisibleAxSize(1)*tRatioHtoW;
      else
	 tMaxVisibleAxSize(1)=tMaxVisibleAxSize(2)/tRatioHtoW;
      end
      axsize=min([axsize; tMaxVisibleAxSize]);
   end
   minFigSize=[745 535];
   figsize = axsize + tPosPanelIntensity([3 2])+[20,10];
   figsize = max([minFigSize; figsize]);
   set(0, 'CurrentFigure', handles.figure1);
   figpos = get(handles.figure1,'Position');
   set(handles.figure1,'Position',[12 75 figsize]);
   oldPos=get(handles.uipanelMovie, 'position');
   set(handles.uipanelMovie, 'position', [oldPos(1:2) figsize(1)-20 oldPos(4)]);
   oldPos=get(handles.slider, 'position');
   set(handles.slider, 'position', [oldPos(1:2) figsize(1)-25-oldPos(1) oldPos(4)]);
   % Adjust the position of all the controls
   % Do the controls on the right column
   tOldImgPos=get(handles.imaxes, 'position');
   set(handles.imaxes,'Position',[tOldImgPos(1:2) axsize]);
       
   % Now set up the image data
   axes(handles.imaxes);
   if(isempty(handles.image))
      himg = image(imstruct(1).data, 'CDataMapping','scaled');
      set(handles.figure1, 'colormap', [ [0 0 1]; gray(254); [1 0 0] ]);
   else
      set(handles.image, 'CData',imstruct(1).data);
      himg=handles.image;
   end
   set(handles.imaxes,'Visible','off');
   % hook the mouse down event of the image:
   set(himg,'EraseMode','none','ButtonDownFcn',@implayer_selectstart);
   % Update the GUI data
   set(handles.frame_of_text,'String',['  of ' num2str(nimages)]);
   set(handles.slider,'Min',1,'Max',nimages,'Value',1,'SliderStep',[1/(nimages-1) 0.1]);
   set(handles.framenumber,'String','1');
   set(handles.figure1, 'name', [basename ' -- normal']);
   
   
function implayer_guistate(handles)
   allhandles = [handles.slider handles.playbutton  handles.framenumber ...
                 handles.frame_of_text  ...
                 handles.adjustcontrast handles.showregions ...
                 handles.showintensities handles.savebutton ];
   if isempty(handles.image)
      % If there is no file loaded, turn off everything but the "Load..."
      % button
      hon = handles.loadbutton;
      hoff = setdiff(allhandles,hon);
      set(hoff,'Enable','off');
      set(hon,'Enable','on');
   else
      if get(handles.playbutton,'Value')
         % Are we playing? Turn off all but:
         %   Play/Stop button
         %   Frame number text
         %   Subtract background checkbox
         hon = [handles.playbutton  handles.framenumber handles.frame_of_text handles.slider];
         hoff = setdiff(allhandles,hon);
         set(hoff,'Enable','off');
         set(hon,'Enable','on');
         % Also set the playbutton text to "Stop"
         set(handles.playbutton,'String','Stop');
      else
         % We're not playing, but we do have a file loaded. Set up states on
         % buttons appropriately.
         set(allhandles,'Enable','on');
         if isempty(handles.regions)
            set([handles.showregions handles.showintensities],'Enable','off');
            set(handles.showregions,'Value',0);
         end
         if ( isempty(handles.regions))
            set(handles.savebutton,'Enable','off');
         end
         % If we're at the last frame, disable the play button
         if (handles.currentframe == length(handles.imstruct))
            set(handles.playbutton,'Enable','off');
         end
      end
   end
   
 
function result = getNextFrame(handles)
   if(get(handles.cbReverse, 'value'))
      set(handles.cbReverse, 'value', 0);
      if(handles.isPlayForward)
         handles.isPlayForward=0;
      else
         handles.isPlayForward=1;
      end
   end
   
   if(handles.isPlayForward)
      next=handles.currentframe+1;
   else
      next=handles.currentframe-1;
   end
   if(get(handles.cbContPlay, 'value'))
      if(next>length(handles.imstruct))
         if(get(handles.cbBackforth, 'value'))
            next =length(handles.imstruct)-1;
            handles.isPlayForward=0;
         else
            next=1;
         end
      end
      if(next<1)
         if(get(handles.cbBackforth, 'value'))
            next=2;
            handles.isPlayForward=1;
         else
            next=length(handles.imstruct);
         end
      end
   end
   handles.currentframe=next;
   result=handles;
   
function currentframe = implayer_playmovie(handles)
   tic;
   fps=str2num(get(handles.editFps, 'string'));
   nFramesPlayed =0;
   while (handles.currentframe>=1 && handles.currentframe <= length(handles.imstruct) ...
          && get(handles.playbutton,'Value'))
      implayer_showframe(handles);
      handles = getNextFrame(handles);
      nFramesPlayed=nFramesPlayed+1;
      if(fps>0)
         curTime=toc;
         supposedTime=nFramesPlayed/fps;
         if(supposedTime>curTime)
            sleep_g(supposedTime-curTime);    
         end
      end
   end
   if(fps==0)
      fps=nFramesPlayed/toc;
      set(handles.editFps, 'string', num2str(round(fps)));
   end
   
   currentframe=handles.currentframe;
   if(handles.currentframe>length(handles.imstruct))
      currentframe = length(handles.imstruct);
   elseif(handles.currentframe<1)
      currentframe = 1;
   end
      
   
function implayer_showframe(handles)
% Drawing is done with xor/none erasing, so first we have to turn off the
% regions if they are showing
   if get(handles.showregions,'Value')
      set(handles.regionshcircles,'Visible','off');
   end
   img = handles.imstruct(handles.currentframe).data;
   % Render the image
   set(handles.image,'CData',img);
   % Turn the regions back on if needed
   if get(handles.showregions,'Value')
      set(handles.regionshcircles,'Visible','on');
   end
   % Update the slider and edit boxes with the current frame number
   set(handles.framenumber,'String',num2str(handles.currentframe));
   set(handles.slider,'Value',handles.currentframe);
   drawnow
   
   
% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over slider.
function slider_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   set(handles.playbutton, 'value', 0);
   
% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
   
   
% --- Executes during object creation, after setting all properties.
function editFps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
   
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
   if ispc
      set(hObject,'BackgroundColor','white');
   else
      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   end
   
   
   
function editFps_Callback(hObject, eventdata, handles)
% hObject    handle to editFps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% Hints: get(hObject,'String') returns contents of editFps as text
%        str2double(get(hObject,'String')) returns contents of editFps as a double
   
   
% --- Executes on button press in btnIncFps.
function btnIncFps_Callback(hObject, eventdata, handles)
% hObject    handle to btnIncFps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fps=str2num(get(handles.editFps, 'string'));
   set(handles.editFps, 'string', num2str(fps+1));
   
% --- Executes on button press in btnDecFps.
function btnDecFps_Callback(hObject, eventdata, handles)
% hObject    handle to btnDecFps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   fps=str2num(get(handles.editFps, 'string'));
   if(fps<2)
      fps=2;    
   end
   set(handles.editFps, 'string', num2str(fps-1));
   

% --- Executes on mouse press over axes background.
function imaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to imaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cbContPlay.
function cbContPlay_Callback(hObject, eventdata, handles)
% hObject    handle to cbContPlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbContPlay


% --- Executes on button press in cbBackforth.
function cbBackforth_Callback(hObject, eventdata, handles)
% hObject    handle to cbBackforth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbBackforth


% --- Executes on button press in cbReverse.
function cbReverse_Callback(hObject, eventdata, handles)
% hObject    handle to cbReverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbReverse
% $$$    if(get(hObject, 'value'))
% $$$       if(handles.isPlayForward)
% $$$          handles.isPlayForward=0;
% $$$       else
% $$$          handles.isPlayForward=1;
% $$$       end
% $$$       guidata(hObject, handles);
% $$$    end
   


% --- Executes during object creation, after setting all properties.
function imaxes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate imaxes


% --- Executes on button press in btnSelectLine.
function btnSelectLine_Callback(hObject, eventdata, handles)
% hObject    handle to btnSelectLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function editCurLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCurLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editCurLine_Callback(hObject, eventdata, handles)
% hObject    handle to editCurLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCurLine as text
%        str2double(get(hObject,'String')) returns contents of editCurLine as a double


% --- Executes on button press in btnDeleteRegions.
function btnDeleteRegions_Callback(hObject, eventdata, handles)
% hObject    handle to btnDeleteRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   userResponse=questdlg('Delete selected or all regions?', 'implayer', ...
                         'All', 'Selected', 'Selected');
   if(strcmp(userResponse, 'All'))
      tIndices=1:length(handles.regionshcircles);
   else
      tIndices=[];
      for i=1:length(handles.regionshcircles)
         tH=handles.regionshcircles(i);
         if(get(tH, 'marker')=='*')
            tIndices=[tIndices i];
         end
      end
   end
   
   delete(handles.regionshcircles(tIndices));
   handles.regionshcircles(tIndices)=[];
   handles.regions(:, tIndices)=[];
   guidata(hObject,handles);
   implayer_guistate(handles);
   refresh_image(hObject, eventdata);


% --- Executes during object creation, after setting all properties.
function editMinRegionRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinRegionRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editMinRegionRadius_Callback(hObject, eventdata, handles)
% hObject    handle to editMinRegionRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinRegionRadius as text
%        str2double(get(hObject,'String')) returns contents of editMinRegionRadius as a double

function result=relRegion2absRegion(aRegions, ip)
   result=aRegions+repmat([ip(1).xrange(1); ip(1).yrange(1);0], ...
			  1, size(aRegions, 2) );

% --- Executes on button press in btnSaveRegions.
function btnSaveRegions_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   [pathstr,name,ext] = fileparts(handles.filename);
   [thefilename,pathname] = uiputfile([pathstr name '.region'],'Save regions to file...');
   if(thefilename==0) 
      return; % no file selected
   end
   tRegions=handles.regions;
   tRegions=relRegion2absRegion(tRegions, handles.ip);
   
   tFilename=handles.filename;
   save([pathname thefilename], 'tRegions', 'tFilename', '-MAT');

   
function result=absRegion2relRegion(aRegions, ip)
   result=aRegions-repmat([ip(1).xrange(1); ip(1).yrange(1);0], ...
			  1, size(aRegions, 2) );
   
% --- Executes on button press in btnLoadRegions.
function btnLoadRegions_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   [imfilename,impathname] = uigetfile({'*.region'},'Pick a region definition file');
   if(imfilename==0)
      return; % no file is selected
   end
   load([impathname imfilename], '-MAT');
   % now we get tRegions and tFilename (unused):
   tRegions=absRegion2relRegion(tRegions, handles.ip);
   
   axes(handles.imaxes);
   for i=1:size(tRegions, 2)
      % Figure out the appropriate color for the new circle and create it
      co = get(handles.imaxes,'ColorOrder');
      colindex = size(handles.regions,2)+i-1;
      col = co(mod(colindex,size(co,1))+1,:);
      options.color = col;
      
      [x,y]=get_circle_points(tRegions(1,i), tRegions(2,i), tRegions(3,i));
      hcircle = line(x,y,'Color',options.color,'EraseMode','none');
      set(hcircle,'UserData',struct('center',tRegions(1:2,i),'radius',tRegions(3,i)));
      handles.regionshcircles(end+1) = hcircle;
      set(hcircle, 'ButtonDownFcn', @regions_callback);
   end
   handles.regions=[handles.regions tRegions];
   guidata(hObject,handles);
   implayer_guistate(handles);
   refresh_image(hObject, eventdata);
   
% cp from selectcircle/selectcircle_draw()
function [x,y] = get_circle_points(ax, ay,radius)
   npts = 100;
   th = linspace(0,2*pi,npts);
   x = ax+radius*cos(th);
   y = ay+radius*sin(th);


% --- Executes on button press in cbWeightedMean.
function cbWeightedMean_Callback(hObject, eventdata, handles)
% hObject    handle to cbWeightedMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbWeightedMean


% --- Executes on button press in cbFilter.
function cbFilter_Callback(hObject, eventdata, handles)
% hObject    handle to cbFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbFilter


% --- Executes during object creation, after setting all properties.
function editFilterScope_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFilterScope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function editFilterScope_Callback(hObject, eventdata, handles)
% hObject    handle to editFilterScope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFilterScope as text
%        str2double(get(hObject,'String')) returns contents of editFilterScope as a double


% --- Executes on button press in btnNextRound.
function btnNextRound_Callback(hObject, eventdata, handles)
% hObject    handle to btnNextRound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   set(handles.comboxPlayWhat, 'value', 2); % raw
   set(handles.figure1, 'waitstatus', 0); % 0: not wait. used by imstimplay()


% --- Executes on selection change in comboxPlayWhat.
function comboxPlayWhat_Callback(hObject, eventdata, handles)
% hObject    handle to comboxPlayWhat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns comboxPlayWhat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from comboxPlayWhat
   contents = get(hObject,'String');
   selection=contents{get(hObject,'Value')};
   if(strcmp(selection,'filtered subtract'))
      filter_hsize=getappdata(handles.figure1, 'filter_hsize');
      filter_sigma=getappdata(handles.figure1, 'filter_sigma');
      smoothfilt = fspecial('gaussian',filter_hsize, filter_sigma); 
   end
   
   switch selection
      case 'raw'
         handles.imstruct=handles.imstruct_raw;
      case {'unfiltered subtract', 'filtered subtract'}
         ip=handles.ip_raw;
         imfilename = ip(1).imfile; % todo:
         [pathstr,basename,extname] = fileparts(imfilename);

         nimages = length(ip);
         for i = 1:nimages
            imstruct(i).data = imphysfetch(ip(i));
            if(strcmp(selection,'filtered subtract'))
               imstruct(i).data = imfilter(single(imstruct(i).data) - ip(i).background,smoothfilt);
            else
               imstruct(i).data = feval(class(handles.imstruct(i).data), ...
                                        single(imstruct(i).data) - ip(i).background);
            end
            imstruct(i).width =  size(imstruct(i).data,2);
            imstruct(i).height = size(imstruct(i).data,1);
            if(mod(i,5)==1 || i==nimages)
               progress_bar(struct('progress', i, 'max', nimages, ...
                                   'what', ['processing ' [basename extname] ' ...']));
            end
         end
         handles.imstruct=imstruct;
   end % switch
   guidata(hObject, handles);
   
   % adjust contrast
   immin = min(min(handles.imstruct(1).data));
   immax = max(max(handles.imstruct(1).data));
   handles.rawcontrastrange = double([immin immax]);
   adjustcontrast_Callback(hObject, [], handles, 0);
   
   implayer_showframe(handles);
      
      


% --- Executes on button press in btnChgFilter.
function btnChgFilter_Callback(hObject, eventdata, handles)
% hObject    handle to btnChgFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   dlgTitle='change filter parameters';
   prompt={'size:', ...
           'sigma:', ...
          };
   filter_hsize=getappdata(handles.figure1, 'filter_hsize');
   filter_sigma=getappdata(handles.figure1, 'filter_sigma');
   defaultValues={num2str(filter_hsize), num2str(filter_sigma)};          
   lineNo=1;
   answer=inputdlg(prompt,dlgTitle,lineNo,defaultValues);
   if(~isempty(answer))
      filter_hsize=str2num(answer{1});
      filter_sigma=str2num(answer{2});
      setappdata(handles.figure1, 'filter_hsize', filter_hsize);
      setappdata(handles.figure1, 'filter_sigma', filter_sigma);
      contents = get(handles.comboxPlayWhat,'String');
      selection=contents{get(handles.comboxPlayWhat,'Value')};
      if(strcmp(selection,'filtered subtract'))
         comboxPlayWhat_Callback(handles.comboxPlayWhat, [], handles);
      end
   end

