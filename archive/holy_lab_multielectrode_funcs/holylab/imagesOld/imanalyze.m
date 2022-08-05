function varargout = imanalyze(varargin)
% IMANALYZE M-file for imanalyze.fig
%      IMANALYZE, by itself, creates a new IMANALYZE or raises the existing
%      singleton*.
%
%      H = IMANALYZE returns the handle to a new IMANALYZE or the handle to
%      the existing singleton*.
%
%      IMANALYZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMANALYZE.M with the given input arguments.
%
%      IMANALYZE('Property','Value',...) creates a new IMANALYZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imanalyze_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imanalyze_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% 20041118: 
%   usage:
%      imanalyze();
%      imanalyze(struct('mayShink', 0)); % default mayShink=1
% 
% See also: GUIDE, GUIDATA, GUIHANDLES
   
% Edit the above text to modify the response to help imanalyze
   
% Last Modified by GUIDE v2.5 05-Mar-2004 12:49:17

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
                      'gui_OpeningFcn', @imanalyze_OpeningFcn, ...
                      'gui_OutputFcn',  @imanalyze_OutputFcn, ...
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


% --- Executes just before imanalyze is made visible.
function imanalyze_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imanalyze (see VARARGIN)
   
% Load in a new image stack
%[handles.imstruct,handles.image] = imanalyzeload(handles);
   if(nargin == 4)
      handles.imanalyzeParams=varargin{1};
   else
      handles.imanalyzeParams=[];
   end
   if(~isfield(handles.imanalyzeParams, 'mayShink'))
      handles.imanalyzeParams.mayShink=1;
   end
   handles.image = [];
   handles.imstruct = struct;
   handles.isPlaying = 0;
   handles.background = [];
   handles.contrastindex = [];
   handles.rawcontrastrange = [];
   handles.subcontrastrange = [];
   handles.colorsaturated = 1;
   handles.saturated = [0 4095];
   handles.regions = zeros(3,0);  % Rows are x,y,r
   handles.regionshcircles = [];
   handles.currentframe = 1;
   handles.isPlayForward=1; % 

   tCurUser=getenv('USER');
   if(strcmp(tCurUser, 'holekamt'))
       cd('/usr/lab/home/holekamt');    
   end
   
   % Choose default command line output for imanalyze
   handles.output = hObject;
   
   % Update handles structure
   guidata(hObject, handles);
   
% UIWAIT makes imanalyze wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imanalyze_OutputFcn(hObject, eventdata, handles)
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
   imanalyze_showframe(handles);
   imanalyze_guistate(handles);
   
% --- Executes on button press in playbutton.
function playbutton_Callback(hObject, eventdata, handles)
% Determine whether we are playing or stopping
   handles.isPlaying = get(hObject,'Value');
   if ~handles.isPlaying
      set(handles.playbutton,'String','Play');
      imanalyze_guistate(handles);
   else
      set(handles.playbutton,'String','Stop');
      imanalyze_guistate(handles);
      handles.currentframe = imanalyze_playmovie(handles);
      % We may have reached the end of the movie without being interrupted
      % if (length(handles.imstruct) == handles.currentframe)
         handles.isPlaying = 0;
         set(handles.playbutton,'Value',0,'String','Play')
      % end
      guidata(hObject,handles);
      imanalyze_guistate(handles);
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
      imanalyze_showframe(handles);
      imanalyze_guistate(handles);
   end
   
   
% --- Executes on button press in computebackground.
function computebackground_Callback(hObject, eventdata, handles)
   framerange = str2num([get(handles.backframemin,'String') ' ' get(handles.backframemax,'String')]);
   fs = zeros(size(handles.imstruct(framerange(1))));
   for i = framerange(1):framerange(2)
      fs = fs + double(handles.imstruct(i).data);
   end
   fs = fs/(diff(framerange)+1);
   handles.background = fs;
   guidata(hObject,handles);
   imanalyze_guistate(handles);
   
   
% Backgroundframemin
% --- Executes during object creation, after setting all properties.
function backframemin_CreateFcn(hObject, eventdata, handles)
   if ispc
      set(hObject,'BackgroundColor','white');
   else
      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   end
   
function backframemin_Callback(hObject, eventdata, handles)
% Nothing to do
   
   
% Backgroundframemax
% --- Executes during object creation, after setting all properties.
function backframemax_CreateFcn(hObject, eventdata, handles)
   if ispc
      set(hObject,'BackgroundColor','white');
   else
      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   end
   
function backframemax_Callback(hObject, eventdata, handles)
% Nothing to do
   
   
% --- Executes on button press in showbackground.
function showbackground_Callback(hObject, eventdata, handles)
% Create a new figure and plot the background with the same colormap used
% for the raw image
   figure
   image(uint16(handles.background));
   colormap(handles.contrastindex);
   set(gca,'Visible','off');
   drawnow
   
   
% --- Executes on button press in subtractbackground.
function subtractbackground_Callback(hObject, eventdata, handles)
% If going into subtraction mode, check to see if the contrast range has
% been set. If not, try to give it sensible defaults
   if get(hObject,'Value')
      if isempty(handles.subcontrastrange)
         imgtmp = double(handles.imstruct(handles.currentframe).data) - handles.background;
         handles.subcontrastrange = [min(min(imgtmp)) max(max(imgtmp))];
      end
      guidata(hObject,handles);
   end
   imanalyze_showframe(handles);
   
   
% --- Executes on button press in adjustcontrast.
function adjustcontrast_Callback(hObject, eventdata, handles)
% First, obtain the values in the image (either raw or
% background-subtracted)
   imgtmp = double(handles.imstruct(handles.currentframe).data);
   rng = handles.rawcontrastrange;
   subtracting = 0;
   if get(handles.subtractbackground,'Value')
      imgtmp = imgtmp - handles.background;
      subtracting = 1;
      rng = handles.subcontrastrange;
   end
   % Now call the range-setting GUI
   [tmprange,tmpcs] = imrangegui(imgtmp,rng,handles.colorsaturated);
   if isempty(tmprange)
      return; % User hit cancel, skip the rest
   end
   handles.colorsaturated = tmpcs;
   if subtracting
      handles.subcontrastrange = tmprange;
   else
      tmprange = round(tmprange);
      handles.rawcontrastrange = tmprange;
      % jason: another hack to solve matlab color map problem on windows
      if(~isa(handles.imstruct(1).data(1), 'uint16'))
         handles.contrastindex = imrange_to_colormap(tmprange,handles.saturated(2)+1);
         if (handles.colorsaturated)
            % Color saturated pixels blue and red
            handles.contrastindex(1,:) = [0 0 1]; % Blue
            handles.contrastindex(end,:) = [1 0 0]; % Red
         end
         colormap(handles.contrastindex);
      end
      
      % jason: a hack to solve matlab colormap problem on windows 20040922
      if(isa(handles.imstruct(1).data(1), 'uint8'))
         colormap(handles.contrastindex(1:256,:));   
      end
      
      % jason: another hack to solve matlab color map problem on windows
      if(isa(handles.imstruct(1).data(1), 'uint16'))
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
   imanalyze_guistate(handles);
   
   
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
      imanalyze_showframe(handles);
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
   % Make sure colors are the same
   for i = 1:nregions
      set(hlines(i),'Color',get(handles.regionshcircles(i),'Color'));
   end
   xlabel('Frame');
   ylabel('Mean Pixel Intensity');
   %title('ROI Intensity');
   set(gcf, 'name', 'ROI Intensity');
   set(gca, 'xlim', [1 nframes-1]); % ignore the last frame
   if(~isempty(handles.stimSeqInFrame))
      % set(gca, 'xtick', handles.stimSeqInFrame(:,2)-handles.nFramesToDiscard);
      % Xs=handles.stimSeqInFrame(:,2)-handles.nFramesToDiscard;
      
      % set(gca, 'xgrid', 'on');
      set(gca, 'xminortick', 'on');
      Ys=get(gca, 'ylim');
      yOffset=(Ys(2)-Ys(1))/16;
      ttY=mean(Ys)-yOffset;
      for idx=1:size(handles.stimSeqInFrame,1)
	 ttX=handles.stimSeqInFrame(idx,2)-handles.nFramesToDiscard;
	 if(ttX <1) continue; end
         h=line([ttX, ttX], [Ys(1), Ys(2)]);
         set(h, 'color', 'red');
         set(h, 'linestyle', ':');
	 text(ttX, ttY+mod(idx,2)*2*yOffset, num2str(handles.stimSeqInFrame(idx,1)));
      end
   end
   
   % Save the data
   handles.intensity = mag;
   guidata(hObject,handles);
   
% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
   [pathstr,name,ext] = fileparts(handles.filename);
   [thefilename,pathname] = uiputfile([pathstr name '.mat'],'Save data to file...')
   regions = handles.regions;
   intensity = handles.intensity;
   filename = handles.filename;
   save([pathname thefilename],'filename','regions','intensity');
% Still to be written


% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
   handles.nFramesToDiscard=2; % set the number of frames to discard
   [handles.imstruct,handles.image,handles.filename, stimSeqInFrame] = imanalyzeload(handles);
   handles.stimSeqInFrame=stimSeqInFrame;
   % Now set a reasonable default contrasts--stretch the min/max of the image
   % into the available range of greyscales
   immin = min(min(handles.imstruct(1).data));
   immax = max(max(handles.imstruct(1).data));
   handles.rawcontrastrange = double([immin immax]);
   % jason: 20050321: a hack to load Andor camera's data
   if(~isa(handles.imstruct(1).data(1), 'uint16'))
      handles.contrastindex = imrange_to_colormap([immin immax],handles.saturated(2)+1);
      if (handles.colorsaturated)
         % Color saturated pixels blue and red
         handles.contrastindex(1,:) = [0 0 1]; % Blue
         handles.contrastindex(end,:) = [1 0 0]; % Red
      end
      colormap(handles.contrastindex);
   end
   
   % jason: a hack to solve matlab colormap problem on windows 20040922
   if(isa(handles.imstruct(1).data(1), 'uint8'))
      colormap(handles.contrastindex(1:256,:));   
   end
   
   % jason: another hack to solve matlab color map problem on windows
   if(isa(handles.imstruct(1).data(1), 'uint16'))
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
         set(gca, 'clim', handles.rawcontrastrange); % gca = handles.imaxes
   end

   
   % clear old region info:
   % delete(handles.regionshcircles);
   handles.regionshcircles=[];
   handles.regions=zeros(3,0);  % Rows are x,y,r
   
   guidata(hObject,handles);
   % Put GUI in consistent state
   imanalyze_guistate(handles);
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refresh_image(hObject, eventdata)
   handles = guidata(hObject);
   set(handles.image,'CData',[]); % This seems needed to force a redraw
   imanalyze_showframe(handles);

% Function to handle region selection
function regions_callback(hObject, eventdata)
   if(get(hObject, 'marker')=='*')
      set(hObject, 'marker', 'none');
   else
      set(hObject, 'marker', '*');
   end
   refresh_image(hObject, eventdata);
   
function imanalyze_selectstart(hObject, eventdata)
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
   imanalyze_guistate(handles);
   imanalyze_showframe(handles);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% A function to handle loading of a new image file, and updating the GUI
% accordingly.
function [imstruct,himg,imfilename, stimSeqInFrame] = imanalyzeload(handles)
   [imfilename,impathname] = uigetfile({'*.*';'*.tif';'*.stk';'*.mat'},'Pick an image stack');
   if ~isstr(imfilename)
      set(handles.imaxes,'CData',uint16(zeros(size(get(handles.imaxes,'CData')))));
      set(handles.subtractbackground,'Value',0);
      set([handles.showbackground handles.subtractbackground],'Enable',0);
      return
   end
   [pathstr,basename,extname] = fileparts(imfilename);
   if strcmp(extname,'.stk')
      [imstruct,nimages] = tiffread([impathname imfilename]);
      % todo: discard the first 2 frames
   elseif strcmp(extname,'.tif')
      tNumToDiscard=handles.nFramesToDiscard;
      iminfs =  imfinfo([impathname imfilename]);
      iminfs = iminfs(tNumToDiscard+1:end); % jason: discard the first tNumToDiscard frames 4831
      nimages = length(iminfs);
      for i = 1:nimages
          imstruct(i).data = imread([impathname imfilename],i+tNumToDiscard);
          imstruct(i).width = iminfs(i).Width;
          imstruct(i).height = iminfs(i).Height;
	  if(mod(i,5)==1 || i==nimages)
	     progress_bar(struct('progress', i, 'max', nimages, 'what', ['loading ' imfilename ' ...']));
          end 
      end % for,
   else
      % load raw data here using imphys
      tDataFile=fullfile(impathname, imfilename);
      tHeaderFile=[tDataFile '.txt'];
      if(~fileexist(tHeaderFile))
         error(['Cannot find header file: ' tHeaderFile]);
      end
      ip=imphysfrom2d(tDataFile, tHeaderFile);
      tNumToDiscard=handles.nFramesToDiscard;
      ip=ip(tNumToDiscard+1:end-1); % throw away the last frame also
      % Show an image to get crop region
      hCropRect = imselrect(imphysfetch(ip(1)));
      ip = imcroprect(ip,hCropRect);
      close(hCropRect);
      
      nimages = length(ip);
      for i = 1:nimages
          imstruct(i).data = imphysfetch(ip(i));
          imstruct(i).width = diff(ip(i).xrange)+1;
          imstruct(i).height = diff(ip(i).yrange)+1;
	  if(mod(i,10)==1 || i==nimages)
	     progress_bar(struct('progress', i, 'max', nimages, 'what', ['loading ' imfilename ' ...']));
	  end
      end % for,
      
      % load([impathname imfilename]);
      % imstruct = s;
      % nimages = length(imstruct);
   end
   
   % load the stimulus file if any
   descFileName=[impathname imfilename '.txt'];
   if(~fileexist(descFileName))
      descFileName=[impathname pathstr basename '.txt'];
      if(~fileexist(descFileName))
	 descFileName=[];
      end
   end
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
   axmin = [281 261];
   gap = 20;
   rcol = 141;
   axsize = max([axmin;imstruct(1).width imstruct(1).height]);
   if(handles.imanalyzeParams.mayShink)
      tMaxVisibleAxSize=get(0, 'screensize');
      tMaxVisibleAxSize=tMaxVisibleAxSize(3:4)-[20 100]-[3*gap+rcol,6*gap];
      tRatioHtoW=imstruct(1).height/imstruct(1).width;
      if(tMaxVisibleAxSize(1)*tRatioHtoW<=tMaxVisibleAxSize(2))
	 tMaxVisibleAxSize(2)=tMaxVisibleAxSize(1)*tRatioHtoW;
      else
	 tMaxVisibleAxSize(1)=tMaxVisibleAxSize(2)/tRatioHtoW;
      end
      axsize=min([axsize; tMaxVisibleAxSize]);
   end
   figsize = axsize + [3*gap+rcol,6*gap];
   figpos = get(gcbf,'Position');
   set(gcbf,'Position',[12 75 figsize]);
   % Adjust the position of all the controls
   % Do the controls on the right column
   set(handles.imaxes,'Position',[gap 5*gap axsize]);
   allrcolh = [handles.computebackground ...
               handles.text_usingframes ...
               handles.text_to ...
               handles.showbackground ...
               handles.subtractbackground ...
               handles.adjustcontrast ...
               handles.showregions ...
               handles.showintensities ...
               handles.savebutton ...
               handles.loadbutton];
       
   % jason: a hack to solve clipped display problem on windows
   oldaxsize=axsize; axsize(2)=axsize(2)/2+20;
   
   rcolvert = [0 1 2 3 5 7 11 12 14 16];
   for i = 1:length(allrcolh)
      pos = get(allrcolh(i),'Position');
      set(allrcolh(i),'Position',[2*gap+axsize(1) axsize(2)-(rcolvert(i)-3)*gap pos(3:4)]);
   end
   % Still right column; the background frame range edit boxes are aligned
   % differently, so require special treatment
   specialrcolh = [handles.backframemin handles.backframemax];
   for i = 1:length(specialrcolh)
      pos = get(specialrcolh(i),'Position');
      set(specialrcolh(i),'Position',[7*gap+axsize(1) axsize(2)-(i-3)*gap pos(3:4)]);
   end
   
   % jason: a hack to solve clipped display problem on windows
   axsize=oldaxsize;
   
   % Now do the controls below the image
   set(handles.slider,'Position',[gap 3.5*gap axsize(1) gap]);
   tOldPos=get(handles.playbutton, 'position');
   set(handles.playbutton, 'position', [tOldPos(1), 2*gap tOldPos(3), tOldPos(4)]);
   tOldPos=get(handles.editFps, 'position');
   set(handles.editFps, 'position', [tOldPos(1), 2*gap, tOldPos(3:4)]);
   tOldPos=get(handles.textFps, 'position');
   set(handles.textFps, 'position', [tOldPos(1), 2*gap, tOldPos(3:4)]);
   tOldPos=get(handles.btnDecFps, 'position');
   set(handles.btnDecFps, 'position', [tOldPos(1), 2*gap, tOldPos(3:4)]);
   tOldPos=get(handles.btnIncFps, 'position');
   set(handles.btnIncFps, 'position', [tOldPos(1), 2*gap+tOldPos(4)-2, tOldPos(3:4)]);
   
   set(handles.frametext,'Position',[axsize(1)-6*gap 2*gap 2.2*gap gap]);
   set(handles.framenumber,'Position',[axsize(1)-4*gap 2*gap 2*gap gap]);
   set(handles.frame_of_text,'Position',[axsize(1)-2*gap 2*gap 3*gap gap]);
   % Now set up the image data
   axes(handles.imaxes);
   himg = image(imstruct(1).data);
   set(handles.imaxes,'Visible','off');
   % hook the mouse down event of the image:
   set(himg,'EraseMode','none','ButtonDownFcn',@imanalyze_selectstart);
   % Update the GUI data
   set(handles.frame_of_text,'String',['  of ' num2str(nimages)]);
   set(handles.slider,'Min',1,'Max',nimages,'Value',1,'SliderStep',[1/(nimages-1) 0.1]);
   set(handles.backframemin,'String','1');
   % set(handles.backframemax,'String',num2str(nimages));
   if(isempty(stimSeqInFrame))
      set(handles.backframemax,'String', num2str(round(nimages/9)) );
   else
      tRow=find(stimSeqInFrame(:,1)>0);
      set(handles.backframemax, 'string', num2str(stimSeqInFrame(tRow(1),2)-2-handles.nFramesToDiscard));
   end
   set(handles.framenumber,'String','1');
   set(handles.figure1, 'name', [basename ' -- normal']);
   
   
function imanalyze_guistate(handles)
   allhandles = [handles.slider handles.playbutton handles.frametext handles.framenumber ...
                 handles.frame_of_text handles.computebackground ...
                 handles.text_usingframes handles.backframemin handles.text_to handles.backframemax ...
                 handles.showbackground handles.subtractbackground handles.adjustcontrast handles.showregions ...
                 handles.showintensities handles.savebutton handles.loadbutton];
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
         hon = [handles.playbutton handles.frametext handles.framenumber handles.frame_of_text handles.slider];
         if ~isempty(handles.background)
            hon = [hon handles.subtractbackground];
         end
         hoff = setdiff(allhandles,hon);
         set(hoff,'Enable','off');
         set(hon,'Enable','on');
         % Also set the playbutton text to "Stop"
         set(handles.playbutton,'String','Stop');
      else
         % We're not playing, but we do have a file loaded. Set up states on
         % buttons appropriately.
         set(allhandles,'Enable','on');
         if isempty(handles.background)
            set([handles.showbackground handles.subtractbackground],'Enable','off');
            set(handles.subtractbackground,'Value',0);
         end
         if isempty(handles.regions)
            set([handles.showregions handles.showintensities],'Enable','off');
            set(handles.showregions,'Value',0);
         end
         if (isempty(handles.background) && isempty(handles.regions))
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
   
function currentframe = imanalyze_playmovie(handles)
   tic;
   fps=str2num(get(handles.editFps, 'string'));
   nFramesPlayed =0;
   while (handles.currentframe>=1 && handles.currentframe <= length(handles.imstruct) ...
          && get(handles.playbutton,'Value'))
      imanalyze_showframe(handles);
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
   
   if(handles.currentframe>length(handles.imstruct))
      currentframe = length(handles.imstruct);
   elseif(handles.currentframe<1)
      currentframe = 1;
   end
      
   
function imanalyze_showframe(handles)
% Drawing is done with xor/none erasing, so first we have to turn off the
% regions if they are showing
   if get(handles.showregions,'Value')
      set(handles.regionshcircles,'Visible','off');
   end
   img = handles.imstruct(handles.currentframe).data;
   if get(handles.subtractbackground,'Value')
      % Subtract the background. This requires that we first identify saturated
      % pixels (if applicable), do the subtraction, and finally adjust the
      % contrast range
      lowindx = []; highindx = [];
      if (handles.colorsaturated)
         lowindx = find(img == handles.saturated(1));
         highindx = find(img == handles.saturated(2));
      end
      % Subtract the background
      imgtmp = double(img) - handles.background;
      % Now scale to the contrast range
      imgtmp = (imgtmp - handles.subcontrastrange(1)) * diff(handles.rawcontrastrange) / ...
               diff(handles.subcontrastrange) + handles.rawcontrastrange(1);
      % Convert back to uint16
      img = uint16(imgtmp);
      % Restore information about saturated pixels
      img(lowindx) = handles.saturated(1);
      img(highindx) = handles.saturated(2);
   end
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
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
   
% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
   usewhitebg = 1;
   if usewhitebg
      set(hObject,'BackgroundColor',[.9 .9 .9]);
   else
      set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
   end
   
   
% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
   
   
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
   userResponse=questdlg('Delete selected or all regions?', 'imanalyze', ...
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
   imanalyze_guistate(handles);
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
   tFilename=handles.filename;
   save([pathname thefilename], 'tRegions', 'tFilename', '-MAT');

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
   imanalyze_guistate(handles);
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


