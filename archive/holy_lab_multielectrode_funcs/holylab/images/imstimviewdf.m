function imstimviewdf(varargin)
% IMSTIMVIEWDF: show and analyze deltaf changes upon stimulation
%
% This function allows you to specify a number of IMPHYS structures which
% serve different roles. The valid roles are:
%   'raw': images (all of them) relative to which ROIs are to be defined
%   'movie': frames to use in playing movies (if desired)
%   'dfdraw': used for drawing ROIs
%   'dfcalc': measure changes in fluorescence upon stimulation (if
%     absent, uses dfdraw)
%   'dfthumb': thumbnails of drdraw (if absent, uses dfdraw)
%
% Note there must be a 1-1 correspondence between all the "df" types;
% dfcalc(i) must correspond to dfdraw(i), etc.  You will get
% unpredictable behavior if this is not true.
%   
% Syntax:
%   imstimviewdf(role1,value1,role2,value2,...)
%   imstimviewdf(...,options)
%
% Note that you can specify the same images for more than one role: e.g.,
% the following syntax is valid:
%   imstimviewdf('raw',raw,'movie',raw,'dfdraw',dfsmooth)
%
% options is a structure which may contain the following fields:
%   valvenum: a vector of valve numbers to show, or ':' for
%     all (default ':');
%   valvelabel: a cell array, valvelabel{j} is the label used for valve 
%     number valvenum(j) (by default uses the valvenum);
%   clim: color limits for contrast mapping (default: autoscales each
%     axis);
%   invert: if true, inverts the contrast (default: false)
%   extraframes: if set, continues the time associated with each valve
%     transition for 'extraframes' frames beyond the offset of the valve
%     (default: 0).
  
  % Parse arguments role/value arguments
  clim = [-0.1 .15];
  nrvpairs = floor(nargin/2);  % number of role/value pairs
  role = varargin(1:2:2*nrvpairs);
  value = varargin(2:2:2*nrvpairs);
  good_roles = {'raw','dfdraw','dfcalc','dfthumb','movie'};
  for i = 1:nrvpairs
    if ~isempty(strmatch(role{i},good_roles,'exact'))
      ip.(role{i}) = value{i};
    else
      error(['Unrecognized role ' role{i}]);
    end
  end
  if ~isfield(ip,'dfcalc')
    ip.dfcalc = ip.dfdraw;
  end
  if ~isfield(ip,'dfthumb')
    ip.dfthumb = ip.dfdraw;
  end

  % Parse options
  if (2*nrvpairs == nargin)
    options = struct;
  else
    options = varargin{end};
  end
  options = isvdf_options(options);
  

  % Collect the different stimuli together
  ustim = unique([ip.dfdraw.stimulus]);
  if (ischar(options.valvenum) & strmatch(options.valvenum,':'))
    vlvs = ustim;
  else
    vlvs = intersect(options.valvenum,ustim);
  end
  if isempty(vlvs)
    warning('No valves selected, exiting');
    return;
  end
  
  %
  % Show results for each stimulus valve in separate figures, plus one
  % summary figure
  %
  % First set up the summary figure
  dims_sum = CalcSubplotDims(length(vlvs));
  hfigsum = figure('NumberTitle','off','Name','Summary images');
  sz = size(ip.dfthumb(1).image);
  
  % a workaround of a bug(?) in vimage
  if(sz(end)==0)
     sz(end)=[];
  end
     
  haxsum = imsubplot(dims_sum(1),dims_sum(2),sz,0.15);
  
  if options.invert
    colormap(1-gray(256));
  else
    colormap(gray(256));
  end
  
  % save data as summury figure's appdata
  setappdata(hfigsum, 'ip', ip);
  
  valveFigures=zeros(1,length(vlvs))-1;
  
  % Loop over the individual valves, displaying each in a separate window
  % and then filling in the appropriate axis in the summary figure
  for i = 1:length(vlvs)
    %progress_bar(struct('progress',i,'max',length(vlvs),'what',...
    %                    ['Getting thumbnails']));
    if isfield(options,'valvelabel')
      indx = find(options.valvenum == vlvs(i));
      wintitle = options.valvelabel{indx};
    else
      wintitle = ['Valve ' num2str(vlvs(i))];
    end
    figCurVlv=figure('NumberTitle','off','Name',wintitle, 'toolbar', 'figure'); % BackingStore off?
    valveFigures(i)=figCurVlv;
    
    % save "parent" figure (the summary figure)
    setappdata(figCurVlv, 'figSummary', hfigsum);
    setappdata(figCurVlv, 'valve', vlvs(i));
    
    % Find the individual trials
    indx = find([ip.dfthumb.stimulus] == vlvs(i));
    % Create the axes
    dims = CalcSubplotDims(length(indx));
    tmp = imsubplot(dims(1),dims(2),sz, 0.2);
    delete(tmp(length(indx)+1:end));
    haxvlv{i} = tmp(1:length(indx));
    if options.invert
      colormap(1-gray(256))
    else
      colormap(gray(256))
    end
    
    setappdata(figCurVlv, 'allAxes', haxvlv{i});
    
    % Plot each trial
    for j = 1:length(indx)
      im = imphysfetch(ip.dfthumb(indx(j)));
      axesCur=haxvlv{i}(j);
      
      % save the "indx(j)" in axesCur's appdata
      setappdata(axesCur, 'index', indx(j));
      setappdata(axesCur, 'trial', j);
      
      axes(axesCur);
      if isfield(options,'clim')
        himg = imagesc(im,options.clim);
        clim = options.clim;
      else
        himg = imagesc(im);
      end
      if isa(ip.dfthumb(indx(j)).image,'vimage')
        setimagehandle(ip.dfthumb(indx(j)).image,himg);
      end
      set(gca,'Visible','off')
      %title(num2str(ip.dfthumb(indx(j)).stimulus))
      title(axesCur, ['trial #' num2str(j)], 'visible', 'on');
      
      %
      % Set up callbacks
      %
      %set(himg,'HitTest','off')
      cmenu = uicontextmenu;
      set(himg,'UIContextMenu',cmenu)
      % ROI drawing
      %set(himg,'ButtonDownFcn',{@roidraw,??});
      uimenu(cmenu,'Label','ROIs','Callback',{@fManipulateROIs, struct('axes', axesCur, 'options_imstimviewdf', options)});
              % or maybe set to call roidraw after de-thumbnailing?
      % (todo: to copy code from here)Movie of valve transition
      stacknum = [ip.dfthumb(indx(j)).stacknumbg ...
                  ip.dfthumb(indx(j)).stacknumfg];
      stacknum = [min(stacknum) max(stacknum)+options.extraframes];  % The range of frames
      setappdata(axesCur, 'movie_stacknum', stacknum); % indices to ip.movie and ip.raw
      setappdata(axesCur, 'movie_filename', ip.dfthumb(indx(j)).filename);
      if isfield(ip,'movie')
        movieindex = find([ip.movie.stacknum] >= stacknum(1) & ...
          [ip.movie.stacknum] <= stacknum(2));
        uimenu(cmenu,'Label','Movie of transition','Callback',{@mplayimphys,ip.movie,...
          struct('fps',3,'busymode','queue','frameindex',movieindex,...
          'clim',himg,'showstimulus',1)});
        % Movie since prev. valve transition
        if (j > 1)
          stacknum = [ip.dfthumb(indx(j-1)).stacknumbg ...
            ip.dfthumb(indx(j)).stacknumfg];
          stacknum = [min(stacknum) max(stacknum)];  % The range of frames
          movieindex = find([ip.movie.stacknum] >= stacknum(1) & ...
            [ip.movie.stacknum] <= stacknum(2));
          uimenu(cmenu,'Label','Movie since prev. transition',...
            'Callback',{@mplayimphys,ip.movie, ...
            struct('fps',30,'busymode','drop','showstimulus',1,...
            'frameindex',movieindex,'clim',himg)});
        end
      end
    end
    axis(haxvlv{i}, 'image');
    % Calculate the average response across trials, and display in a
    % summary window
    tmp = imphysfetch(ip.dfthumb(indx));
    if(~iscell(tmp))
       tmp={tmp};
    end
    im_mean = immean(tmp{:});
    axes(haxsum(i));
    if isfield(options,'clim')
      himg=imagesc(im_mean,options.clim);
    else
      himg=imagesc(im_mean);
    end
    set(haxsum(i),'Visible','off')
    title(wintitle,'Visible','on')
    setappdata(haxsum(i), 'relatedVlvFig', valveFigures(i) );
    install_mouse_event_handler(haxsum(i), 'up', @onClickOnAxes);
    
    cmenu = uicontextmenu;
    set(haxsum(i),'UIContextMenu',cmenu);
    uimenu(cmenu,'Label','Adjust contrast','Callback',{@fAdjustContrast, struct('axes', haxsum(i), 'himg', himg)});
    uimenu(cmenu,'Label','load ROI definitions from file','Callback',{@fLoadRoiDefsFromFile, struct('figSummary', hfigsum)});
    uimenu(cmenu,'Label','save ROI definitions to file','Callback',{@fSaveRoiDefsToFile, struct('figSummary', hfigsum)});
  end  % Loop over valves
  axis(haxsum(:), 'image'); % axis(haxsum, 'image'); should be fine but matlab has a bug.
  
  % save "child figures" in summury figure's appdata
  setappdata(hfigsum, 'valveFigures', valveFigures);

function result=onClickOnAxes(sender, event_args)
   if(is_button_down(sender, 'left'))
      relatedVlvFig=getappdata(sender, 'relatedVlvFig'); % sender: the axes
      if(ishandle(relatedVlvFig) && ~isempty(getappdata(relatedVlvFig, 'valve')))
         figure(relatedVlvFig);
      end
   else
      show_popup_menu(sender);
   end
   
   result=0;
  
function options = isvdf_options(options)
  if ~isfield(options,'valvenum')
    options.valvenum = ':';
  end
  if ~isfield(options,'invert')
    options.invert = 0;
  end
  if ~isfield(options,'extraframes')
    options.extraframes = 0;
  end

function fManipulateROIs(sender, event_data, args)
   manipulate_rois(args);
   
function fAdjustContrast(sender, event_data, args)
   curSummaryAxes=args.axes;
   hImg=args.himg;
   [newclim, newcs] = imrangegui(get(hImg, 'CData'), get(curSummaryAxes, 'clim'), 0);
   set(curSummaryAxes, 'clim', newclim);
   relatedVlvFig=getappdata(curSummaryAxes, 'relatedVlvFig');
   allAxes=getappdata(relatedVlvFig, 'allAxes');
   set(allAxes, 'clim', newclim);
   
function fLoadRoiDefsFromFile(sender, event_data, args)
   figSummary=args.figSummary;
   if(isunix)
      file_filter='*.mat';
   else
      file_filter='*.mat';
   end
   [filename, pathname] = uigetfile(file_filter, 'Pick a *.mat file');
   if(filename==0)
      return;
   end
   filename  =fullfile(pathname, filename);
   
   tt=load(filename);
   setappdata(figSummary, 'rois', tt.rois);
   
   
function fSaveRoiDefsToFile(sender, event_data, args)
   figSummary=args.figSummary;
   if(isunix)
      file_filter='*.mat';
   else
      file_filter='*.mat';
   end
   [filename, pathname] = uiputfile(file_filter, 'Pick a *.mat file');
   if(filename==0)
      return;
   end
   filename  =fullfile(pathname, filename);
   % todo: add ext name when necessary
   
   rois=getappdata(figSummary, 'rois');
   save(filename, 'rois', '-mat');
   
   
   