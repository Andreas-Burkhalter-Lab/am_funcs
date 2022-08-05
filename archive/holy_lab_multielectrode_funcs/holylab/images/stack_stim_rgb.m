function s = stack_stim_rgb(s)
% STACK_STIM_RGB: overlay stimulus responses in an rgb image
% This function allows you to inspect an image stack as either a grayscale
% or as an RGB image, where the intensity in any one of the 3 color
% channels is set by the deltaF/F for a particular stimulus.
%
% Syntax:
%   sout = stack_stim_rgb
% This version prompts the user to select a .imagine file for analysis. It
% returns a structure sout which can be useful for modifying the nature of
% the display (see below).
%
%   sout = stack_stim_rgb(struct('filename','vno1.imagine'))
% Allows you to specify the filename on the command line
%
%   sout = stack_stim_rgb(s)
% (s is a structure) allows you to control many things via the command
% line. Below is a complete list of fields.
%
% In addition to being able to control the display programmatically, you
% have a set of key commands at your disposal:
%    space: toggles between grayscale and rgb
%    r: set the stimulus displayed in the red channel
%    g: set the stimulus displayed in the green channel
%    b: set the stimulus displayed in the blue channel
%    t: set the trial number(s) displayed
%    up arrow: switch to higher frame number (move in z)
%    down arrow: switch to lower frame number (move in z)
%    f: go to a particular frame number
%    m: set a mark (click with mouse; hit return to abort)
%    d: delete a mark (")
%    v: toggle visibility of marks
%    z: set colorization metric
%    c: set clims
%    i: set bg and fg indices
%    e: export to stack_3d
%    s: save the current image stack to a vol3d model
%
% Here are the fields that can be used to programmatically control the
% display:
%   filename: the name of the .imagine file
%   hax: sets the destination axis of the plot.  If s is derived from sout
%     (the return value from an earlier call), this allows you to modify an
%     existing plot.
%   viewmode: 'gray' or 'rgb' (default 'gray'). The middle stack is chosen
%     for the grayscale plot
%   r, g, b: the stimulus name assigned to the red, green, and blue
%     channels respectively
%   trial (default 1): the trial number(s) used for the display. If you
%     specify more than one, it averages.
%   bgindex, fgindex: the frame numbers relative to valve onset used to
%     define the "background" and "foreground" for computing deltaF/F.
%     Default values are bgindex = [-2 -1] and fgindex = [0 1].
%   clim_raw, clim_dfof, clim_z: the "color limits" for the grayscale and deltaF/F
%     views, respectively. By default the raw view scales to the min/max,
%     and the dfof default is [0 0.2]. Default of clim_z is [0 1.5].
%   sigma (default [0 0 0]): smoothing (in pixels) applied before
%     calculating deltaF/F. [0 0 0] corresponds to no smoothing.
%   metric (default 'dfof'): describe metric used to evaluate voxel
%     response magnitude.  Currently supported: 'dfof','zscore'
%
% The defaults are tuned for performance (no spatial smoothing, no
% averaging over trials); in some circumstances you may prefer to adjust
% these to improve the appearance of the images.
%
% See also: stack_roi_rg.

% Copyright 2010 by Timothy E. Holy.

  %% Initialization
  if (nargin == 0)
    s = struct;
  end
  if isempty(fieldnames(s))
    % User supplied no arguments, so this must be a "new" call. Get the
    % filename.
    [filename,pathname] = uigetfile('*.imagine; *.mat','Select data file');
    if isequal(filename,0)
      return
    end
    if ~isempty(pathname)
      s.filename = [pathname filename];
    else
      s.filename = filename;
    end
  end
  smm = stackmm(s.filename);
  sz = smm.size;
  % Determine whether there is already stack data present; if not,
  % load it and calculate stimulus onset times
  if ~isfield(s,'hax')
    % The user did not supply an axis in which to draw, so create a new
    % figure and make it interactive
    figsz = min([1 1 sz(1:2); get(0,'ScreenSize')],[],1);
    %figsz(3:4) = min(figsz(3:4));
    hfig = figure('Position',figsz,...
      'KeyPressFcn',@ssrgb_keypress);
    s.hax = gca;
%     set(s.hax,'DataAspectRatio',[1 1 1]);%,...
      %'ButtonDownFcn',@ssrgb_btndown);
    set(hfig,'HandleVisibility','Callback');
  end
  if ~isappdata(s.hax,'smm')
    header = smm.header;
    onset = find_stimulus_start(header.stim_lookup);
    setappdata(s.hax,'smm',smm);
    setappdata(s.hax,'onset',onset);
    setappdata(s.hax,'stacksz',sz(1:3));
    % Extract the middle stack to show for the grayscale
    stk = smm(:,:,:,ceil(sz(end)/2));
    setappdata(s.hax,'gray',stk);
  end
  
  %% Set defaults
  channelname = {'r','g','b'};
  recalculate_channel = false(1,3);
  for i = 1:3
    recalculate_channel(i) = appdatadefault(s.hax,channelname{i},'',s);
  end
  recalculate_all = appdatadefault(s.hax,'sigma',[0 0 0],s);
  recalculate_all = recalculate_all | appdatadefault(s.hax,'bgindex',[-2 -1],s);
  recalculate_all = recalculate_all | appdatadefault(s.hax,'fgindex',[0 1],s);
  recalculate_all = recalculate_all | appdatadefault(s.hax,'trial',1,s);
  recalculate_all = recalculate_all | appdatadefault(s.hax,'metric','dfof',s);
  if recalculate_all
    recalculate_channel = true(1,3);
  end
  appdatadefault(s.hax,'viewmode','gray',s);
  appdatadefault(s.hax,'clim_raw',[],s);
  appdatadefault(s.hax,'clim_dfof',[0 0.2],s);
  appdatadefault(s.hax,'clim_z',[0 1.5],s);
  appdatadefault(s.hax,'showmarks',true,s);
  appdatadefault(s.hax,'nextlabel',1);
  if ~isappdata(s.hax,'him')
    him = image(zeros(sz(1:2),'single'),'Parent',s.hax);
    setappdata(s.hax,'him',him);
    set(s.hax,'DataAspectRatio',[1 1 1],...
      'TickDir','out',...
      'YDir','normal');
  end
  sz = getappdata(s.hax,'stacksz');
  appdatadefault(s.hax,'frame',ceil(sz(3)/2),s);
  
  %% Calculate the deltaf/f (as needed)
  ssrgb_recalculate_dfof(s.hax,recalculate_channel);
  
  %% Display the image
  ssrgb_display(s.hax);
end

%% Recalculate deltaf/f
function ssrgb_recalculate_dfof(hax,recalculate_channel)
  if any(recalculate_channel)
    sigma = getappdata(hax,'sigma');
    trial = getappdata(hax,'trial');
    onset = getappdata(hax,'onset');
    bgindex = getappdata(hax,'bgindex');
    fgindex = getappdata(hax,'fgindex');
    metric = getappdata(hax,'metric');
    sz = getappdata(hax,'stacksz');
    smm = getappdata(hax,'smm');
    header = smm.header;
    channelname = {'r','g','b'};
    for i = 1:3
      if recalculate_channel(i)
        stimname = getappdata(hax,channelname{i});
        dataname = [channelname{i} 'data'];
        if isempty(stimname)
          dat = [];
        else
          stimindex = findainb(stimname,header.stim_labels);
          stackindex = onset{stimindex}(trial);
          % Do the averaging in a loop to reduce the likelihood of memory
          % overflow
          bg = zeros(sz,'single');
          fg = zeros(sz,'single');
          if ~isempty(strmatch('dfof',metric));
              for j = 1:length(stackindex)
                  bg = bg + mean(single(smm(:,:,:,stackindex(j)+bgindex)),4);
                  fg = fg + mean(single(smm(:,:,:,stackindex(j)+fgindex)),4);
              end
              dat = imfilter_gaussian_mex(fg,sigma) ./ imfilter_gaussian_mex(bg,sigma) - 1;
          elseif ~isempty(strmatch('zscore',metric))
              dat = zeros(sz,'single');
              for j = 1:sz(3)
                  bstacks = [];
                  fstacks = [];
                  tmpdfm = [];
                  tmpdfstd = [];
                  for k = 1:length(bgindex)
                      bstacks(k,:) = stackindex'+bgindex(k);
                  end
                  for k = 1:length(fgindex)
                      fstacks(k,:) = stackindex'+fgindex(k);
                  end
                  for k = 1:length(stackindex)
                      tmpbasis(:,:,:,k) = mean(single(smm(:,:,j,fstacks(:,k))),4)./mean(single(smm(:,:,j,bstacks(:,k))),4)-1;
                  end
                  dat(:,:,j)=mean(tmpbasis,4)./(std(tmpbasis,0,4)+0.02);
              end
              clear tmpdfm tmpdfstd;
              dat = imfilter_gaussian(dat,sigma);      
          else
              dat = [];
          end
        end
        setappdata(hax,dataname,dat);
      end
    end
    % Recalculate the displayed value
    ssrgb_recalculate_mappedrgb(hax,recalculate_channel);
  end
end

%% Recalculate the displayed RGB value
function ssrgb_recalculate_mappedrgb(hax,recalculate_channel)
  metric = getappdata(hax,'metric');
  if ~any(recalculate_channel)
    return
  end
  channelname = {'r','g','b'};
  if ~isempty(strmatch('dfof',metric))
      clim = getappdata(hax,'clim_dfof');
  elseif ~isempty(strmatch('zscore',metric))
      clim = getappdata(hax,'clim_z');
  end
  slope = 1/diff(clim);
  offset = -slope*clim(1);
  for i = 1:3
    if recalculate_channel(i)
      dataname = [channelname{i} 'data'];
      mapname = [channelname{i} 'mapped'];
      dat = getappdata(hax,dataname);
      if ~isempty(dat)
          dat = slope*dat + offset;
          dat(dat < 0) = 0;
          dat(dat > 1) = 1;
          setappdata(hax,mapname,permute(dat,[2 1 3]));
      else
        setappdata(hax,mapname,[]);
      end
    end
  end
end

%% Display
function ssrgb_display(hax)
  him = getappdata(hax,'him');
  viewmode = getappdata(hax,'viewmode');
  showmarks = getappdata(hax,'showmarks');
  frame = getappdata(hax,'frame');
  titlestr = sprintf('Frame %d',frame);
  delete(findobj(hax,'type','line'));  % delete any old marks
  switch viewmode
    case 'gray'
      stk = getappdata(hax,'gray');
      set(him,'CData',stk(:,:,frame)','CDataMapping','scaled');
      clim = getappdata(hax,'clim_raw');
      if isempty(clim)
        set(hax,'CLimMode','auto')
      else
        set(hax,'CLim',clim);
      end
      colormap(hax,gray(256));
    case 'rgb'
      channelname = {'r','g','b'};
      sz = getappdata(hax,'stacksz');
      rgbc = cell(1,3);
      for i = 1:3
        tmp = getappdata(hax,[channelname{i} 'mapped']);
        if ~isempty(tmp)
          rgbc{i} = tmp(:,:,frame);
        else
          rgbc{i} = zeros(sz([2 1]),'single');
        end
      end
      rgb = cat(3,rgbc{:});
      set(him,'CData',rgb);
      % In RGB mode show the stimulus colors
      for i = 1:3
        name = getappdata(hax,channelname{i});
        if ~isempty(name)
          titlestr = [titlestr ', ' channelname{i} ': ' name]; %#ok<AGROW>
        end
      end
  end
  title(hax,titlestr);
  if showmarks
    ssrgb_showmarks(hax);
  end
end

%% Display markers
function ssrgb_showmarks(hax,frame)
  marks = getappdata(hax,'marks');
  if (nargin < 2)
    frame = getappdata(hax,'frame');
  end    
  if ~isempty(marks)
    pos = cat(1,marks.posInPixel);
    keepFlag = pos(:,3) == frame;
    pos = pos(keepFlag,:);
    if ~isempty(pos)
      line(pos(:,2),pos(:,1),...
        'Parent',hax,...
        'MarkerFaceColor','w',...
        'MarkerEdgeColor','k',...
        'MarkerSize',16,...
        'LineWidth',3,...
        'EraseMode','xor',...
        'Marker','x',...
        'LineStyle','none');
    end
  end
end
%% Handle key presses
function ssrgb_keypress(src,event)
  hax = findobj(src,'type','axes');
  switch event.Key
    case {'r','g','b'}
      % Set the stimulus for a given color channel
      smm = getappdata(hax,'smm');
      header = smm.header;
      lbl = header.stim_labels;
      lbl{end+1} = '';  % to turn off a color channel
      currentselection = findainb(getappdata(hax,event.Key),lbl);
      [selection,ok] = listdlg('PromptString',['Pick stimulus for ' event.Key ' channel'],...
        'ListString',lbl,...
        'InitialValue',currentselection,...
        'SelectionMode','single');
      if (ok == 0 || selection == currentselection)
        return
      end
      setappdata(hax,event.Key,lbl{selection});
      recalculate_channel = false(1,3);
      recalculate_channel(findainb(event.Key,{'r','g','b'})) = true;
      ssrgb_recalculate_dfof(hax,recalculate_channel);
      setappdata(hax,'viewmode','rgb');
      ssrgb_display(hax)
    case 'space'
      % Toggle between rgb and grayscale
      modes = {'gray','rgb'};
      viewmode = getappdata(hax,'viewmode');
      indx = findainb(viewmode,modes);
      setappdata(hax,'viewmode',modes{3-indx});
      ssrgb_display(hax);
    case 't'
      % Set the trial number(s)
      onset = getappdata(hax,'onset');
      l = cellfun(@length,onset);
      currentselection = getappdata(hax,'trial');
      [selection,ok] = listdlg('PromptString','Pick trial number(s)',...
        'ListString',cellstr(num2str((1:min(l))')),...
        'InitialValue',currentselection,...
        'SelectionMode','multiple');
      if (ok == 0 || isequal(selection,currentselection) || isempty(selection))
        return
      end
      setappdata(hax,'trial',selection)
      ssrgb_recalculate_dfof(hax,true(1,3));
      ssrgb_display(hax)
    case 'uparrow'
      % Increase the frame number
      frame = getappdata(hax,'frame');
      sz = getappdata(hax,'stacksz');
      if (frame < sz(end))
        setappdata(hax,'frame',frame+1);
        ssrgb_display(hax);
      end
    case 'downarrow'
      % Decrease the frame number
      frame = getappdata(hax,'frame');
      if (frame > 1)
        setappdata(hax,'frame',frame-1);
        ssrgb_display(hax);
      end
    case 'f'
      % Choose a particular frame number
      sz = getappdata(hax,'stacksz');
      answer = inputdlg(sprintf('Enter frame number (1:%d)',sz(end)));
      frame = getappdata(hax,'frame');
      newframe = round(str2double(answer));
      if isnan(newframe)
        return
      end
      if (newframe > sz(end))
        newframe = sz(end);
      elseif (newframe < 1)
        newframe = 1;
      end
      if (newframe ~= frame)
        setappdata(hax,'frame',newframe);
        ssrgb_display(hax);
      end
    case 'm'
      % Mark a point
      marks = getappdata(hax,'marks');
      nextlabel = getappdata(hax,'nextlabel');
      frame = getappdata(hax,'frame');
      [y,x] = ginput(1);  % switched because of permute/transpose
      if ~isempty(x)
        if ~isempty(marks)
          marks(end+1).posInPixel = [x y frame];
        else
          marks = struct('posInPixel',[x y frame]);
        end
        marks(end).label = nextlabel;
        setappdata(hax,'marks',marks);
        setappdata(hax,'nextlabel',nextlabel+1);
        ssrgb_display(hax);
      end
    case 'd'
      % Delete a marker
      marks = getappdata(hax,'marks');
      if isempty(marks)
        return
      end
      % Find the markers in the current frame
      frame = getappdata(hax,'frame');
      pos = cat(1,marks.posInPixel);
      keepFlag = pos(:,3) == frame;
      pos = pos(keepFlag,:);
      if isempty(pos)
        return
      end
      [y,x] = ginput(1);  % switched because of permute/transpose
      if ~isempty(x)
        % Find the closest marker
        [~,minIndex] = mindist([x;y],pos(:,1:2)');
        allIndex = find(keepFlag);
        % Delete
        marks(allIndex(minIndex)) = [];
        setappdata(hax,'marks',marks);
        ssrgb_display(hax);
      end
    case 'v'
      % Toggle showmarks
      showmarks = getappdata(hax,'showmarks');
      if showmarks
        delete(findobj(hax,'type','line'));
        setappdata(hax,'showmarks',false);
      else
        setappdata(hax,'showmarks',true);
        ssrgb_showmarks(hax);
      end
    case 'z'
      % Set the metric
      metric = getappdata(hax,'metric');
      l = {'dfof','zscore'}; % add more choices here as they are added
      [selection] = questdlg('Select preferred metric','Select metric',...
        l{1},l{2},...
        metric);
      if isempty(selection)
        return
      end
      if ~isempty(strmatch(selection,'zscore'))
          trial = getappdata(hax,'trial');
          onset = getappdata(hax,'onset');
          if length(trial)<2
              trial = 1:length(onset{1}); % this is a dumb way to do it
              setappdata(hax,'trial',trial);
          end
      end
      setappdata(hax,'metric',selection)
      ssrgb_recalculate_dfof(hax,true(1,3));
      ssrgb_display(hax);
    case 'c'
      clim_raw = getappdata(hax,'clim_raw');
      if isempty(clim_raw)
          clim_raw = get(hax,'clim');
      end
      temp = inputdlg({'raw clim low:','raw clim high:'},...
                       'raw clims',...
                       [1; 1],...
                       {num2str(clim_raw(1)), num2str(clim_raw(2))});
      if isempty(temp); return; end;
      for i = 1:length(temp)
          if ~isempty(temp{i})
              clim_raw(i) = str2num(temp{i});
          end
      end;  
      setappdata(hax,'clim_raw',clim_raw);
      clim_dfof = getappdata(hax,'clim_dfof');
      temp = inputdlg({'dfof clim low:','dfof clim high:'},...
                       'dfof clims',...
                       [1; 1],...
                       {num2str(clim_dfof(1)), num2str(clim_dfof(2))});
      if isempty(temp); return; end;
      for i = 1:length(temp)
          if ~isempty(temp{i})
              clim_dfof(i) = str2num(temp{i});
          end
      end;
      setappdata(hax,'clim_dfof',clim_dfof);
      clim_z = getappdata(hax,'clim_z');
      temp = inputdlg({'zscore clim low:','zscore clim high:'},...
                       'zscore clims',...
                       [1; 1],...
                       {num2str(clim_z(1)), num2str(clim_z(2))});
      if isempty(temp); return; end;
      for i = 1:length(temp)
          if ~isempty(temp{i})
              clim_z(i) = str2num(temp{i});
          end
      end;
      setappdata(hax,'clim_z',clim_z);
      ssrgb_recalculate_mappedrgb(hax,[1 1 1]);
      ssrgb_display(hax);
    case 'i'
      bgindex = getappdata(hax,'bgindex');
      bginput = [min(bgindex) max(bgindex)];
      fgindex = getappdata(hax,'fgindex');
      fginput = [min(fgindex) max(fgindex)];
      temp = inputdlg({'background first stack:','background last stack:',...
                       'test first stack:','test last stack:'},...
                       'Stacks to analyze',...
                       [1; 1; 1; 1],...
                       {num2str(bginput(1)),num2str(bginput(2)),num2str(fginput(1)),num2str(fginput(2))});
      if isempty(temp); return; end;             
      for i = 1:2
          bginput(i) = str2num(temp{i});
      end;  
      for i = 1:2
          fginput(i) = str2num(temp{i+2});
      end;
      if diff(bginput)<0 || diff(fginput)<0 
          errordlg('bad stack values given');
          return; 
      end;
      bgindex = bginput(1):bginput(2);
      fgindex = fginput(1):fginput(2);
      setappdata(hax,'bgindex',bgindex);
      setappdata(hax,'fgindex',fgindex);
      ssrgb_recalculate_dfof(hax,true(1,3));
      ssrgb_display(hax);
    case 'e'
      [model sz0] = prep_vol3d_model(hax);
      
      newfig = figure('position',[1 1 sz0(1:2)],'color','k','renderer','opengl');
      model.parent = gca;
      set(gca,'color','k');
      h = vol3d(model);
      daspect([1 1 1])
      camup([1 0 0])
      axis tight
      set(gca,'Color','k','CLim',[0 1])
      
%       % TODO set up scalebars
%       scalebarpos = [250 -250 -200];
%       scalebardims = [-100 100 0];
    case 's'
      [output.file output.path output.fidx] = uiputfile('.mat','Select vol3d model save file name');
      if (output.file~=0)
          model = prep_vol3d_model(hax);
          save([output.path output.file],'model');
      end
    otherwise
        
          
      display(event.Key)
  end
  end

function imout = s3d_shrink(im,max_size)
  imout = im;
  restrictFlag = size(im) > max_size;
  while any(restrictFlag)
    imout = array_restrict(imout,restrictFlag);
    restrictFlag = size(imout) > max_size;
  end
end

function [model sz0] = prep_vol3d_model(hax)
  pixel_spacing = getappdata(hax,'pixel_spacing');
  if isempty(pixel_spacing)
      % gather input on pixel spacing from user
      temp = inputdlg({'um/pix (x):','um/pix (y):', 'um/pix (z)'},...
          'Pixel Dimensions (in um)',...
          [1; 1; 1],...
          {'0.71', '0.71', '8'});
      if isempty(temp); return; end;
      for i = 1:length(temp)
          pixel_spacing(i) = str2num(temp{i});
      end;
      setappdata(hax,'pixel_spacing',pixel_spacing);
  end
  
  sz = getappdata(hax,'stacksz');
  max_size = getappdata(hax,'max_size');
  if isempty(max_size)
      temp = inputdlg({'max_size (x):','max_size (y):', 'max_size (z)'},...
          'Max pixels',...
          [1; 1; 1],...
          {'512', '512', num2str(sz(3))});
      if isempty(temp); return; end;
      for i = 1:length(temp)
          max_size(i)=str2num(temp{i});
      end;
      setappdata(hax,'max_size',max_size);
  end
      
  smm = getappdata(hax,'smm');
  n_stacks = smm.size; n_stacks = n_stacks(4);
  base_stack = getappdata(hax,'base_stack');
  if isempty(base_stack)
      temp = inputdlg({'base_stack'},...
          'Choose base stack',...
          [1],...
          {num2str(floor(n_stacks/2))});
      if isempty(temp); return; end;
      for i = 1:length(temp)
          base_stack(i)=str2num(temp{i});
      end;
      setappdata(hax,'base_stack',base_stack);
  end
  
  % assemble parameters for vol3d
  % CDATA
  colorchans = {'r','g','b'};
  for i = 1:length(colorchans)
      temp = getappdata(hax,[colorchans{i} 'mapped']);
      if ~isempty(temp)
          model.cdata(:,:,:,i) = s3d_shrink(temp,max_size);
      else
          temp = zeros(sz,'single');
          model.cdata(:,:,:,i) = s3d_shrink(temp,max_size);
      end
  end
  
  % set up alpha mask
  temp = model.cdata;
  nanmask = find(isnan(temp(:,:,:,1)));
  for i = 1:3
      tmp = temp(:,:,:,i);
      tmp(nanmask) = 0;
      temp(:,:,:,i)=tmp;
  end; clear tmp;
  noactivity = find(sum(temp,4)==0);
  background_alpha = 0.02;
  active_alpha_min = 0.05;
  active_alpha_max = 0.9;
  model.alpha = max(model.cdata,[],4);
  model.alpha(model.alpha>0)=model.alpha(model.alpha>0)+active_alpha_min;
  model.alpha(model.alpha>active_alpha_max)=active_alpha_max;
  model.alpha(nanmask)=0;
  
  rawimg = s3d_shrink(single(smm(:,:,:,base_stack)),max_size);
  rawimg(isnan(rawimg))=0;
  rawimg = rawimg-median(median(median(rawimg,3),2),1);
  rawimg = rawimg/max(max(max(rawimg,[],3),[],2),[],1);
  for i = 1:sz(3)
      tmp(:,:,i)=rot90(rawimg(:,end:-1:1,i),1);
  end
  rawimg = tmp; clear tmp;
  for i = 1:3
      tmp = model.cdata(:,:,:,i);
      tmp(noactivity)=rawimg(noactivity);
      model.cdata(:,:,:,i)=tmp;
  end;
  model.alpha(noactivity)=background_alpha;
  model.alpha(rawimg<0.05)=0;
  
  sz0 = size(model.cdata(:,:,:,1));
  pixel_spacing = pixel_spacing .* (sz0 ./ sz);
  coordlim = pixel_spacing.*sz;
  model.xdata = [-1 1]*coordlim(1)/2;
  model.ydata = [-1 1]*coordlim(2)/2;
  model.zdata = [-1 1]*coordlim(3)/2;
  model.handles = [];
  model.texture = '3d';
  model.cdata = double(model.cdata);
end