function stackviewer_wm(smm)
% stackviewer_wm: volumetric image data browser using the wheel mouse
%
% Syntax:
%   stackviewer_wm(smm)
% where smm is either a 3- or 4-d array or a stackmm object.
%
% Controls:
%   Right-click on image: opens a context menu that allows you to adjust
%     contrast, go to a particular frame/stack, etc.
%   scroll (wheel mouse): move forward/backward in time
%   shift-scroll: move up/down in frame # (within a stack)
%
% When not using java:
%   up/down arrows: move up/down in frame # (within a stack)
%
% See also: stackviewer.

% Future extensions:
%    Type in a stack # and go directly to it, or maybe even better browse
%      to a particular stimulus via a popup menu?

% Copyright 2011 by Timothy E. Holy

  tryjava = usejava('jvm');

  hfig = figure('Tag','figMain');
  pos = [1 1];       % frame# stack#
  
  % Plot the frame and set up handles to all relevant objects
  im = smm(:,:,pos(1),pos(2));
  clim = [0 max(im(:))];
  if any(isnan(clim))
    % Read a whole stack
    stk = smm(:,:,:,pos(2));
    clim = [0 max(stk(:))];
    if any(isnan(clim))
      clim = [0 1000];
    end
  end
  himg = imshowsc(im,clim);
  set(himg,'Tag','imFrame');
  hax = get(himg,'Parent');
  set(hax,'Tag','axesFrame');
  htit = title(sprintf('Frame %d, stack %d',pos));
  set(htit,'Tag','textTitle');
%   uicontrol('Style','pushbutton','Position',[20 20 80 20],'String','Contrast...','Tag','buttonCLim','Callback',@swm_clim);
  handles = guihandles(hfig);
  guidata(hfig,handles);
  
  % Store application data
  setappdata(hfig,'position',pos);
  setappdata(hfig,'smm',smm);
  if isnumeric(smm)
    % It's a numeric array, so we can get the size directly
    sz = size(smm);
    sz(end+1:4) = 1;
    setappdata(hfig,'sz',sz);
  else
    % We assume this is a stackmm object
    sz = smm.size;
    setappdata(hfig,'sz',sz);
    % Store stimulus info
    header = smm.header;
    if (length(header.stim_lookup) >= sz(4))
      setappdata(hfig,'stim_lookup',header.stim_lookup);
      setappdata(hfig,'stim_labels',header.stim_labels);
    end
  end
  
  % Set up the context menu
  hcm = uicontextmenu('Parent',hfig);
  set(himg,'UIContextMenu',hcm);
  hcontrast = uimenu(hcm,'Label','Contrast...','Callback',@swm_clim);
  hgoto = uimenu(hcm,'Label','Go to...','Callback',@swm_goto);
%   hrotneg = uimenu(hcm,'Label','Rotate clockwise','Callback',@swm_rotneg);
%   hrotpos = uimenu(hcm,'Label','Rotate counterclockwise','Callback',@swm_rotpos);
  hflipx = uimenu(hcm,'Label','Flip horizontal','Callback',@swm_flipx);
  hflipy = uimenu(hcm,'Label','Flip vertical','Callback',@swm_flipy);
%   set([hcm hcontrast hgoto hrotneg hrotpos hflipx hflipy],'HandleVisibility','off');
  set([hcm hcontrast hgoto hflipx hflipy],'HandleVisibility','off');
  
  
  % Set up the interaction
  if tryjava
    % This version uses java, including the JavaFrame property that may be
    % disabled in future releases
    jf = get(handle(hfig),'JavaFrame');
    jfh = handle(jf,'CallbackProperties');
    jax = get(jfh,'AxisComponent');
    jaxh = handle(jax,'CallbackProperties');
    set(jaxh,'MouseWheelMovedCallback',@(obj,evt) swm_wheel(obj,evt,hfig))
  else
    % Version that doesn't use java, and uses the up/down buttons for
    % changing the frame
    set(hfig,'WindowScrollWheelFcn',@swm_scroll,...
      'KeyPressFcn',@swm_key)
  end
end

% Function to update the image and "location" text
function swm_update(hfig,smm,pos)
  handles = guidata(hfig);
  titlestr = sprintf('Frame %d, stack %d',pos);
  if isappdata(hfig,'stim_lookup')
    stim_lookup = getappdata(hfig,'stim_lookup');
    stim_labels = getappdata(hfig,'stim_labels');
    if (stim_lookup(pos(2)) == 0)
      stimstr = 'Flush';
    else
      stimstr = stim_labels{stim_lookup(pos(2))};
    end
    titlestr = [titlestr ', ' stimstr];
  end
  set(handles.textTitle,'String',titlestr);
  im = smm(:,:,pos(1),pos(2));
  set(handles.imFrame,'CData',im);
end

%% Context menu functions
% Contrast...
function swm_clim(obj,~)
  hfig = ancestor(obj,'figure');
  handles = guidata(hfig);
  clim = get(handles.axesFrame,'CLim');
  im = get(handles.imFrame,'CData');
  clim = imrangegui(im,clim);
  if ~isempty(clim)
    set(handles.axesFrame,'CLim',clim);
  end
end

% Go to...
function swm_goto(obj,~)
  hfig = ancestor(obj,'figure');
  pos = getappdata(hfig,'position');
  sz = getappdata(hfig,'sz');
  response = inputdlg({'Frame #','Stack #'},'Choose location',1,{num2str(pos(1)),num2str(pos(2))});
  pos = [str2double(response{1}), str2double(response{2})];
  setappdata(hfig,'position',pos)
  smm = getappdata(hfig,'smm');
  swm_update(hfig,smm,pos)
end

% Rotate
function swm_rotneg(obj,~)
  hfig = ancestor(obj,'figure');
  handles = guidata(hfig);
  [az,el] = view(handles.axesFrame);
  az = az-90;
  view(handles.axesFrame,az,el);
end
function swm_rotpos(obj,~)
  hfig = ancestor(obj,'figure');
  handles = guidata(hfig);
  [az,el] = view(handles.axesFrame);
  az = az+90;
  view(handles.axesFrame,az,el);
end

% Flip
function swm_flipx(obj,~)
  hfig = ancestor(obj,'figure');
  handles = guidata(hfig);
  dirstr = get(handles.axesFrame,'XDir');
  if strcmp(dirstr,'normal')
    set(handles.axesFrame,'XDir','reverse');
  else
    set(handles.axesFrame,'XDir','normal');
  end
end
function swm_flipy(obj,~)
  hfig = ancestor(obj,'figure');
  handles = guidata(hfig);
  dirstr = get(handles.axesFrame,'YDir');
  if strcmp(dirstr,'normal')
    set(handles.axesFrame,'YDir','reverse');
  else
    set(handles.axesFrame,'YDir','normal');
  end
end

%% Navigation
% Function handling the wheel mouse when using java
function swm_wheel(~,evt,hfig)
  rotdir = get(evt,'WheelRotation');
  shiftdown = get(evt,'ShiftDown');
  sz = getappdata(hfig,'sz');
  pos = getappdata(hfig,'position');
  if strcmp(shiftdown,'on')
    pos(1) = pos(1)+rotdir;
  else
    pos(2) = pos(2)+rotdir;
  end
  pos = max([1 1],min(sz(3:4),pos));
  setappdata(hfig,'position',pos)
  smm = getappdata(hfig,'smm');
  swm_update(hfig,smm,pos)
end



% Non-java scroll method to change the stack #
function swm_scroll(src,evt)
  hfig = get_parent_fig(src);
  smm = getappdata(hfig,'smm');
  pos = getappdata(hfig,'position');
  pos(2) = pos(2) + evt.VerticalScrollCount;
  if (pos(2) < 1)
    pos(2) = 1;
  end
  sz = getappdata(hfig,'sz');
  if (pos(2) > sz(4))
    pos(2) = sz(4);
  end
  setappdata(hfig,'position',pos)
  swm_update(hfig,smm,pos)
end

% Non-java method to change the frame #
function swm_key(src,evt)
  hfig = get_parent_fig(src);
  pos = getappdata(hfig,'position');
  switch evt.Key
    case 'downarrow'
      smm = getappdata(hfig,'smm');
      pos(1) = pos(1) - 1;
      if (pos(1) < 1)
        pos(1) = 1;
      end
      setappdata(hfig,'position',pos)
      swm_update(hfig,smm,pos);
    case 'uparrow'
      smm = getappdata(hfig,'smm');
      sz = getappdata(hfig,'sz');
      pos(1) = pos(1) + 1;
      if (pos(1) > sz(3))
        pos(1) = sz(3);
      end
      setappdata(hfig,'position',pos)
      swm_update(hfig,smm,pos);
  end
end