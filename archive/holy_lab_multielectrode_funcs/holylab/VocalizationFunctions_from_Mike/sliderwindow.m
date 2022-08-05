function hfigo = sliderwindow(action,hfig)
% SLIDERWINDOW: create a window for interactive plot zooming
% Syntax:
%   sliderwindow(axishandles)
%   sliderwindow(axishandles,options)
%   hfig = sliderwindow(...)
%
%   These axes become controlled by the new window, the slider window,
%   which allows the user to zoom in on particular portions of the axes.
%   The slider window displays the objects contained in axishandles(1),
%   along with a selection rectangle, and intercepts both mouse and key
%   commands to control the display limits of the original axes.
%
%   A single axis may be manipulated by supplying a single axis handle;
%   several axes may be simultaneously manipulated by supplying a vector
%   of axis handles.  When several axes are supplied, they will generally
%   have the same range along the x-axis, and this x-coordinate will
%   refer to a common set of events.  For example, one axis could be
%   displaying a measured signal, and the other axis the spectrogram of that
%   signal.  This function allows the user to zoom in on specific time
%   periods in both axes simultaneously.
%
%   If axishandles is a scalar, then click-dragging a selection rectangle
%   specifies the X- and Y-ranges in the original plot.
%   Shift-click-dragging the rectangle results in only X-axis zooming.
%
%   If axishandles is a vector, all axes zoom simultaneously.  X-zooming
%   is enabled on all axes, but only axishandles(1) can allow Y-zooming.
%
%   The selection rectangle itself can be dragged; holding down the shift
%   key before dragging permits only x-translation.
%
%   The "Select All" button makes the selection rectangle encompass the
%   whole slider window, restoring the axes to their original conditions.
%
%   In addition to using the mouse, several key commands can be used to
%     control the axis limits:
%        The arrow keys shift the selection rectangle horizontally by its
%          full width. 
%        The keys '<' and '>' shift the selection rectangle horizontally
%          by 10% of its width.
%
%   sliderwindow(...,options) allows one to customize the behavior of
%   this function.  options is a structure, in which the following fields
%   may be set: 
%       position:  the position for the slider figure window (in pixels).
%         A default value is chosen if this parameter is absent or empty.
%       restoreonclose: if true, the view axes are returned to their 
%         original limits when the slider window is closed.  (Default: false.)
%       axisupdatefcn: a function handle (or cell array of function
%         handles, if axishandles is a vector) used to update the axes
%         with new axis limits.  If absent, the axes are updated by
%         set(axishandles,'XLim',...) calls.  A user-supplied function
%         must accept calls of the form
%               axisupdatefcn(axishandle,'XLim',xlim)
%               axisupdatefcn(axishandle,'XLim',xlim,'YLim',ylim)
%         where
%               axishandle is the handle of the axis
%               xlim/ylim contain the new limits (as a 2-vector)
%         This option is useful in cases where the displayed data
%         must be calculated from other data, and the calculation
%         depends on the axis limits.  Examples might include the
%         display of data which consume a great deal of memory,
%         e.g., displaying a sparse matrix as an image.  See
%         SETSPIMG as an example.
%         SET itself may be used for axes which require no special
%         handling.
%   future: allow multiple axes to appear in the slider window? (control
%       through option field)
%
% hfig = sliderwindow(...) returns the figure handle of the interactive
% plot window.
%
% See also: ZOOM, SETSPIMG.

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>

% Determine whether this is the initial call or a callback service
if all(ishandle(action))
  % Initial call, set up the slider window
  hax = action;
  if (nargin > 1)
    swoptions = hfig;
  else
    swoptions = struct;
  end
  if (~isfield(swoptions,'position') | isempty(swoptions.position))
    ss = get(0,'ScreenSize');
    swoptions.position = ss .* [1 1 1 0.2];
  end
  hfig = figure('Position',swoptions.position,...
                'ButtonDownFcn','sliderwindow Select',...
                'KeyPressFcn','sliderwindow Key',...
                'HandleVisibility','callback');
  if (isfield(swoptions,'restoreonclose') & swoptions.restoreonclose)
    set(hfig,'DeleteFcn','sliderwindow RestoreAndKill');
  end
  if isfield(swoptions,'axisupdatefcn')
    setappdata(hfig,'axupdtfcn',swoptions.axisupdatefcn);
  else
    ad = cell(1,length(hax));
    for i = 1:length(ad)
      ad{i} = @set;
    end
    setappdata(hfig,'axupdtfcn',ad);
  end
  hnewax = axes('Parent',hfig,'Position',[0.05 0.11 0.8 0.78]);
  % Set up the buttons
  hbsel = uicontrol(hfig,'Style','PushButton','String','Select All',...
                    'units','normalized','Position',[0.88 0.7 0.09 0.19],...
                    'Callback','sliderwindow SelectAll');
  hbapply = uicontrol(hfig,'Style','PushButton','String','Apply',...
                      'units','normalized',...
                      'Position',[0.88 0.41 0.04 0.19], ...
                      'Callback','sliderwindow PlotTop'); 
  hbnew = uicontrol(hfig,'Style','PushButton','String','New',...
                    'units','normalized',...
                    'Position',[0.93 0.41 0.04 0.19],...
                    'Callback','sliderwindow New');
  icons = swicons;
  lmarg = 0.008;
  hbleft = uicontrol(hfig,'Style','PushButton',...
                     'CData',icons.left,...
                     'units','normalized',...
                     'Position',[0.86+lmarg 0.11 0.025 0.19],...
                     'Callback','sliderwindow Left');
  hbleftsm = uicontrol(hfig,'Style','PushButton',...
                       'CData',icons.leftsm,...
                       'units','normalized',...
                       'Position',[0.89+lmarg 0.11 0.025 0.19],...
                       'Callback','sliderwindow Leftsm');
  hbrightsm = uicontrol(hfig,'Style','PushButton',...
                        'CData',icons.rightsm,...
                        'units','normalized',...
                        'Position',[0.92+lmarg 0.11 0.025 0.19],...
                        'Callback','sliderwindow Rightsm');
  hbright = uicontrol(hfig,'Style','PushButton',...
                      'CData',icons.right,...
                      'units','normalized',...
                      'Position',[0.95+lmarg 0.11 0.025 0.19],...
                      'Callback','sliderwindow Right');

  % For optimal performance, one has to be careful about the order
  % in which the objects appear in the list of children of the slider axes.
  % Otherwise graphs that take a long time to draw can get drawn twice!
  % The strategy: make the selection box the first of the children, and
  % then add on the children of the source axis.
  hc = get(hax(1),'Children');        % Get all children
  lims = get(hax(1),{'XLim','YLim'});
  xlim = lims{1}; ylim = lims{2};
  % Set & plot the selection rectangle
  selrectx = [xlim(1) xlim(2) xlim(2) xlim(1) xlim(1)];
  selrecty = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)];
  hselrect = line(selrectx,selrecty,'LineStyle',':','Color','k',...
                  'ButtonDownFcn','sliderwindow Slide',...
                  'Tag','HSelRect',...
                  'Parent',hnewax,...
                  'EraseMode','xor');
  % Remember which axes to update
  setappdata(hfig,'UpdAx',hax);
  setappdata(hfig,'FullXLim',xlim);
  setappdata(hfig,'FullYLim',ylim);
  % Now copy all the children of the original axis to the new one
  hcnew = copyobj(hc(end:-1:1),hnewax);
  set([hcnew;hnewax],'HitTest','off');
  set(hnewax,'XLim',xlim,'YLim',ylim);

  if nargout
    hfigo = hfig;
  end
  return
elseif ischar(action)
  % Callback service
  if (nargin < 2)
    hfig = gcbf;
  end
  switch(action)
   % This does the actual rendering
   case 'PlotTop'
    hselrect = findobj(hfig,'Tag','HSelRect');
    xd = get(hselrect,'XData');
    yd = get(hselrect,'YData');
    xlim = [min(xd) max(xd)];
    ylim = [min(yd) max(yd)];
    if (diff(xlim)>0 & diff(ylim)>0)
      haxupd = getappdata(hfig,'UpdAx');
      axupdtfcn = getappdata(hfig,'axupdtfcn');
      if all(ishandle(haxupd))
        % Adjust y scaling on only the first axis
        if iscell(axupdtfcn)
          feval(axupdtfcn{1},haxupd(1),'XLim',xlim,'YLim',ylim);
        else
          feval(axupdtfcn,haxupd(1),'XLim',xlim,'YLim',ylim);
        end
        % Adjust x scaling on all
        for i = 2:length(haxupd)
          feval(axupdtfcn{i},haxupd(i),'XLim',xlim);
        end
      end
    end
   % This creates a selection rectangle
   case 'Select',
    hselrect = findobj(hfig,'Tag','HSelRect');
    theRect = GetSelRect;
    % Check to see if selection rectangle is outside the axis
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    if (theRect(1) > xlim(2) | theRect(1)+theRect(3) < xlim(1) ...
        | theRect(2) > ylim(2) | theRect(2)+theRect(4) < ylim(1))
      return;
    end
    % Check to see if it was just a single click
    if ~any(theRect(3:4))  % See if zero width, height
      return
    end
    theRect(1) = min(max(theRect(1),xlim(1)),xlim(2));
    theRect(3) = min(theRect(3),xlim(2)-theRect(1)); % width
    xd = [0 theRect(3) theRect(3) 0 0] + theRect(1);
    set(hselrect,'XData',xd);
    seltype = get(hfig,'SelectionType');
    if (~strcmp(seltype,'extend'))
      theRect(2) = max(theRect(2),ylim(1));
      theRect(4) = min(theRect(4),ylim(2)-theRect(2));
      yd = [0 0 theRect(4) theRect(4) 0] + theRect(2);
      set(hselrect,'YData',yd);
    end
    sliderwindow('PlotTop',hfig);
   % The next 3 allow dragging the selection rectangle
   case 'Slide',
    set(hfig,'WindowButtonMotionFcn','sliderwindow Move');
    set(hfig,'WindowButtonUpFcn','sliderwindow Stop');
    % Compute the relative position of the pointer and the selection rectange,
    % and store as offset
    currPt = get(gca,'CurrentPoint');
    hrect = findobj(hfig,'Tag','HSelRect');
    xd = get(hrect,'XData');
    yd = get(hrect,'YData');
    offset = currPt(1,1:2)-[min(xd) min(yd)];
    setappdata(hfig,'Offset',offset);
   case 'Move',
    currPt = get(gca,'CurrentPoint');
    hrect = findobj(hfig,'Tag','HSelRect');
    offset = getappdata(hfig,'Offset');
    %offset = [0 0];
    xd = get(hrect,'XData');
    xlimrect = [min(xd) max(xd)];
    % Make sure stay in bounds
    xlimabs = get(gca,'XLim');
    width = min(xlimrect(2)-xlimrect(1),xlimabs(2)-xlimabs(1));
    x = max(currPt(1,1) - offset(1),xlimabs(1));
    x = min(x,xlimabs(2)-width);
    xd = [0 width width 0 0] + x;
    set(hrect,'XData',xd);
    seltype = get(hfig,'SelectionType');
    if (~strcmp(seltype,'extend'))
      yd = get(hrect,'YData');
      ylimrect = [min(yd) max(yd)];
      % Make sure stay in bounds
      ylimabs = get(gca,'YLim');
      height = min(ylimrect(2)-ylimrect(1),ylimabs(2)-ylimabs(1));
      y = max(currPt(1,2) - offset(2),ylimabs(1));
      y = min(y,ylimabs(2)-height);
      yd = [0 0 height height 0] + y;
      set(hrect,'YData',yd);
    end
   case 'Stop',
    set(hfig,'WindowButtonMotionFcn','');
    set(hfig,'WindowButtonUpFcn','');
    sliderwindow('PlotTop',hfig);
   % Handle paging (left & right)
   case {'Left','Leftsm','Rightsm','Right'}
    hselrect = findobj(hfig,'Tag','HSelRect');
    xd = get(hselrect,'XData');
    if strcmp(action,'Right')
      shift = xd(2)-xd(1);
    elseif strcmp(action,'Left')
      shift = -(xd(2)-xd(1));
    elseif strcmp(action,'Rightsm')
      shift = (xd(2)-xd(1))/10;
    else
      shift = -(xd(2)-xd(1))/10;
    end
    xlim = get(gca,'XLim');
    if (shift > 0)
      shift = min(shift,xlim(2)-xd(2));
    else
      shift = max(shift,xlim(1) - xd(1));
    end
    set(hselrect,'XData',xd+shift);
    sliderwindow('PlotTop',hfig);
   case 'Key',
    c = get(hfig,'CurrentCharacter');
    if (isempty(c) | isempty(double(c)))
      return;    % User must have hit shift
    end
    if (double(c) == 29)
      sliderwindow('Right',hfig);
    elseif (double(c) == 28)
      sliderwindow('Left',hfig);
    elseif strcmp('>',c)
      sliderwindow('Rightsm',hfig);
    elseif strcmp('<',c)
      sliderwindow('Leftsm',hfig);
    else
      return;
    end
    
   case 'ZoomIn'
    hselrect = findobj(hfig,'Tag','HSelRect');
    xd = get(hselrect,'XData');
    yd = get(hselrect,'YData');
    xlim = [min(xd) max(xd)];
    ylim = [min(yd) max(yd)];
    hselax = get(hselrect,'Parent');
    set(hselax,'XLim',xlim,'YLim',ylim);
   case 'ZoomOut'
    xlim = getappdata(hfig,'FullXLim');
    ylim = getappdata(hfig,'FullYLim');
    hselax = findobj(hfig,'Type','axes');
    set(hselax,'XLim',xlim,'YLim',ylim);
   case 'SelectAll',
    hselrect = findobj(hfig,'Tag','HSelRect');
    xlim = get(gca,'XLim');
    xd = [xlim(1) xlim(2) xlim(2) xlim(1) xlim(1)];
    set(hselrect,'XData',xd);
    ylim = get(gca,'YLim');
    yd = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)];
    set(hselrect,'YData',yd);
    sliderwindow('PlotTop',hfig);
   case 'New'
    haxupd = getappdata(hfig,'UpdAx');
    axupdtfcn = getappdata(hfig,'axupdtfcn');
    position = get(hfig,'Position');
    if iscell(axupdtfcn)
      sliderwindow(haxupd,struct('position',position,...
        'axisupdatefcn',{axupdtfcn}));
    else
      sliderwindow(haxupd,struct('position',position,...
        'axisupdatefcn',axupdtfcn));
    end
   case 'RestoreAndKill',
    sliderwindow('ZoomOut',hfig);
    sliderwindow('SelectAll',hfig);
    delete(hfig);
   otherwise,
    action
    error('callback action not recognized')
  end
else
  action
  error('First argument not recognized');
end

% A helper function to draw the paging icons
function icons = swicons
  iconsize = [17 17];
  bkgrnd = get(0,'DefaultUicontrolBackgroundColor');
  blank = zeros(iconsize) + bkgrnd(1);
  
  right = blank;
  index = mask(iconsize,1,1);
  right(index) = 0;
  right = repmat(right,[1 1 3]);

  rightsm = blank;
  index = mask(iconsize,0.5,1);
  rightsm(index) = 0;
  rightsm = repmat(rightsm,[1 1 3]);

  left = blank;
  index = mask(iconsize,1,-1);
  left(index) = 0;
  left = repmat(left,[1 1 3]);
  
  leftsm = blank;
  index = mask(iconsize,0.5,-1);
  leftsm(index) = 0;
  leftsm = repmat(leftsm,[1 1 3]);
  
  icons.right = right;
  icons.rightsm = rightsm;
  icons.left = left;
  icons.leftsm = leftsm;
  
% This function takes parameters from swicons
% to do the actual icon "drawing"
function index = mask(iconsize,width,direction)
  [row,col] = find(ones(iconsize));
  % Make them go from 0 to 1
  row = (row-1)/(iconsize(1)-1);
  col = (col-1)/(iconsize(2)-1);
  if (direction > 0)
    index = find(2*abs(row-0.5) <= (width-col)/width);
  else
    index = find(2*abs(row-0.5) + (1 - width)/width <= col/width);
  end
  return;
