function tags_out = tag_points(action,varargin)
% TAG_POINTS: a utility to tag data points with a label
% Syntax:
%   tag_points on
%   tag_points('on')
% turns on tagging for the current axis. Right click a location on the
% axes, and a context menu will pop up. By default, each click will be
% associated with the closest point on any line in the plot (but this
% behavior can be adjusted). Choose a numeric label for each click. You
% can re-use labels if you want to classify several points the same way.
%
%   tag_points off
%   taginfo = tag_points off
%   taginfo = tag_points('off');
% turns off tagging and deletes all current tags. The second two syntaxes
% pass back the complete tag information, in case you want to resume it
% later.
%
%   tag_points('on',taginfo)
% Allows you to resume tagging after you've turned it off, or add tags to
% an existing taginfo structure.
%
%   tag_points('on',options)
% allows you to control the behavior in greater detail.
% Here are the fields of options:
%   next_label (default 1): the integer label used for the next new class
%     of tagged object
%   on_line (default true): if true, the tag will be associated with a
%     point on a line object that is closest to the clicked point.
%   hit_off (default true): turns off hit testing on all children of the
%     axis
%   pointdir (default 0): specify the direction for drawing labels:
%     1 draws above the point
%     -1 draws below the point
%     0 will draw above if the click is at a positive-y location, and below
%       for a negative-y location.
% These are all fields of the taginfo structure, too, so you can alter
% the behavior when you resume tagging.
%
%   tag_points(action,hax,...)
%   tag_points(action,...,hax)
% sets up tagging in the axis specified by the axis handle hax.
%
% Users who wish to manipulate the available labels may be interested in
% different actions that can return handles to the internal subfunctions of
% tag_points. For example, tag_points('getlabelfunc') returns a handle to
% the function that draws the text label in the axes. See the code for
% further examples.
%
% Example:
% figure; plot(rand(20,2))
% tag_points on
% % Now right click near some points on the lines; you could, for example,
% % label some points on the blue curve "1" and some points on the green
% % curve "2".
% [pos,label] = tag_points_fetch
% tp = tag_points('off')
% % Now look at your plot: the labels are gone! And the context menu is gone
% slidwincmenu(gca) % to put something into a context menu
% tag_points('on',tp)
% % Now note that your old points are redrawn, and the tagging items have
% % been appended to the context menu
%
% See also: TAG_POINTS_FETCH.

% Copyright 2007 by Timothy E. Holy

  if ~ischar(action)
    error('Must indicate action (''on'' or ''off'')');
  end
  hax = gca;
  for i = 1:length(varargin)
    if ishandle(varargin{i})
      hax = varargin{i};
      if ~strcmp(get(hax,'type'),'axes')
        error('Object handle is not an axis!');
      end
      break
    end
  end
  options = struct;
  for i = 1:length(varargin)
    if isstruct(varargin{i})
      options = varargin{i};
      break
    end
  end
  
  switch(action)
    case 'on'
      % Here we have to cope with either brand-new calls, or calls in which
      % taginfo is supplied
      options = default(options,'next_label',1);
      options = default(options,'on_line',true);
      options = default(options,'hit_off',true);
      options = default(options,'pointdir',0);

      tags = options; % taginfo starts out as the options structure
      % If needed, set up empty fields for data
      tags = default(tags,'pos',[]);  % the position of a tag
      tags = default(tags,'label',[]); % the label associated with a tag
      tags = default(tags,'htext',[]); % the handle for the text label

      if isappdata(hax,'point_tags')
        warning('point_tags:overwrite',...
          'The given axis already has tagging on, but the taginfo is being overwritten');
        old_tags = get(hax,'point_tags');
        delete(old_tags.htext);
      end

      cmenu = get(hax,'UIContextMenu');
      if isempty(cmenu)
        cmenu = uicontextmenu('Parent',get_parent_fig(hax));
        set(hax,'UIContextMenu',cmenu);
      end

      % Add the tagging menu items
      add_separator = ~isempty(get(cmenu,'Children'));
      separator_text = {'off','on'};
      separator_text = separator_text{add_separator+1};
      % The "delete tag" option
      uimenu(cmenu,'Label','Del','Callback',@wst_del,...
        'Tag', 'tagcmenudel','Separator',separator_text);
      % Populate with any existing labels
      ulabel = unique(tags.label);
      for i = 1:length(ulabel)
        wst_newlabel(cmenu,ulabel(i));
      end
      % Show the labels
      n_labels = length(tags.label);
      tags.htext = zeros(1,n_labels);
      for i = 1:n_labels
        tags.htext(i) = tp_textlabel(tags.pos(:,i),tags.label(i),tags.pointdir);
      end
      % Add one extra label in case the user wants to define a new class
      wst_newlabel(cmenu,tags.next_label);

      if options.hit_off
        hc = get(hax,'Children');
        set(hc,'HitTest','off');
      end

      setappdata(hax,'point_tags',tags);



    case 'off'
      tags = getappdata(hax,'point_tags');
      if ~isempty(tags)
        delete(tags.htext);
        tags.htext = [];
        rmappdata(hax,'point_tags');
      end
      % Clear out items from the context menu
      cmenu = get(hax,'UIContextMenu');
      hobj = findobj(cmenu,'Tag','tagcmenudel');
      delete(hobj);
      hobj = findobj(cmenu,'Tag','tagcmenulabel');
      delete(hobj);
      % Delete the context menu if it's empty
      hobj = get(cmenu,'Children');
      if isempty(hobj)
        delete(cmenu);
      end
      if (nargout > 0)
        tags_out = tags;
      end
    case 'getlabelfunc'
      tags_out = @tp_textlabel;
    case 'getnewlabelfunc'
      tags_out = @tp_newlabel;
    case 'getorganizelabelsfunc',
      tags_out = @tp_organizelabels;
    otherwise
      error('Action must be ''on'' or ''off''');
  end


function wst_newlabel(cmenu,labelnum)
  uimenu(cmenu,'Label',num2str(labelnum),...
    'Callback',@wst_add,'Tag','tagcmenulabel');

       
function wst_add(obj, event)
  hax = gca;
  cp = get(hax,'CurrentPoint');
  cp = cp(1,1:2);
  cmenu = get(obj,'Parent');
  tags = getappdata(hax,'point_tags');
  label = str2double(get(obj,'Label'));
  % If necessary, add a new menu item for next time
  if (isempty(tags.label) || label == tags.next_label)
    tags.next_label = tags.next_label+1;
    wst_newlabel(cmenu,tags.next_label);
  end
  % Now gravitate towards point, if necessary
  [pos,pointdir] = wst_findpoint(hax,cp,tags);
  % See if that position is already occupied
  occupied = false;
  if ~isempty(tags.pos)
    dx = tags.pos - repmat(pos,1,length(tags.label));
    dist = sum(dx.^2,1);
    if any(dist == 0)
      tags.label(dist == 0) = label;
      set(tags.htext(dist == 0),'String',num2str(label));
      occupied = true;
    end
  end
  if ~occupied
    % This is a new position, add it in
    tags.pos(:,end+1) = pos;
    tags.label(end+1) = label;
    tags.htext(end+1) = tp_textlabel(pos,label,pointdir);
  end
  setappdata(hax,'point_tags',tags);

  
function wst_del(obj, event)
  hax = gca;
  cp = get(hax,'CurrentPoint');
  cp = cp(1,1:2);
  tags = getappdata(hax,'point_tags');
  % Now gravitate towards point, if necessary
  pos = wst_findpoint(hax,cp,tags);
  % See if that position is already occupied
  dx = tags.pos - repmat(pos,1,length(tags.label));
  dist = sum(dx.^2,1);
  delIndex = dist == 0;
  if any(delIndex)
    tags.label = tags.label(~delIndex);
    tags.pos = tags.pos(:,~delIndex);
    delete(tags.htext(delIndex));
    tags.htext = tags.htext(~delIndex);
  end
  setappdata(hax,'point_tags',tags);

  
function [pos,pointdir] = wst_findpoint(hax,cp,tags)
  pos = cp(:);
  pointdir = tags.pointdir; % decide which direction to plot label
  if (pointdir == 0)
    pointdir = sign(cp(2)); % if > 0, seek a max; <0, a min
  end
  if tags.on_line
    pos = findpoint(cp,hax);
  end
  
function htext = tp_textlabel(pos,label,pointdir)
  if (nargin < 3)
    pointdir = 1;
  end
  htext = text(pos(1),pos(2),num2str(label));
  if (pointdir > 0)
    set(htext,'VerticalAlignment','bottom');
  else
    set(htext,'VerticalAlignment','top');
  end
  set(htext,'HitTest','off');

function tp_organizelabels(cmenu)
  % Reorder the children to keep them in numeric order---this can be
  % necessary if outside functions have added an extra label to the context
  % menu, but not in numeric order.
  hc = get(cmenu,'Children');
  hlabels = findobj(hc,'flat','Tag','tagcmenulabel');
  label_strings = get(hlabels,'Label');
  labels = zeros(1,length(hlabels));
  if iscell(label_strings)
    for i = 1:length(label_strings)
      labels(i) = str2num(label_strings{i});
    end
  else
    labels = str2num(label_strings);
  end
  [slabels,sortIndex] = sort(labels);
  labelIndices = findainb(hlabels,hc);
  hc(labelIndices) = hlabels(sortIndex);
  set(cmenu,'Children');
