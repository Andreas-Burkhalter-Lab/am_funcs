function [haxo,hlineo] = magnifyplot(shape,haxin)
% MAGNIFYPLOT: draw a "magnified viewing area" for a region of a graph
% This function supports two modes of operation:
%   'external' mode is used when the magnified view is in a separate axis
%     already defined by the user;
%   'internal' mode is used when you want to have a blow-up which overlies
%     the current axes (e.g., in a blank space in the current axes)
%
% Syntax:
%   hline = magnifyplot(shape)         % 'external' mode syntax
%   [hax,hline] = magnifyplot(shape)   % 'internal' mode syntax
%   ... = magnifyplot(shape,haxin)
%
% shape is a structure; it has the following fields for either mode:
%   mode: specifies the mode ('external' or 'internal'); if absent it
%     tries to guess on the basis of fields 'haxmag' or 'endbox';
%   startbox: a vector [left bottom width height], all expressed in the
%       units of the data, specifying the region to be magnified;
%   xconnect (optional): if true, draw the top and bottom edges of the
%     start box (default true);
%   yconnect (optional): if true, draw the left and right edges of the
%     start box (default true);
%   boxmapping: a 2-by-n matrix specifying the vertices of the start
%     and end box to be connected by lines.  1 corresponds to lower
%     left, 2 to lower right, 3 to upper right, and 4 to upper
%     left. The vertex on the startbox specified boxmapping(1,i) is
%     connected to the vertex on the endbox specified by
%     boxmapping(2,i);
% In 'external' mode, shape has the additional fields:
%   haxmag: an axis handle in which the magnified view will be placed;
% In 'internal' mode, shape also has these fields instead:
%   endbox: a vector [left bottom width height] (also in data units)
%     giving the region occupied by the magnified view
%
% Note that 'startbox' doesn't actually have to include the whole region
% to be magnified; for example, when magnifying a region along
% the x-axis you might want to draw vertical ticks rather than a box; you
% can accomplish this with shape.xconnect = 0, shape.yconnect = 1, and
% use the startbox to specify the tick separation and height.
%
% haxin is the handle of the axis to which this should be applied
%   (default gca).
%
% On output,
%   hax is an axis handle in which you can draw the exploded view
%     ('internal' mode only)
%   hline is a vector of line handles for the lines which connect the
%     normal and magnified views.

% Copyright 2003 by Timothy E. Holy <holy@wustl.edu>
% Changelog: TEH 2004-09-09 added 'external' mode

  if (nargin < 2)
    haxin = gca;
  end
  if ~isfield(shape,'mode')
    if isfield(shape,'haxmag')
      shape.mode = 'external';
    elseif isfield(shape,'endbox')
      shape.mode = 'internal';
    else
      error(['shape does not contain sufficient parameters, can''t ' ...
             'determine the mode.']);
    end
  end
  if ~isfield(shape,'xconnect')
    shape.xconnect = 1;
  end
  if ~isfield(shape,'yconnect')
    shape.yconnect = 1;
  end
  
  % For external mode, convert coordinates of magnfied axis into data
  % coordinates for the current axis (i.e., create an endbox)
  if strcmp(shape.mode,'external')
    if (get(haxin,'Parent') ~= get(shape.haxmag,'Parent'))
      error('axis and magnified axis must be in the same figure');
    end
    % Get position of axes in figure-normalized units
    axunits = get([haxin shape.haxmag],'Units');
    set([haxin shape.haxmag],'Units','normalized');
    axpos = get([haxin shape.haxmag],'Position');
    set([haxin shape.haxmag],{'Units'},axunits);
    % Convert magnified axes coordinates to data units
    xlim = get(haxin,'XLim');
    ylim = get(haxin,'YLim');
    shape.endbox = ...
        [diff(xlim)*(axpos{2}(1)-axpos{1}(1))/axpos{1}(3) + xlim(1), ...
         diff(ylim)*(axpos{2}(2)-axpos{1}(2))/axpos{1}(4) + ylim(1), ...
         diff(xlim)*axpos{2}(3)/axpos{1}(3), ...
         diff(ylim)*axpos{2}(4)/axpos{1}(4)];
  end
  
  % Convert rectangle 4-vectors to (x,y) coordinates
  xslim = shape.startbox(1)+[0 shape.startbox(3)];
  yslim = shape.startbox(2)+[0 shape.startbox(4)];
  xs = xslim([1 2 2 1]);
  ys = yslim([1 1 2 2]);
  xelim = shape.endbox(1)+[0 shape.endbox(3)];
  yelim = shape.endbox(2)+[0 shape.endbox(4)];
  xe = xelim([1 2 2 1]);
  ye = yelim([1 1 2 2]);
  
  % Draw the lines connecting the axis and magnified axis
  haxlines = ghostaxis(haxin,'fill','fill');
  hline = zeros(0,1);
  if shape.xconnect
    bx = [1 3;2 4];
    hline(end+1:end+2) = line(xs(bx),ys(bx),'Parent',haxlines);
  end
  if shape.yconnect
    bx = [1 2;4 3];
    hline(end+1:end+2) = line(xs(bx),ys(bx),'Parent',haxlines);
  end
  for i = 1:size(shape.boxmapping,2)
    hline(end+1) = ...
        line([xs(shape.boxmapping(1,i)) xe(shape.boxmapping(2,i))],...
             [ys(shape.boxmapping(1,i)) ye(shape.boxmapping(2,i))],...
             'Parent',haxlines);
  end
  set(hline,'Color',[1 0 0],'Visible','on');

  % For internal mode, set up the magnified axis.
  if strcmp(shape.mode,'internal')
    haxo = ghostaxis(haxin,xelim,yelim);
    set(haxo,'XLimMode','auto','YLimMode','auto');
    hlineo = hline;
    axes(haxo);  % Set focus to exploded-view axis
  else
    haxo = hline;
  end
