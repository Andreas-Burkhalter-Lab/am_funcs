function [pos,linehandle,vertexIndex] = findpoint(cp,hax)
% FINDPOINT: find nearby point on a line
% Syntax:
%   pos = findpoint(cp)
%   [pos,linehandle,vertexIndex] = findpoint(cp,hax)
% where
%   cp is a 2-element vector (in data units)
%   hax is an optional axis handle (default gca)
% and
%   pos is the position (a 2-element vector) of the closest vertex of any
%     line in the plot;
%   linehandle is the handle for the line
%   vertexIndex is the index of the data point closest to cp.

% Copyright 2007 by Timothy E. Holy

  if (nargin < 2)
    hax = gca;
  end
  % Find the point on a line that's nearest to the point cp
  hline = findobj(hax,'type','line');  % get all lines
  min_dist = Inf;
  % Find the closest line, and the closest point on that line
  % Note this only uses the vertices, not the closest-approach
  % Note also that this is done in terms of distance on the screen,
  % rather than distance as defined by the data values.
  units = get(hax,'Units');
  set(hax,'Units','pixels');
  axpos = get(hax,'Position');
  set(hax,'Units',units);
  scalex2 = (axpos(3)/diff(xlim))^2;
  scaley2 = (axpos(4)/diff(ylim))^2;
  for i = 1:length(hline)
    xd = get(hline(i),'XData');
    yd = get(hline(i),'YData');
    dist = scalex2*(xd-cp(1)).^2 + scaley2*(yd-cp(2)).^2;
    [md,index] = min(dist);
    if (md < min_dist)
      min_dist = md;
      pos = [xd(index);yd(index)];
      linehandle = hline(i);
      vertexIndex = index;
    end
  end
