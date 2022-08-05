function [pos,label] = tag_points_fetch(hax)
% TAG_POINTS_FETCH: retrieve the tagged points from an axis
% Syntax:
%   [pos,label] = tag_points_fetch
%   [pos,label] = tag_points_fetch(hax)
% These fetch the position (pos, a 2-by-ntags matrix) and label (a
%   vector, one entry per tag) from the axis (default gca).
%
% See also: TAG_POINTS.
  
% Copyright 2007 by Timothy E. Holy
  
  if (nargin < 1)
    hax = gca;
  end
  tags = getappdata(hax,'point_tags');
  if ~isempty(tags)
    pos = tags.pos;
    label = tags.label;
  else
    pos = zeros(2,0);
    label = [];
  end
  
