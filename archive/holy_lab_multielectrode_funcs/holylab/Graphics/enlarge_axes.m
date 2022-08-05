function enlarge_axes(hax,factor)
% ENLARGE_AXES: swell or shrink an axis in a figure
%
% The axis expands or contacts around its center point, by the
% magnification factor you supply. Factors bigger than 1 make the axis
% grow, factors smaller than 1 make it shrink.
%
% Syntax:
%   enlarge_axes(hax,factor)

% Copyright 2011 by Timothy E. Holy

  pos = get(hax,'Position');
  wh = pos(3:4);
  new_wh = factor*wh;
  diff_wh = (factor-1)*wh;
  set(hax,'Position',[pos(1:2)-diff_wh/2 new_wh]);
  