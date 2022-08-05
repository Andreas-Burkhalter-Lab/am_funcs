function hax_inset = create_inset(hax,pos)
% CREATE_INSET: make an inset axis
% Syntax:
%   hax_inset = create_inset(hax,pos)
% where
%   hax is the handle of the original axis
%   pos is the [left bottom width height] of the location of the inset,
%     in normalized units relative to the original axis
% and
%   hax_inset is the handle for the inset axis
  
% Copyright 2007 by Timothy E. Holy
  
  units_parent = get(hax,'Units');
  set(hax,'Units','normalized');
  pos_parent = get(hax,'Position');
  pos_inset = [pos_parent(1:2)+pos(1:2).*pos_parent(3:4),...
	       pos_parent(3:4).*pos(3:4)];
  hax_inset = axes('Units','normalized','Position',pos_inset);
  set(hax,'Units',units_parent);
  