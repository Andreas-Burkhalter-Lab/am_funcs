function tickleneq(hax,len,units)
% TICKLENEQ: make tick marks of equal absolute length across axes
% In matlab there is no "units" property for ticklength, and the length
% is normalized to the size of the axis.  This makes it hard to have
% equal ticklengths in multipanel figures.
%
% Syntax:
%   tickleneq(hax,len,units)
% where
%   hax is a vector of axis handles;
%   len is a scalar specifying tick length in units of ...
%   units is a string giving the units for len.
%
% This works only for 2D plots. I can't figure out how to get it to work
% for 3D plots, as the help says tick lengths are "normalized relative to
% the longest of the visible X-, Y-, or Z-axis annotation lines," and I
% can't figure out an easy way to determine the length of these lines.
  
% Copyright 2004 by Timothy E. Holy
  
  for i = 1:prod(size(hax))
    u = get(hax(i),'Units');
    set(hax(i),'Units',units)
    pos = get(hax(i),'Position');
    maxl = max(pos(3:4));
    tl = get(hax(i),'TickLength');
    set(hax(i),'TickLength',[len/maxl tl(2)],...
               'Units',u);
  end
  