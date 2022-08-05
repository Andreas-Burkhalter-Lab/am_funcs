function [x,y] = norm2data(xn,yn,hax)
% NORM2DATA: convert normalized coordinates to data coordinates
%
% Syntax:
%   [xd,yd] = norm2data(xn,yn,hax)
% where
%   xn,yn are the coordinates in figure-normalized units;
%   hax is the axis handle (default gca).
%   xd,yd are the coordinates in data units;
%
% See also: DATA2NORM.

  if (nargin < 3)
    hax = gca;
  end
  % Get position of axis in figure-normalized units
  axunits = get(hax,'Units');
  set(hax,'Units','normalized');
  axpos = get(hax,'Position');
  set(hax,'Units',axunits);
  % Get data limits
  xlim = get(hax,'XLim');
  if strcmp(get(hax,'XDir'),'reverse')
    xlim = xlim([2 1]);
  end
  ylim = get(hax,'YLim');
  if strcmp(get(hax,'YDir'),'reverse')
    ylim = ylim([2 1]);
  end
  
  % Make this also work on log scale
  if(strcmp(get(hax, 'xscale'), 'log'))
     xlim=log10(xlim);
     xn=log10(xn);
  end
  if(strcmp(get(hax, 'yscale'), 'log'))
     ylim=log10(ylim);
     yn=log10(yn);
  end
  
  % Do the conversion
  x = diff(xlim)*(xn - axpos(1))/axpos(3) + xlim(1);
  y = diff(ylim)*(yn - axpos(2))/axpos(4) + ylim(1);

  