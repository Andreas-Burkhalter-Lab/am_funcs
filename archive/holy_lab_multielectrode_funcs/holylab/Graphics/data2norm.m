function [x,y] = data2norm(xd,yd,hax)
% DATA2NORM: convert data coordinates to normalized coordinates
%
% Syntax:
%   [xn,yn] = data2norm(xd,yd,hax)
% where
%   xd,yd are the coordinates in data units;
%   xn,yn are the coordinates in figure-normalized units;
%   hax is the axis handle (default gca).
%
% See also: NORM2DATA.   

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
  
  % make this also works on log scale
  if(strcmp(get(hax, 'xscale'), 'log'))
     xlim=log10(xlim);
     xd=log10(xd);
  end
  if(strcmp(get(hax, 'yscale'), 'log'))
     ylim=log10(ylim);
     yd=log10(yd);
  end
  
  % Do the conversion
  x = axpos(3)*(xd - xlim(1))/diff(xlim) + axpos(1);
  y = axpos(4)*(yd - ylim(1))/diff(ylim) + axpos(2);

  