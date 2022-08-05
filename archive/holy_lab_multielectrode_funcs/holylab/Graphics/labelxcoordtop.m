function hghost = labelxcoordtop(hax,xcoord,tickIndices,fmtstr)
% LABELXCOORDTOP: label top of axis with new coordinate
% Syntax:
%   labelxcoordtop(hax,xcoord)
%   labelxcoordtop(hax,xcoord,tickIndices)
%   labelxcoordtop(hax,xcoord,tickIndices,fmtstr)
%   hghost = labelxcoordtop(hax,xcoord,...)
% where
%   hax is the handle of the axis to which you want to add
%     x-coordinate labelling along the top of the axis;
%   xcoord is a vector which gives the values of the x-coordinate;
%   tickIndices (optional) specifies the index of the xcoord values
%     you want to generate ticks & labels for (default: it will
%     choose the same number of ticks as appear along the bottom of
%     the axis, evenly spaced);
%   fmtstr (optional) is a string used to control how the tick
%     labels appear (default: '%.0f');
% and
%   hghost is a handle to the "ghost axis" created to carry the labelling.
%
% See also: GHOSTAXIS.
  
% Copyright 2005 by Timothy E. Holy
  
  haxtmp = gca;  % Save the current axis
  if (nargin < 4)
    fmtstr = '%.0f';
  end
  if (nargin < 3)
    nticks = length(get(hax,'XTick'));
    tickIndices = round(1,length(xcoord),nticks);
  end
  xtl = cell(1,length(tickIndices));
  for i = 1:length(tickIndices)
    xtl{i} = sprintf(fmtstr,xcoord(tickIndices(i)));
  end
  hfig = get(hax,'Parent');
  hghost = axes('Position',get(hax,'Position'),...
                'Visible','on',...
                'Parent',hfig,...
                'XAxisLocation','top',...
                'YTick',[],...
                'Color','none',...
                'XLim',[1 length(xcoord)],...
                'XTick',tickIndices,...
                'XTickLabel',xtl,...
                'TickDir','out');
  if 0 %strcmp(get(hax,'Layer'),'top')
    % To make sure the axis lines show, we have to swap the drawing order
    % of these two axes
    hchildren = get(hfig,'Children');
    indx = findainb([hax hghost],hchildren);
    hchildren(indx) = [hghost hax];
    set(hfig,'Children',hchildren);
    set(hghost,'Layer','top')
  end

  set(hghost,'HandleVisibility','off')
  axes(haxtmp);   % restore the original current axis