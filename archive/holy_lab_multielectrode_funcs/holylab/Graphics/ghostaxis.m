function hgax = ghostaxis(haxin,xlim,ylim)
% GHOSTAXIS: add a layered axis (useful for scalebars, annotation)
% Often annotation needs to be added just outside of the plotting region
% of an axis. However, this material cannot be added to the axis, because
% it will be cut off. This function provides an invisible axis with an
% overlapping coordinate system, which can be made larger or smaller than
% the reference axis. Consequently, annotation can easily be added.
%
% hgax = ghostaxis(haxin,xlim,ylim)
% where
%   haxin is the handle of the reference axis (i.e., your plot)
%   xlim,ylim are the limits of the new ghost axis. These determine the
%     position of the ghost axis, being set so that the coordinates of
%     the reference and ghost axes are consistent. If left empty, these
%     parameters default to the value of the reference axis. If set to
%     'fill', the ghost axis will fill the figure window.
%   hgax is the handle of the ghost axis. The following axis properties
%     are set:
%         'Visible'->'off'
%         'HandleVisibility'->'off'
%         'NextPlot'->'replacechildren'
%     You will need to set 'Visible'->'on' on any object that you parent
%     to the ghost axis.
raxxlim = get(haxin,'XLim');
raxylim = get(haxin,'YLim');
if (nargin < 3)
  ylim = raxylim;
end
if (nargin < 2)
  xlim = raxxlim;
end
axu = get(haxin,'Units');
set(haxin,'Units','normalized');
axpos = get(haxin,'Position');
if ischar(xlim)
  if strcmp(xlim,'fill')
    stretch = diff(raxxlim)/axpos(3);
    xlim = raxxlim(1) - stretch*axpos(1);
    xlim(2) = xlim(1) + diff(raxxlim)/axpos(3);
  else
    error(['String ' xlim ' not recognized']);
  end
end
if ischar(ylim)
  if strcmp(ylim,'fill')
    stretch = diff(raxylim)/axpos(4);
    ylim = raxylim(1) - stretch*axpos(2);
    ylim(2) = ylim(1) + diff(raxylim)/axpos(4);
  else
    error(['String ' ylim ' not recognized']);
  end
end
set(haxin,'Units',axu); % Restore the original units

    
% Calculate the limits of the ghost axis in figure coordinates
gxf = axpos(3)*(xlim - raxxlim(1))/diff(raxxlim) + axpos(1);
gyf = axpos(4)*(ylim - raxylim(1))/diff(raxylim) + axpos(2);
% Create the new axis
hgax = axes('Position', [gxf(1) gyf(1) diff(gxf) diff(gyf)], ...
            'XLim', xlim, 'YLim', ylim, 'Visible', 'off', ...
            'HandleVisibility', 'off', 'NextPlot', 'replacechildren');
