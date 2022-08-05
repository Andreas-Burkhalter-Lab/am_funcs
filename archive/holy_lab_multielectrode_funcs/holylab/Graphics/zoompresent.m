function [hax_out,hbox] = zoompresent(params,h,rect)
% ZOOMPRESENT: make figure pairs that facilitate zoom animations in presentations
%
% Powerpoint, OpenOffice, and perhaps other presentation software have
% (limited) ability to simulate a "zoom in" or "zoom out" on a portion of a
% figure. This is not a "real zoom" transition as, say, you might perform
% with a video camera; if that's what you require, you may need to prepare
% the transition as a movie clip. However, presentation software can give an
% approximation to this ideal by animating two separate graphics objects,
% corresponding to "zoomed in" and "zoomed out" versions of your figure,
% in a way that visually indicates the relationship between them.
%
% To insure that the resolution is acceptable in both zoomed in & zoomed out
% configurations, it's best to prepare two separate figure files.  (Ideally,
% these figures are of the same physical size so that one obscures the other
% before the animation begins.)  Within Matlab, starting from a "zoomed out"
% figure that the user prepares, this function creates 2 or 3 equal-sized
% figures at the two viewing scales (zoomed in and zoomed out).  Optionally,
% it also draws a box around the zoom region in these figures so that the
% correspondence between them is clearer.
%
% Syntax:
%   [hax = zoompresent(params)
%   [hax,hbox] = zoompresent(params,hfig)
%   [hax,hbox] = zoompresent(params,hax,rect)
% where
%   params is a structure with the following fields:
%     xlim, ylim: specify the x- and/or y-limits of the zoom region. If
%       you want to preserve the aspect ratio of your plot, specify only
%       one of these and then specificy the corresponding xcenter,
%       ycenter as described below.
%     xcenter, ycenter: if you specify only xlim above, then set ycenter
%       to be the center of the zoomed-in region along the yaxis.
%       Alternatively, if you specify only ylim above, set xcenter.
%     draw_box (default true): set to true if you want a box around the
%       zoomed region.  It's a good idea to set 'TickDir' to 'out' in
%       your original plot if you are using this feature, at least unless
%       you choose to animate the version produced by "noaxes" below.
%     box_linewidth (default 1): controls the width of the zoom box
%     box_color (default 'r'): controls the color of the zoom box
%     noaxes (default true): creates a 2nd zoomed-in variant that
%       removes all the axis information (labels, ticks, etc). You might
%       want to apply the "complex" animations to this image, and use the
%       other two as initial/final endpoints.
%   hfig allows you to specify a figure other than the current figure;
%   hax allows you to specify zooming for a particular axis in a figure;
%   rect gives the position in normalized coordinates of the zoomed-in axis
%     in the figure;
% and
%   hax is a vector of handles to the zoomed axes in their separate
%     figures;
%   hbox is a vector of handles to the zoomed-in box and zoomed-out box,
%     respectively. Using these handles, you can manipulate further
%     properties of the boxes.
%
% The resulting figures can then be exported as separate files for
% presentation.  It is best not to resize the figures before exporting them
% (unless you resize them identically); instead, resize your original figure
% before calling this function.
%
% Example:
%   x = linspace(0,5*pi,1000);
%   y = sin(x.^2);
%   plot(x,y)
%   xlabel('Time (s)'); ylabel('Voltage (V)')
%   set(gca,'TickDir','out','Box','off')
%   zoompresent(struct('xlim',[13 15],'ycenter',0.1))
%
% Tips for using these figures in your presentation software: briefly,
% combine a "faded zoom" entrance/exit with a custom motion path, and
% optionally a "fade" on the final image, all triggered on a single
% mouse click. Don't forget to adjust the animation speeds.  For more
% detail, see tutorials available on the web.
% Minimum versions: MSOffice 2003 or later, OpenOffice 2.3 or later.

% Copyright 2008 by Timothy E. Holy
  
  % Parse the arguments, and figure out what mode we're running in
  copy_whole_figure = true;
  if (nargin < 2)
    hfig = gcf;
  else
    if ~ishandle(h)
      error('Second argument must be a handle to a graphics object.');
    end
    type = get(h,'type');
    switch type
      case 'axes'
        hax = h;
        hfig = get(hax,'Parent');
        copy_whole_figure = false;
      case 'figure'
        hfig = h;
      otherwise
        error('Object handle type not recognized');
    end
  end
  
  % Set the default behavior
  params = default(params,'draw_box',true,...
    'box_linewidth',1,...
    'box_color','r',...
    'noaxes',true);
  
  if copy_whole_figure
    hax = childaxis(hfig);
  end
  hax0 = hax;
  xlim = get(hax,'XLim');
  ylim = get(hax,'YLim');
  if ~isfield(params,'xlim')
    % Assume the user intends to preserve the aspect ratio
    % (Among other things, that makes layering in powerpoint more
    % straightforward)
    xwidth = diff(xlim)/diff(ylim) * diff(params.ylim);
    params.xlim = params.xcenter + [-1 1]*xwidth/2;
  elseif ~isfield(params,'ylim')
    ywidth = diff(ylim)/diff(xlim) * diff(params.xlim);
    params.ylim = params.ycenter + [-1 1]*ywidth/2;
  end

  % Devise an overall zoom factor as the geometric mean of the zoom on
  % each axis (which may be uniform if aspect ratio is "locked").  This
  % is used to determine the thickness of zoombox lines
  zoomfactor = sqrt(diff(xlim)/diff(params.xlim) * ...
    diff(ylim)/diff(params.ylim));
  
  % Make copies of the axis/figure, with all the same settings
  % The first will be for the zoomed-in plot
  hzoom = copyobj(hfig,0);
  if params.noaxes
    % The 2nd will be for the axis-free zoomed-in version
    hzoom_noax = copyobj(hfig,0);
  end
  if ~copy_whole_figure
    % We're copying just 1 axis from what is probably a multipanel
    % figure. We'll keep the figures we've already created so that
    % properties of the figure are preserved, but we first have to clear
    % the children from the plot.  Note this assumes there are no hidden
    % handles.
    clf(hzoom)
    copyobj(hax,hzoom);
    if params.noaxes
      clf(hzoom_noax)
      copyobj(hax,hzoom_noax);
    end
  end
  hax_out = correspondingaxis(hzoom,hfig,hax);
  if params.noaxes
    hax_out(end+1) = correspondingaxis(hzoom_noax,hfig,hax);
  end

  
  % The last figure copy is full-size, but needed to avoid corrupting the
  % original if we're drawing a zoom box
  if params.draw_box
    hfull = copyobj(hfig,0);
    hax_out(end+1) = correspondingaxis(hfull,hfig,hax);
  end
  
  % Do the zooming in
  hax = childaxis(hzoom);
  set(hax,'XLim',params.xlim,'YLim',params.ylim);
  if params.noaxes
    hax = childaxis(hzoom_noax);
    set(hax,'XLim',params.xlim,'YLim',params.ylim,'Visible','off');
  end
  
  % Draw the boxes
  if params.draw_box
    hboxtmp = line(params.xlim([1 2 2 1 1]),params.ylim([1 1 2 2 1]),...
      'Parent',hax,'LineWidth',zoomfactor*params.box_linewidth,...
      'Color',params.box_color);
    hax = correspondingaxis(hfull,hfig,hax0);
    hboxtmp(2) = line(params.xlim([1 2 2 1 1]),params.ylim([1 1 2 2 1]),...
      'Parent',hax,'LineWidth',params.box_linewidth,...
      'Color',params.box_color);
    if (nargout > 0)
      hbox = hboxtmp;
    end
  end
  
  if (nargin > 2)
    set(hax_out(1:end-1),'Position',rect);
  end

function hax = childaxis(hfig)
  hax = findobj(hfig,'type','axes');
  if (length(hax) ~= 1)
    error('Axis is ambiguous, use only one (handle-visible) axis in the figure');
  end

function haxout = correspondingaxis(figtarget,figsource,axsource)
  hc = get(figsource,'Children');
  flag = hc == axsource;
  hc = get(figtarget,'Children');
  haxout = hc(flag);