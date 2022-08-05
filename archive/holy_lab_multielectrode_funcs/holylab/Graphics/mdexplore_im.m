function hax = mdexplore_im(im,options)
% MDEXPLORE_IM: explore multidimensional data interactively from image
% This displays a matrix as an image. When the user clicks on a pixel, a
% new plotting function is called to display a more detailed view of the
% data underlying that pixel.
%
% Syntax:
%   hax = mdexplore_im(im,params)
% where
%   im is a m-by-n image
%   params must have the following required fields:
%     plotfunc: a function handle for creating the more detailed
%       view. plotfunc must have the following syntax:
%         plotfunc(col,row)
%       where col is the column number and row is the row number (i.e., x
%       and y) of the clicked pixel in the image. You should make use of
%       function handles in specifying any additional data for your plot
%       function, for example:
%         params.plotfunc = @(col,row) my_plotfunc(col,row,plotting_data)
%       where plotting_data consists of one or more arguments that your
%       plotting function needs to generate its plots. MDEXPLORE_IM_EPHYS
%       provides an example of such a function.
%
%   params may also have the following optional fields:
%     clim: color limits to apply to the image. For more sophisticated
%       handling of color scaling, simply call COLORMAP after you call
%       this function.  Note that the function COLORMAP_PM is especially
%       useful if you want to encode positive/negative values with
%       different colors.
%     hax: the axis to be used for displaying the image
%     xlabel,ylabel: cell arrays of strings to be used for labeling along
%       the x- and y- axes
%
%  On output, hax is the axis handle of the deltarate image plot.
%
% See also: MDEXPLORE, MDEXPLORE_IM_EPHYS, DRGUI, COLORMAP_PM.

% Copyright 2007 by Timothy E. Holy
  
  if (nargin < 2)
    error('Must have 2 input arguments');
  end
  if ~isstruct(options) || ~isfield(options,'plotfunc')
    error('You must supply a plotting function');
  end
  
  if isfield(options,'hax')
    hax = options.hax;
  else
    hfig = figure;
    hax = gca;
    set(hfig,'HandleVisibility','callback');
  end
  
  if isfield(options,'clim')
    him = image(im,'Parent',hax,'CDataMapping','scaled');
    set(hax,'CLim',options.clim);
  else
    him = imagesc(im,'Parent',hax);
  end
  
  [nr,nc] = size(im);
  % Draw the labels
  if isfield(options,'xlabel')
    set(hax,'XTick',1:nc,'XTickLabel',options.xlabel)
  end
  if isfield(options,'ylabel')
    set(hax,'YTick',1:nr,'YTickLabel',options.ylabel)
  end
  set(hax,'TickDir','out')
  
  % Set up the clicking callback
  set(him,'ButtonDownFcn',{@mdim_click,options});
  
function mdim_click(hobject,eventdata,params)
  % Get the current mouse position
  hax = get(hobject,'Parent');
  cp = get(hax,'CurrentPoint');
  cp = round(cp(1,1:2));
  xlim = get(hobject,'XData');
  ylim = get(hobject,'YData');
  if (cp(1) > 0 && cp(1) <= xlim(end) && ...
      cp(2) > 0 && cp(2) <= ylim(end))
    params.plotfunc(cp(1),cp(2));
  end
