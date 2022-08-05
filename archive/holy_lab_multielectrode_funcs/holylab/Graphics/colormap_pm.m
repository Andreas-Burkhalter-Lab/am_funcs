function cm = colormap_pm(im,options)
% COLORMAP_PM: colormaps for indicating positive/negative values
% Syntax:
%  cm = colormap_pm(im)
% Here, im can be an image, or it can be a 2-vector giving min/max values.
%  cm = colormap_pm(im,'symmetrical')
% The first syntax will generate a red (for positive)/green (for
%   negative) colormap, where the red/green are scaled to their maxima, the
%   second so that the intensity is symmetric in value around white
%
% More flexibly,
%   cm = colormap_pm(im,options)
% where options is a structure with the following fields:
%   background: a 3-vector (default [0 0 0]), use [1 1 1] for when neutral
%     is white
%   plus: a 3-vector, default [1 0 0], the color for positive values
%   minus: a 3-vector, default [0 1 0], the color for negative values
%   symmetrical: a boolean, default false

% Copyright 2007-2010 by Timothy E. Holy

  if (nargin > 1)
    if ischar(options)
      flag = options;
      options = struct;
      options.symmetrical = strcmp(flag,'symmetrical');
    end
  end
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'background',[1 1 1],'plus',[1 0 0],'minus',[0 1 0],'symmetrical',false);
 
  if options.symmetrical
    minmax = [-1 1];
  else
    minmax = [min(im(:)) max(im(:))];
    minmax = minmax/max(abs(minmax));
  end
  colval = linspace(minmax(1),minmax(2),256);
  cm = zeros(length(colval),3);
  for i = 1:length(colval)
    f = colval(i);
    if (f <= 0)
      cm(i,:) = -f*options.minus + (1+f)*options.background;
    else
      cm(i,:) = f*options.plus + (1-f)*options.background;
    end
  end
    