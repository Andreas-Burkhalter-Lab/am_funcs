function split = SplitAxesEvenly(n,gapfrac)
% SPLIT_AXES_EVENLY: generate split points for equally-spaced axes
% Syntax:
%   split = SplitAxesEvenly(n,gapfrac)
% where
%   n is the number of axes along the given dimension (horizontal or
%     vertical)
%   gapfrac is the fraction of the axis size used for the gap between axes
%     (i.e., the width of the gap is gapfrac*width(axis))
% and
%   split is a row vector containing split points that can be passed to
%     SplitHoriz, SplitVert, or SplitGrid.
% 
% See also: SPLITHORIZ, SPLITVERT, SPLITGRID.

% Copyright 2007 by Timothy E. Holy

  if (gapfrac == 0)
    gapfrac = 1e-3;
  end
  width_axes = 1/(n + (n-1)*gapfrac);
  width_gap = gapfrac*width_axes;
  coef_axes = [1:n-1; 1:n-1]; coef_axes = coef_axes(:);
  coef_gap = [0:n-2; 1:n-1]; coef_gap = coef_gap(:);
  split = coef_axes*width_axes + coef_gap*width_gap;
  split = split';

    