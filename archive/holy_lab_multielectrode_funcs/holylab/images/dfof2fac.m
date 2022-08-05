function fac = dfof2fac(dfof,clim_dfof)
% dfof2fac: normalize deltaf/f values
% Syntax:
%   fac = dfof2fac(dfof,clim_dfof)
% where
%   dfof is an array containing deltaf/f values for each pixel;
%   clim_dfof is a three-vector, [sat_neg neutral sat_pos], containing the
%     saturating values for negative and positive responses as well as the
%     "neutral" response. For example, [-0.05 0 0.2] would map deltaf/f
%     values of -0.05 or lower to -1, 0 to 0, and 0.2 or higher to 1. You
%     can use [-Inf 0 0.2] to truncate all negative values.
% and
%    fac is an array of the size of dfof, mapped to the range [-1 1].
%
% See also: colorize_pixels.

% Copyright 2012 by Timothy E. Holy

  if length(clim_dfof) ~= 3 || ~issorted(clim_dfof)
    error('dfof must be a sorted 3-vector');
  end
  fac = zeros(size(dfof));
  % Do the first interval
  flag = dfof < clim_dfof(2);
  m = 1/diff(clim_dfof(1:2));
  b = -m*clim_dfof(2);
  fac(flag) = m*dfof(flag)+b;
  % Do the second interval
  flag = ~flag;
  m = 1/diff(clim_dfof(2:3));
  b = -m*clim_dfof(2);
  fac(flag) = m*dfof(flag)+b;
  % Truncate
  fac(fac<-1) = -1;
  fac(fac>1) = 1;
end