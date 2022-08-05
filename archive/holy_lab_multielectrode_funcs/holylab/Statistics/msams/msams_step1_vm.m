function [ynew,neighborIndex,n,msd] = msams_step1(x,y,options)
% MSAMS_STEP1: move one landmark by minimally-significant adaptive mean shift
% This function takes a single step of adaptive mean shift, for one
% landmark. The step size is determined by the criterion of minimal
% significance, meaning that the step should only barely exceed the
% statistical uncertainty in the center of mass of the landmark's
% neighbors.
%
% This function also supports variable-metric MSAMS (VMAMS), in cases
% where one may not have an a priori metric. One can scale the coordinates
% through an optional input.
%
% Syntax:
%   ynew = msams_step1(x,y)
% where
%   x is a d-by-N matrix of data points in d-dimensional space;
%   y is a d-by-1 vector, the position of the landmark;
% and
%   ynew is the new position of the landmark, after a single mean shift
%     step.
%
% The syntax
%   [ynew,neighborIndex,n,msd] = msams_step1(x,y,options)
% provides for much finer-grained control. "options" is a structure which
% may have the following fields:
%   scale, a d-by-1 vector containing the factor by which each coordinate
%     should be scaled before performing MSAMS (note the results are scaled
%     back again before returning). Default scale is all 1s. This is used
%     for variable metric MSAMS.
%   factor (default 3): the number of times that the step size needs to
%     exceed the standard error by in order to be considered
%     statistically significant;
%   min_to_check (default 2): the minimum number of neighbors to use when
%     initially considering whether the step is statistically
%     significant. (With backtracking, see below, it's possible for the
%     final number of neighbors to be smaller than this, but only if
%     statistical significance is first triggered beyond this minimum
%     number.)
%   backtrack (default true): once statistical significance with "factor"
%     is triggered, examine collections of closer points to see where the
%     step first exceeds the uncertainty, by a factor of 1. This is useful
%     in "purifying" the collections of contributing points, to reduce the
%     likelihood that points from other clusters will contribute.
% Also,
%   neighborIndex lists the data points in increasing distance from y;
%   n is a scalar, the number of data points that contributed to the mean
%     shift step;
%   msd is the vector of mean squared displacements (one for each
%     coordinate) of points that contribute to the new mean position. This
%     can be useful in variable metric scenarios, where the scaling of each
%     coordinate is adjusted until all coordinates have unity msd.
%
% See also: MSAMS, VMAMS_STEP1.
  
% Copyright 2006 by Timothy E. Holy

if (nargin < 3)
  options = struct;
end
if ~isfield(options,'factor')
  options.factor = 3;
end
if ~isfield(options,'min_to_check')
  %options.min_to_check = 10;
  options.min_to_check = 2;
end
if ~isfield(options,'backtrack')
  options.backtrack = true;
end
% if ~isfield(options,'plot')
%   options.plot = false;
% end
do_scaling = isfield(options,'scale');
if ~isfield(options,'anisotropy_protection')
  options.anisotropy_protection = do_scaling;
end


[d,N] = size(x);
if do_scaling
  % Rescale the coordinates, to allow for variable metric MSAMS
  x = x .* repmat(options.scale(:),1,N);
  y = y .* options.scale(:);
end
% Find the neighbors of y
dx = x - repmat(y,1,N);
sd = sum(dx.^2,1);
[ssd,neighborIndex] = sort(sd);
% Calculate the mean shift step and its uncertainty, separately in each
% coordinate
dxs = dx(:,neighborIndex);
dxsum = cumsum(dxs,2);   % Sum over points; except for /n, this is the step
dx2sum = cumsum(dxs.^2,2); % Sum of squares
% Require that the criterion for statistical significance of the step
% be satisfied for at least one coordinate
exceeds_stderr_by_factor = dxsum.^2 > options.factor^2 * dx2sum;
if options.anisotropy_protection
  coordOK = (dx2sum == repmat(max(dx2sum,[],1),[d 1]));
else
  coordOK = true(d,N);
end
is_significant = any(exceeds_stderr_by_factor & coordOK,1);
% Don't be fooled by small numbers of points
is_significant(1:options.min_to_check) = false;
% Find the place where we first cross the threshold of statistical
% significance
n = find(is_significant,1,'first');
if isempty(n)
  n = N;
end
% else
%   n = n-1;  % Get the last non-significant one
% end
% Now do backtracking, to get back to the last point for which the
% minimally-significant criterion with alpha = 1 is true.
if options.backtrack
  exceeds_stderr = dxsum.^2 > dx2sum; % note lack of factor
  is_minimal = any(exceeds_stderr,1); % deliberately don't require coordOK
  n = find(~is_minimal(1:n),1,'last');
  if isempty(n)
    n = 2;
  else
    n = min(n+1,N);  % get the last one where step exceeds uncertainty
  end
end
% Report the new position
ynew = y + dxsum(:,n)/n;
if do_scaling
  ynew = ynew ./ options.scale(:);
end
% Report other variables of interest
%indexContrib = neighborIndex(1:n);
msd = dx2sum(:,n)/(n-1);  % Note we deliberately do _not_ undo scaling
%msd = dx2sum./repmat(1:N,[d 1]);  % Note we deliberately do _not_ undo scaling

% if options.plot
%   hold off
%   plot(x(1,:),x(2,:),'.')
%   hold on
%   plot(x(1,sort_order(1:n)),x(2,sort_order(1:n)),'r.')
%   axis equal
% end