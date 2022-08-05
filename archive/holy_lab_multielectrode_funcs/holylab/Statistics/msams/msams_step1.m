function [ynew,neighborIndex,n,msd] = msams_step1(x,y,options)
% MSAMS_STEP1: move one landmark by minimally-significant adaptive mean shift
% This function takes a single step of adaptive mean shift, for one
% landmark. The step size is determined by the criterion of minimal
% significance, meaning that the step should only barely exceed the
% statistical uncertainty in the center of mass of the landmark's
% neighbors.
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
%   factor (default 3): the number of times that the step size needs to
%     exceed the standard error by in order to be considered
%     statistically significant;
%   min_to_check (default 3): the minimum number of neighbors to use when
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
%   metric_weights: If supplied, then each coordinate in the distance
%     function will be weighted by the weight matrix:
%        d2(x_i,y) = sum_j w_{ij} (x_{ij} - y_j)^2
%   any_coordinate (default false): If true, the MSAMS criterion is
%     applied to each coordinate individually (euclidean metric), and the
%     neighborhood is set as long as any one passes.  If false, we use
%     the total displacement (summing across all coordinates).
%   metric (NOT IMPLEMENTED defaults to euclidean): if you specify a function handle for
%     this field, call it @sqrdistfun, then this is used to compute the
%     pairwise distances between the landmark y and all data points x.
%     The syntax of this function has the following requirements:
%         [sd,extras] = sqrdistfun(y,x)
%     returns the square distances (in sd) and any "extra" information
%     needed for the next computation;
%         centers = sqrdistfun(y,x,extras)
%     returns a d-by-N matrix containing the centroids of regions that
%     progressively contain more neighbors.
% Also,
%   neighborIndex lists the data points in increasing distance from y; if
%     lminfo is supplied this will only include the contributing points
%     rather than all the points.
%   n is a scalar, the number of data points that contributed to the mean
%     shift step.
%   msd is the vector of mean squared displacements (one for each
%     coordinate) of points that contribute to the new mean position. This
%     can be useful in variable metric scenarios, where the scaling of each
%     coordinate is adjusted until all coordinates have unity msd.
%
% See MSAMS_STEP1_VM for a function also supports variable-metric MSAMS
% (VMAMS), in cases where one may not have an a priori metric. One can scale
% the coordinates through an optional input. This older function doesn't
% have the speed optimizations that this one has.
%
% See also: MSAMS, VMAMS_STEP1, MSAMS_STEP1_VM.
  
% Copyright 2006-2009 by Timothy E. Holy

if (nargin < 3)
  options = struct;
end
options = default(options,'factor',3,'min_to_check',3,'backtrack',true,'any_coordinate',false,'plot',false);

[d,N] = size(x);
if isfield(options,'metric_weights')
  weights = options.metric_weights;
  if any(size(weights) ~= size(x))
    error('metric weights do not agree with the dimensionality of x');
  end
else
  weights = ones(size(x));
end
dx = x - repmat(y,1,N);
% Find the neighbors of y
sd = sum(weights .* dx.^2,1);
[ssd,neighborIndex] = sort(sd);
% Find the displacements
dxs = dx(:,neighborIndex);
dxsum = cumsum(weights .* dxs,2);
wsum = cumsum(weights,2);
dx2sum = cumsum(weights.^2 .* dxs.^2,2);

% Implement the statistical criterion
step = zeros(size(dxsum));
uncert = zeros(size(dxsum));
isnz = wsum > 0;
step(isnz) = dxsum(isnz) .^2 ./ wsum(isnz).^2;
uncert(isnz) = dx2sum(isnz) ./ wsum(isnz).^2;
if ~options.any_coordinate
  step = sum(step,1);
  uncert = sum(uncert,1);
end
exceeds_stderr_by_factor = step > options.factor^2 * uncert;
is_significant = any(exceeds_stderr_by_factor,1);
% Don't be fooled by small numbers of points
is_significant(1:options.min_to_check) = false;
% Find the place where we first cross the threshold of statistical
% significance
n = find(is_significant,1,'first');
if isempty(n)
  n = N;
end
% Now do backtracking, to get back to the last point for which the
% minimally-significant criterion with alpha = 1 is true.
if options.backtrack
  exceeds_stderr = step > uncert; % note lack of factor
  is_minimal = any(exceeds_stderr,1);
  n = find(~is_minimal(1:n),1,'last');
  if isempty(n)
    n = 2;
  else
    n = min(n+1,N);  % get the last one where step exceeds uncertainty
  end
end
% Report the new position & mean square distance
%ynew = y;
msd = zeros(size(y));
isnz = isnz(:,n);
%ynew(isnz) = ynew(isnz) + dxsum(isnz,n)./wsum(isnz,n);
msd(isnz) = msd(isnz) + dx2sum(isnz,n)./wsum(isnz,n).^2;
ynew = mean(x(:,neighborIndex(1:n)),2);

if options.plot && d == 2
  hold off
  plot(x(1,:),x(2,:),'k.')
  hold on
  plot(x(1,neighborIndex(1:n)),x(2,neighborIndex(1:n)),'r.')
  plot(y(1),y(2),'bx')
  plot([y(1) ynew(1)],[y(2) ynew(2)],'b')
  axis equal
end
