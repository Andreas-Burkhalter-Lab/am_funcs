function [slope,offset,r] = linregress(x,y)
% LINREGRESS: perform linear regression
% Fit a set of observations (vectors x and y) with a model
%    y = slope*x + offset
%
% Syntax:
%   [slope,offset,r] = linregress(x,y)
% where
%   x and y are vectors containing the two sets of observations;
% and
%   slope is the best-fit slope;
%   offset is the best-fit offset;
%   r is the correlation coefficient (spans [-1 1]).
%
% Alternatively, x can be a vector and Y a matrix, or X and Y are both
% matrices.  In either case, regression is performed on each column of Y.
% Slope, offset, and r are vectors, with one entry for each column of Y.
  
% Copyright 2004 by Timothy E. Holy
  
  if isvector(x)
    x = x(:);
  end
  if isvector(y)
    y = y(:);
  end
  N = size(x,1);
  if (N ~= size(y,1))
    error('X and Y are mismatched');
  end
  if (N == 1)
      error('Need at least two points to do regression');
  end
  if (size(x,2) == 1)
    x = repmat(x,1,size(y,2));
  end
  xdiff = x - repmat(mean(x),N,1);
  slope = sum(xdiff .* y)./sum(xdiff.^2);
  offset = (sum(y) - sum(x).*slope)/N;
  if (nargout > 2)
    ydiff = y - repmat(mean(y),N,1);
    r = nan(size(slope));
    sumsq = sum(xdiff.^2).*sum(ydiff.^2);
    dotprod = sum(xdiff.*ydiff);
    nzindx = find(sumsq);
    r(nzindx) = dotprod(nzindx)./sqrt(sumsq(nzindx));
  end
