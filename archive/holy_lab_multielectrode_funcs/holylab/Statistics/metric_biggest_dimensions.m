function [sd,weights] = metric_biggest_dimensions(y,x,nd)
% METRIC_BIGGEST_DIMENSIONS: compare points based on largest coordinates
% This metric is useful for "sparse" problems, in which you may have many
% coordinates but want to pay attention only to the largest ones (to,
% e.g., avoid adding noise).  This metric computes the distance between
% any two points by finding their Nd mutually-largest dimensions (in an
% absolute-value sense) and then computing the Euclidean distance between
% them in terms of those coordinates.
%
% Syntax:
%   [sd,weights] = metric_biggest_dimensions(y,x,nd)
% where
%   y is a d-by-1 vector
%   x is a d-by-N matrix
%   nd is the # of dimensions to use
% and
%   sd is a 1-by-N vector of square distances
%   weights is a d-by-N matrix of zeros or ones, each indicating whether
%     a particular coordinate was used or not.
%
% Copyright 2009 by Timothy E. Holy
  
  [d,N] = size(x);
  if (size(y,1) ~= d)
    error('Mismatch in dimensionality');
  end
  if (size(y,2) ~= 1)
    error('The first input must consist of a single point');
  end
  if (nd > d)
    error(['The number of "biggest" dimensions must be less than or equal' ...
	   ' to the number of dimensions']);
  end
  
  yrep = repmat(y,1,N);
  % Compute the weights
  mx = max(abs(yrep),abs(x));
  [smx,sortOrder] = sort(mx,1,'descend');
  sortOrder = sortOrder(1:nd,:);
  ptIndex = repmat(1:N,nd,1);
  indx = sub2ind([d N],sortOrder(:),ptIndex(:));
  weights = zeros(d,N);
  weights(indx) = 1;
  % Compute the square distances
  dy = yrep - x;
  sd = sum(weights .* dy.^2,1);
end
