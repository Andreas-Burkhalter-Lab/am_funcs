function [region,regionIndex] = split_into_contiguous_regions(index,irange)
% SPLIT_INTO_CONTIGUOUS_REGIONS: snippet so that overlapping regions are combined
%
% Given a set of (one-dimensional) indices, break the set of all indices up
% into regions that are buffered by a given width around each index. When
% two (or more) such regions overlap, combine them into one.
%
% Syntax:
%   [region,regionIndex] = split_into_contiguous_regions(index,width)
%   [region,regionIndex] = split_into_contiguous_regions(index,irange)
% where
%   index is a vector of integers
%   width is a scalar specifying the width on either side of index(i) to
%     include in the region (so a single index would have a region of size
%     2*width+1)
% OR
%   irange is a 2-vector, [edge_left edge_right], specifying an asymmetric
%     region (so [-5 10] would mean 5 points before and 10 after)
% and
%   region is a 2-by-n_regions matrix containing the first and last points
%     within each region
%   regionIndex is a 1-by-length(index) array giving the region number
%     associated with each element of index.

  n_points = length(index);
  [sindex,sortOrder] = sort(index);
  if isscalar(irange)
    irange = [-1 1]*abs(irange);
  end
  width = diff(irange);
  breaks = [0 find(diff(sindex) > width) n_points];
  n_intervals = length(breaks)-1;
  region = zeros(2,n_intervals);
  for i = 1:n_intervals
    region(:,i) = [sindex(breaks(i)+1) + irange(1); sindex(breaks(i+1)) + irange(2)];
  end
  breakflag(breaks+1) = 1;
  breakflagcum = cumsum(breakflag); breakflagcum = breakflagcum(1:end-1);
  [tmp,regionIndex] = sort(sortOrder);
  regionIndex = breakflagcum(regionIndex);
