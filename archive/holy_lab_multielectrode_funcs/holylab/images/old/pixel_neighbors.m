function [memoffset,coordoffset,dist,n_values,sz_spatial] = pixel_neighbors(sz,maxdist,P)
% PIXEL_NEIGHBORS: precompute indexing for neighboring pixels in an image
%
% Note: GRID_NEIGHBORS is recommended instead, this function is deprecated.
%
% This function is a useful utility for any function that processes local
% blocks of pixels. It returns indexing information for pixels in groups,
% where each group corresponds to a set of pixels at the same distance
% from a "base" pixel.
%
% Syntax:
%   [memoffset,coordoffset] = pixel_neighbors(sz,maxdist,P)
%   [memoffset,coordoffset,dist] = pixel_neighbors(sz,maxdist,P)
% where
%   imsz is the size of the underlying image, which can be scalar-valued,
%     RGB, or (for example) a temporal sequence of images; in the latter two
%     cases, the first coordinate of im corresponds to the "value".  For
%     example, a three-dimensional movie-stack would have size [n_t n_x n_y
%     n_z]. An RGB image would have size [3 h w]. In fact, it is possible
%     to have more than one value dimension; the number of spatial
%     dimensions is determined by P (see below).
%     Note that having the color-coordinate first is not the usual MATLAB
%     convention, but improves cache performance. You can use permute to
%     change the order of coordinates.
%   maxdist is the maximum distance, in physical units, which can be
%     considered to be part of a given pixel's neighborhood
%   P is either a vector, in which case it specifies the pixel spacing
%     along each coordinate, or a matrix, with each column of P specifying
%     the n_dimensional displacement (in physical units) between adjacent
%     pixels along the given coordinate. Only spatial coordinates are
%     supplied here; if the image is multiply-valued, you ignore the first
%     coordinate.
% and
%   memoffset is a cell array, each element containing the index offset
%     needed to address all pixels at a given distance from a "base pixel"
%   coordoffset is a cell array of the same size as memoffset, each
%     containing an n_spatial_dims-by-N matrix listing the coordinate
%     offsets corresponding to the memory offsets in memoffset
%   dist is a vector of distances to each group
%
% See also: grid_neighbors, permute.

% Copyright 2010 by Timothy E. Holy

% note: P must be supplied, because it help assess whether the image is
% scalar-valued or multiply-valued. So don't try to create a default
% behavior for it.

  %% Parse the inputs
  if isvector(P)
    P = diag(P);
  end
  n_dims = length(sz);
  n_spatial_dims = size(P,1);
  n_value_dims = n_dims-n_spatial_dims;
  if (n_value_dims > 0)
    n_values = prod(sz(1:n_value_dims));
  else
    n_values = 1;
  end
  sz_spatial = sz(end-n_spatial_dims+1:end);
  
  %% Generate all possible in-range neighbors (in the + quadrant/octant/etc)
  p = diag(P,0);
  n = ceil(maxdist ./ p);
  xc = cell(1,n_spatial_dims);
  for dimIndex = 1:n_spatial_dims
    xc{dimIndex} = 0:n(dimIndex);
  end
  Xc = cell(1,n_spatial_dims);
  [Xc{:}] = ndgrid(xc{:});
  for dimIndex = 1:n_spatial_dims
    Xc{dimIndex} = Xc{dimIndex}(:);
  end
  X = cat(2,Xc{:})';
  
  %% Keep those neighbors that are within range of the base point
  % Compute their distances from the base point in physical units
  PX = P*X;
  d2 = sum(PX.^2,1)/maxdist^2;  % normalized square distance
  % Keep only those within range
  keepFlag = d2 <= 1;
  X = X(:,keepFlag);
  d2 = d2(keepFlag);
  
  %% Collect into groups of equal distance
  % make sure grouping is not messed up by roundoff errors
  fac = 100*eps(class(d2));
  d2 = fac*round(d2/fac);
  [dist,~,lbl] = unique(d2);
  dist = sqrt(dist)*maxdist;
  pixel_groups = agglabel(lbl);
  n_groups = length(pixel_groups);
  
  %% Add remaining quadrants/octants/etc, calculate memory offsets
  pixel_skip = [1 cumprod(sz_spatial)]*n_values;
  memoffset = cell(1,n_groups);
  coordoffset = cell(1,n_groups);
  for groupIndex = 1:n_groups
    coordoffset{groupIndex} = zeros(n_spatial_dims,0);
  end
  for i = 0:2^n_spatial_dims-1
    bits = bitget(i,1:n_spatial_dims);
    sgn = (-1).^bits;
    this_pixel_skip = sgn .* pixel_skip(1:end-1);
    this_memoffset = sum(diag(this_pixel_skip)*X,1);
    this_coordoffset = diag(sgn)*X;
    for groupIndex = 1:n_groups
      memoffset{groupIndex} = [memoffset{groupIndex} this_memoffset(pixel_groups{groupIndex})];
      coordoffset{groupIndex} = [coordoffset{groupIndex} this_coordoffset(:,pixel_groups{groupIndex})];
    end
  end
  % Keep only non-redundant entries
  for groupIndex = 1:n_groups
    [memoffset{groupIndex},keepIndex] = unique(memoffset{groupIndex});
    coordoffset{groupIndex} = coordoffset{groupIndex}(:,keepIndex);
  end
