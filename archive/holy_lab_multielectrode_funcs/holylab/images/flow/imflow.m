function map = imflow(im)
% IMFLOW: associate each pixel in an image with its brightest neighbor
% This function creates a "map" which describes how each pixel should "flow
% uphill" to a pixel of higher brightness. The idea is that this can then
% define local "peaks" in an image.
%
% Syntax:
%   map = imflow(im)
% where
%   im is an image or stack
% and
%   map is an array of the same size as im, holding for each pixel the
%   integer index (see IND2SUB) of the most-uphill pixel.
%
% Given map, you can then find each pixel's local peak by iterating map:
%    map = map(map)
% until it no longer changes.
%
% Note: if using the MEX version of this function (which is much faster),
% all edge points will have map = 1. If using the M-file version, all the
% edge pixels point to themselves. You can modify the edges using
% KILLEDGES.
%
% See also: KILLEDGES.

% Copyright 2009 by Timothy E. Holy

  if (exist('imflow_mex') == 3) %#ok<EXIST>
    % Use the much-faster MEX version
    map = imflow_mex(im);
    map = killedges(map,1); % fix the edges so map iteration doesn't cause error
    return
  end

  % Get dimensionality and size information
  map = zeros(size(im));
  n_pixels = numel(im);
  im = squeeze(im);  % get rid of unity dimensions
  sz = size(im);
  szfix = sz(sz > 1);  % properly handle "1-dimensional" images
  n_dims = length(szfix);
  coords = cell(1,n_dims);
  cumskip = cumprod([1 szfix]); % the integer-offset for each dimension
  
  map(1:n_pixels) = 1:n_pixels; % starting, each pixel maps to itself

  % Generate the neighbor offset vector (so we can quickly index all
  % neighbors).  The result is the vector nbr_offset. For a pixel indexed
  % by i, the indices of all its nearest-neighbors is i + nbr_offset.
  n_nbrs = 3^n_dims;
  nbr_offset = zeros(1,n_nbrs);
  threes = repmat(3,[1 n_dims]);  % a vector of 3s
  coef_lookup = [-1 0 1];   % to convert [1 2 3] -> [-1 0 1] for the neighbors
  for nbrIndex = 1:n_nbrs
    [coords{:}] = ind2sub(threes,nbrIndex); % convert neighbor index to coordinates
    c = coef_lookup(cat(2,coords{:}));      % convert coordinate to coordinate-offset
    offset = 0;    % this will hold the total integer offset given the size of the input image
    for dimIndex = 1:n_dims
      offset = offset + c(dimIndex)*cumskip(dimIndex);
    end
    nbr_offset(nbrIndex) = offset;
  end
  
  c = ones(1,n_dims);
  nmod = 0;
  if (n_pixels > 1e5)
    nmod = round(n_pixels/20);
  end  
  for pixelIndex = 1:n_pixels
    if (nmod && mod(pixelIndex,nmod) == 0)
      fprintf('...%d%% done',round(100*pixelIndex/n_pixels));
    end
%     [coords{:}] = ind2sub(sz,pixelIndex); % convert to coords
%     c = cat(2,coords{:});
    if (all(c > 1) && all(c < szfix))
      % if it is an interior point (i.e., not on edge so has all neighbors)
      imnbr = im(pixelIndex + nbr_offset);
      [mx,mxIndex] = max(imnbr);
      map(pixelIndex) = pixelIndex + nbr_offset(mxIndex);
    end
    dimIndex = 1;
    c(dimIndex) = c(dimIndex)+1;
    while (c(dimIndex) > szfix(dimIndex) && dimIndex < n_dims)
      c(dimIndex) = 1;
      dimIndex = dimIndex+1;
      c(dimIndex) = c(dimIndex)+1;
    end
  end
  if (nmod)
    fprintf('\n');
  end
end
