function [n,imfilteredOut] = bn_image_filter(imdata,varargin)
% BN_IMAGE_FILTER: balanced neighborhood spatial filtering of images
%
% This function takes an input image (which may be multivalued, i.e., RGB
% or a temporal sequence) and returns the size of the balanced
% neighborhood at each pixel and, optionally, a filtered version of the
% image.
%
% Syntax:
%  [n,imfiltered] = bn_image_filter(im,memoffset,coordoffset)
%  [n,imfiltered] = bn_image_filter(imdata,imfiltered,memoffset,coordoffset)
%  [n,imfiltered,n_triggered,ratio] = bn_image_filter(...,options)
% where
%   imdata is the underlying image, which can be scalar-valued, RGB, or
%     (for example) a temporal sequence of images; in the latter two cases,
%     the first coordinate of im corresponds to the "value".  For example, an
%     RGB image would have size [3 h w], and a three-dimensional
%     movie-stack would have size [n_t n_x n_y n_z].  Note that having the
%     color-coordinate first is not the usual MATLAB convention, but it
%     improves cache performance. You can use permute to change the order
%     of coordinates.
%   memoffset and coordoffset are specified as in the outputs of
%     pixel_neighbors.
%   options has the following fields:
%     z (default 3): the threshold for triggering the balanced
%       neighborhood criterion.
% and
%   n is an array of a size given by the spatial dimensions of the image,
%     giving the number of pixels within the balanced neighborhood of
%     each pixel.
%   imfiltered is a filtered version of the image.
%   n_triggered is the value of n that triggered the balanced
%     neighborhood criterion (0 if the criterion was never triggered)
%   ratio is the average ratio of the displacement to its standard
%     deviation after adding each new pixel group
%
% Note that using the "memoffset" and "coordoffset" inputs (described
% below), you have a great deal of control as to how the filtering is
% performed.  While pixel_neighbors supplies settings for spatial
% filtering, temporal-recency filtering (or spatio-temporal filtering)
% can easily be performed with suitably-specified memoffset &
% coordoffset.
%
% See also: pixel_neighbors, permute.

% Copyright 2010 by Timothy E. Holy

  %% Parse the inputs
  curarg = 1;
  if isnumeric(varargin{1}) && isequal(size(imdata),size(varargin{1}))
    imfilteredIn = varargin{1};
    curarg = curarg+1;
  else
    imfilteredIn = imdata;
  end
  memoffset = varargin{curarg};
  coordoffset = varargin{curarg+1};
  curarg = curarg+2;
  if (curarg > length(varargin))
    options = struct;
  else
    options = varargin{curarg};
  end
  options = default(options,'z',3);
  z2 = options.z^2;

  %% Parse info about image dimensions, pixel neighbors
  sz = size(imdata);
  n_spatial_dims = size(coordoffset{1},1);
  n_value_dims = length(sz) - n_spatial_dims;
  n_value_dims = max(n_value_dims,1);
  sz_spatial = sz(end-n_spatial_dims+1:end);
  n_pix = prod(sz_spatial);
  n_groups = length(memoffset);
  n_spatial_dims = length(sz_spatial);
  n_nbrs = cellfun(@length,memoffset);
  cum_n_nbrs = cumsum(n_nbrs);
  % Prepare for edge treatment
  n_edge = zeros(n_spatial_dims,n_groups);
  for groupIndex = 1:n_groups
    n_edge(:,groupIndex) = max(coordoffset{groupIndex},[],2);
  end
  n_l = 1+n_edge';
  n_r = repmat(sz_spatial,n_groups,1) - n_edge';
  
  %% Apply the bn condition
  sz_out = sz_spatial;
  if (length(sz_out) == 1)
    sz_out(2) = 1;
  end
  n = zeros(sz_out);
  imfilteredOut = zeros(sz);
  value_skip = (0:n_value_dims-1)';
  coordc = cell(1,n_spatial_dims);
  for pixelIndex = 1:n_pix
    % Do any coordinate preparation
    [coordc{:}] = ind2sub(sz_out,pixelIndex);
    coord = cat(2,coordc{:});
    pixelLookup = (pixelIndex-1)*n_value_dims+1;
    % Look up the value of the current coordinate
    I0 = imfilteredIn(pixelLookup + value_skip);
    % Employ the balanced neighborhood criterion to find the neighbors
    dIsum = zeros(size(I0));
    dIsumsq = zeros(size(I0));
    groupIndex = 1;
    isdone = false;
    while ~isdone
      % extract values of this block of neighbors
      this_memoffset = repmat(memoffset{groupIndex},n_value_dims,1) + ...
        repmat(value_skip,1,n_nbrs(groupIndex));
      I = imdata(pixelLookup+this_memoffset);
      dI = I - repmat(I0,1,size(I,2));
      % evaluate the bn criterion
      dIsum = dIsum + sum(dI,2);
      dIsumsq = dIsumsq + sum(dI.^2,2);
      groupIndex = groupIndex+1;
      isdone = (sum(dIsum.^2) > z2*sum(dIsumsq) || groupIndex > n_groups);
      isdone = isdone || any(coord < n_l(groupIndex,:)) || any(coord > n_r(groupIndex,:));
    end
    % backtrack to find the correct neighbor group
    if (groupIndex <= n_groups)
      n_target = (cum_n_nbrs(groupIndex-1)-z2)/2;
    else
      n_target = cum_n_nbrs(end);
      %n_target = cum_n_nbrs(end)/2;   % we never satisfied the criterion, so don't subtract z2
    end
    [~,groupIndex] = min(abs(n_target - cum_n_nbrs));
    n(pixelIndex) = cum_n_nbrs(groupIndex);
    % evaluate the filtered value
    this_memoffset = repmat(cat(2,memoffset{1:groupIndex}),n_value_dims,1) + ...
      repmat(value_skip,1,n(pixelIndex));
    I = imdata(pixelLookup+this_memoffset);
    imfilteredOut(pixelLookup+value_skip) = mean(I,2);
  end