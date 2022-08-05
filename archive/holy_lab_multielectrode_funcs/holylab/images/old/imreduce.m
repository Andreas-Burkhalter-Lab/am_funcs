function imout = imreduce(im, decfactors, oversmooth)
% IMREDUCE: create smoothed, lower-resolution versions of images
%
% Unlike Matlab's IMRESIZE, this function works to shrink over dimensions
% higher than 2.
%
% Syntax:
%   imout = imreduce(im, decfactors)
%   imout = imreduce(im, decfactors, oversmooth)
% where
%   im is the image (may be multidimensional);
%   decfactors is a vector of integers, yielding the factor by which
%     reduction is to be applied in each coordinate. For example,
%             decfactors = [4 4 1]
%     would result in 4-fold binning along the first and second
%     coordinates, and no binning along the third.
%   oversmooth is a vector that specifies a factor of extra smoothing
%     applied to the images. Supply values > 1 along a particular
%     coordinate for extra smoothing, and values between 0 and 1 for
%     undersmoothing (but you'll get aliasing) (default: [1 1 1]).
% and
%   imout is the smoothed, sub-sampled version of the image.
%
% See also: PYRAMIDS_IMAGINE, IMFILTER_GAUSSIAN.
  
% Copyright 2006 by Timothy E. Holy
  
  if isscalar(decfactors)
    decfactors = repmat(decfactors,[1 ndims(im)]);
  end
  if (ndims(im) > length(decfactors))
    error('Number of dimensions do not agree');
  end
  if (ndims(im) < length(decfactors))
    % The final dimensions are "1"
    decfactors = decfactors(1:ndims(im));
  end
  if ~isequal(decfactors,round(decfactors))
    error('decfactors must be a vector of integers');
  end
  % If user supplied 0 for any coordinate, assume they meant 1
  decfactors(decfactors == 0) = 1;
  if (nargin < 3)
    oversmooth = ones(size(decfactors));
  end
  % Set the smoothing length scale
  sigma = decfactors/2 .* oversmooth;
  % If we're keeping all the points in any dimension, don't
  % smooth in this dimension
  sigma(decfactors == 1) = 0;
  %imout = imfilter_gaussian(im,sigma,'replicate');
  imout = imfilter_gaussian(im,sigma);
  % Now subsample
  coords = cell(1,length(decfactors));
  for i = 1:length(decfactors)
    coords{i} = 1+round((decfactors(i)-1)/2):decfactors(i):size(im,i);
    if isempty(coords{i})
      coords{i} = round(mean([1 size(im,i)]));
    end
  end
  imout = imout(coords{:});
  