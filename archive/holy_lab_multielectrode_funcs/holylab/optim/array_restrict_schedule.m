function pyramid = array_restrict_schedule(sz,options)
% ARRAY_RESTRICT_SCHEDULE: sizes of an array under successive restriction
%
% This function is used to predict the size of an array as you successively
% call array_restrict.  This creates an "array pyramid" of coarser and
% coarser versions of the original array. You can control the sequence of
% restrictions along different array dimensions in a variety of ways.
%
% Syntax:
%   pyramid = array_restrict_schedule(sz)
%   pyramid = array_restrict_schedule(sz,options)
% where
%   sz is the size of the array you want to restrict
%   options may have the following fields:
%     min_pixels (default 2): the minimum number of "pixels" along any
%       dimension.
%     truncate (default false): once the minimum number of pixels is
%       reached along any given dimension, if truncate=true we quit making
%       smaller elements of the pyramid. If false, we stop restricting
%       along this dimension, but we keep restricting along the other
%       dimensions until none of the dimensions are large enough for
%       further restriction.
%     restrict_schedule: a logical array of size n_levels-by-n_dims,
%       specifying the precise series of restrictions. For example,
%       if restrict_schedule(i,:) = [1 1 0], then the i+1st grid will be
%       generated from the ith one by restricting along the first two
%       dimensions, but not the 3rd.
%       If supplied, this overrides pixel_spacing in terms of controlling
%       the restriction schedule. The process of restriction stops once
%       this array is exhausted.
%     pixel_spacing: a vector giving the spacing between adjacent "pixels"
%       along each coordinate.  For example, if your pixels are spaced 0.25
%       microns in x and y, and 5 microns in z, this would be [0.25 0.25
%       5]. Default: all 1s.  Include only dimensions of size bigger than
%       1 in this vector.
%       The process of making coarser-resolution images will selectively
%       restrict along dimensions that are of substantially-higher
%       (>sqrt(2)-fold) resolution, so that as you progress to coarser
%       grids all of the dimensions ultimately have similar pixel spacing.
%       Of course, if you want to bias this process (for example, if your
%       volume is much smaller, in physical space, along one axis than any
%       of the others), you can safely "lie" about this so that your
%       coarser-resolution arrays do not "collapse" along the short
%       dimension.
%       This field is used to make restriction decisions only if
%       restrict_schedule is not specified.
% and
%   pyramid is a structure array, with one entry per level of
%     restriction. The first entry corresponds to the top level (finest
%     resolution).  Each entry has a number of fields:
%       sz: describes the size of the array at the given level
%       restrict: a logical vector that specifies how this level was
%         derived from the next finer level.  For example, [1 1 0] would
%         indicate that the first two dimensions were restricted, and the
%         3rd was not.
%       pixel_spacing: the spacing between adjacent pixels at the given
%         level. 
%
% A convenient way to get all the sizes as an array is
%   sz = cat(1,pyramid.sz);
%
% See also: ARRAY_RESTRICT.
  
% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  dimFlag = sz > 1;
  n_dims = sum(dimFlag);
  options = default(options,'min_pixels',2,...
    'truncate',false,...
    'pixel_spacing',ones(1,n_dims));
  
  pixel_spacing = options.pixel_spacing;
  if length(pixel_spacing) ~= n_dims
    error('Include only dimensions of size > 1 in pixel_spacing');
  end
  pixel_spacing_full = zeros(1,length(sz));
  pixel_spacing_full(dimFlag) = pixel_spacing;
  restrict_flag_full = false(1,length(sz));
  
  pyramid = [];
  while true
    %% Create the pyramid entry
    tmppyramid = struct('sz',sz,...
      'restrict',restrict_flag_full,...
      'pixel_spacing',pixel_spacing_full);
    if isempty(pyramid)
      pyramid = tmppyramid;
    else
      pyramid(end+1) = tmppyramid; %#ok<AGROW>
    end
    
    %% Determine which dimensions should be restricted
    if isfield(options,'restrict_schedule')
      % Use the supplied schedule
      if length(pyramid) > size(options.restrict_schedule,1)
        break
      end
      restrict_flag = options.restrict_schedule(length(pyramid),:);
    else
      % Use the pixel spacing to determine the next restriction pattern
      min_spacing = min(pixel_spacing);
      restrict_flag = pixel_spacing < sqrt(2)*min_spacing;
      if ~options.truncate
        % Avoid restricting dimensions that have reached their smallest size
        sz_new = ceil((sz(dimFlag)+1)/2);  % works for sizes > 2
        restrict_flag = restrict_flag & sz_new >= options.min_pixels & sz(dimFlag) > options.min_pixels;
      end
      if ~any(restrict_flag)
        break  % we're all done, quit
      end
    end
    
    %% Calculate the sizes and pixel spacing of the next grid
    restrict_flag_full(dimFlag) = restrict_flag;
    sz_old = sz(restrict_flag_full);
    sz_new = ceil((sz_old+1)/2);
    sz_new(sz_old == 2) = 1;  % fix the problem with sz = 2
    sz(restrict_flag_full) = sz_new;
    pixel_spacing(restrict_flag) = pixel_spacing(restrict_flag)*2;  % I think 2 is more appropriate than sz_old./sz_new, but this is debatable
    pixel_spacing_full(dimFlag) = pixel_spacing;
    if (options.truncate && any(sz(dimFlag) < options.min_pixels))
      break
    end
    if any(sz(dimFlag) < options.min_pixels)
      error('The next element is too small along at least one dimension. Did you supply an incorrect restriction schedule?');
    end
  end
