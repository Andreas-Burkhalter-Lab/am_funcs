function [dx,val] = array_findpeak(A,options)
% array_findpeak: location of the peak of an array
%
% This function finds the location of the peak of a function defined on a
% grid. Optionally, one can:
%   1. Smooth the array before finding the peak
%   2. Do subpixel interpolation
%   3. Limit the search to just the center region of the array.
%
% Syntax:
%   [dx,val] = array_findpeak(A)
%   [dx,val] = array_findpeak(A,options)
% where
%   A is an array (of arbitrary dimensionality)
%   options is a structure which may have the following fields:
%     smooth: an optional function used to smooth the array. This has the syntax
%         As = ops.smooth(A)
%       In many applications, a good choice is
%         h = fspecial('gaussian',5,1);
%         ops.smooth = @(im) imfilter(im,h);
%     subpixel_mode (default 'interp'): controls whether the peak is
%       defined off-grid or not. Choices are '' (peak is constrained to be
%       at a grid point), 'interp' (use quadratic interpolation to find the
%       maximum), 'com') (compute the center of mass of a block 1 "pixel"
%       on either side of the on-grid peak).
%     dx_max: if supplied, only searches for maximum within a region of the
%       center of the array.  For example, dx_max = [3 5 2] would search in
%       a rectangular region that is 3 "pixels" from the center along the
%       first coordinate, 5 along the second, and 2 along the third (so the
%       size of the rectangular region would be 7-by-11-by-5).
% and
%   dx is the estimated location of the peak. This is reported relative to
%     the array center, ceil((size(A)+1)/2).
%   val is the value of the array at the peak. When options.subpixel_mode =
%     'interp', this includes the effects of quadratic interpolation.
%     
% A primary application of this function is for registration by phase
% correlation.
%
% See also: register_translate.

% Copyright 2011 by Timothy E. Holy

  sz = size(A);
  n_dims = length(sz);
  szFlag = sz > 1;
  n_dims1 = sum(szFlag);
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'subpixel_mode','interp');
  
  center = ceil((sz+1)/2);
  if isfield(options,'smooth') && ~isempty(options.smooth)
    A = options.smooth(A);
  end
  if isfield(options,'dx_max')
    % Snip out a region around the center
    dx_max = min((size(A)-1)/2,abs(options.dx_max));
    dx_max = ceil(dx_max);
    xsnip = cell(1,n_dims);
    for i = 1:n_dims
      xsnip{i} = center(i) + (-dx_max(i):dx_max(i));
    end
    A = A(xsnip{:});
    sz = size(A);
    center = ceil((sz+1)/2);
  end
  % We have to start from a point that has no undefined neighbors
  Avalid = isfinite(A);
  nhood = true(repmat(3,1,ndims(A)));
  Avalid = imerode(Avalid,nhood);
  Avalid = killedges(Avalid,false);
  Atmp = A;
  Atmp(~Avalid) = -inf;
  [val,indx] = max(Atmp(:));
  coords = cell(1,n_dims1);
  [coords{:}] = ind2sub(sz,indx);
  dx = zeros(1,length(sz));
  if ~isempty(options.subpixel_mode) && all(sz(szFlag) > 2)
    % We're looking for the peak on a resolution finer than the grid spacing.
    % Snip out the region one pixel on either side of the peak.
    rng = cell(1,n_dims1);
    for dimIndex = 1:n_dims1
      if szFlag(dimIndex)
        rng{dimIndex} = coords{dimIndex} + [-1 0 1];
      else
        rng{dimIndex} = 1;
      end
    end
    Asnip = A(rng{:});
    switch options.subpixel_mode
      case 'interp'
        % Find the peak of the quadratic interpolation
        minops = optimset('GradObj','on','Display','off');
        func = @(g) imqinterp(g,-double(Asnip));
        one = ones(1,n_dims1);
        [dx,val] = fmincon(func,one+1,[],[],[],[],1.501*one,2.499*one,[],minops); % fractional shift
        dx = dx - 2;
        val = -val;
      case 'com'
        % Find the center of mass
        dx = array_com(Asnip);
      otherwise
        error('subpixel_mode not recognized');
    end
  end
  dx = cat(2,coords{:}) - center(szFlag) + dx;
