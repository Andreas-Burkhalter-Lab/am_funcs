function xnew = immeanshift(xo,im,params)
% IMMEANSHIFT: shift points to the local center of mass of an image
% Syntax:
%   xnew = immeanshift(xo,im,params)
% where
%   xo, xnew are d-by-npts
%   im is a d-dimensional image
%   params has the following fields:
%     size: a 1-by-d vector giving the range in pixels to use along each
%       axis (i.e., [3 3 5] would use [-3 -2 -1 0 1 2 3] relative to the
%       first & second coords, and -5:5 for the 3rd)
%     sigma: a 1-by-d vector giving the width of the gaussian along each
%       axis

% Copyright 2007 by Timothy E. Holy

  [d,npts] = size(xo);
  if (ndims(im) ~= d)
    error('Dimensionality of points and image do not match');
  end
  sz = size(im);
  % Make the coordinate index for the weights, but put them in their
  % appropriate dimension
  wcoords = cell(1,d);
  for dIndex = 1:d
    sz1 = ones(1,d);
    sz1(dIndex) = 2*params.size(dIndex)+1;
    wcoords{dIndex} = reshape(-params.size(dIndex):params.size(dIndex),sz1);
  end
  sigma2 = 2*params.sigma.^2;
  
  xnew = nan(size(xo));
  xsnip = cell(1,d);
  for ptIndex = 1:npts
    x = xo(:,ptIndex);
    % First check to see if the current point will go over the edge; if
    % so, discard it by filling the output with NaNs
    isvalid = true;
    dIndex = 1;
    while (dIndex <= d && isvalid)
      if (x(dIndex) - params.size(dIndex) < 1 || ...
          x(dIndex) + params.size(dIndex) > sz(dIndex))
        isvalid = false;
      end
      dIndex = dIndex+1;
    end
    if ~isvalid
      continue
    end
    % The point is valid. Construct the gaussian weight matrix (with
    % sub-pixel resolution)
    xrnd = round(x);
    xfrac = x-xrnd;
    w = ones(2*params.size + 1);
    for dIndex = 1:d
      sz1 = 2*params.size+1;
      sz1(dIndex) = 1;
      w = w .* repmat(exp(-(wcoords{dIndex}-xfrac(dIndex)).^2/sigma2(dIndex)),sz1);
    end
    % Snip out the image
    for dIndex = 1:d
      xsnip{dIndex} = xrnd(dIndex) + wcoords{dIndex};
    end
    imsnip = single(im(xsnip{:}));
    % Compute the denominator
    I0 = imsnip .* w;
    I0 = sum(I0(:));
    % Compute the terms of the numerator and then the mean
    for dIndex = 1:d
      sz1 = 2*params.size+1;
      sz1(dIndex) = 1;
      tmp = repmat(xsnip{dIndex},sz1) .* w .* imsnip;
      xnew(dIndex,ptIndex) = sum(tmp(:))/I0;
    end
  end
  