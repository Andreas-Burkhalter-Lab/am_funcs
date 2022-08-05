function [snips,snipcoords] = imsnip(xo,im,params)
% IMSNIP: snip out box-shaped volumes in image
% Syntax:
%   snips = imsnip(x,im,params)
%   [snips,snipcoords] = imsnip(x,im,params)
% where
%   x is d-by-npts, specifying the center positions for each snip
%   im is a d-dimensional image
%   params has the following fields:
%     size: a 1-by-d vector giving the range in pixels to use along each
%       axis (i.e., [3 3 5] would use [-3 -2 -1 0 1 2 3] relative to the
%       first & second coords, and -5:5 for the 3rd)
% and
%   snips is a d+1 dimensional matrix, snips(...,i) is the image in the
%     regions of xo(:,i)
%   snipcoords is a cell array, snipcoords{i} is a cell array containing
%     the image coordinates from which the snippet was taken.
%     im(snipcoords{i}{:}) would return the section of the image
%     equal to the ith snippet.

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
  colons = repmat({':'},1,d);
  
  szsnips = 2*params.size+1;
  snips = nan([szsnips npts]);
  snipcoords = cell(1,npts);
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
    % Construct the snipranges
    xrnd = round(x);
    xsnip = cell(1,d);
    for dIndex = 1:d
      xsnip{dIndex} = xrnd(dIndex) + wcoords{dIndex};
      xsnip{dIndex}(xsnip{dIndex} < 1) = [];
      xsnip{dIndex}(xsnip{dIndex} > sz(dIndex)) = [];
    end
    snipcoords{ptIndex} = xsnip;
    % Snip out image
    if isvalid
      snips(colons{:},ptIndex) = im(xsnip{:});
    end      
  end
  