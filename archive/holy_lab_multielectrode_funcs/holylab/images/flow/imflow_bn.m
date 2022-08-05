function [dx,map] = imflow_bn(imdata,n,memoffset,coordoffset)
% [dx,map] = imflow_bn(imdata,n,memoffset,coordoffset)
% where
%   imdata is the underlying image;
%   n is the first output from bn_image;
%   memoffset and coordoffset are the outputs from pixel_neighbors.
%
% See also: bn_image_filter, pixel_neighbors.

% Copyright 2010 by Timothy E. Holy

  %% Parse the inputs
  if iscell(memoffset)
    memoffset = cat(2,memoffset{:});
  end
  if iscell(coordoffset)
    coordoffset = cat(2,coordoffset{:});
  end
  n_spatial_dims = size(coordoffset,1);
  sz = size(imdata);
  n_dims = length(sz);
  n_value_dims = n_dims-n_spatial_dims;
  if n_value_dims > 0
    n_values = prod(sz(1:n_value_dims));
  else
    n_values = 1;
  end
  sz_spatial = sz(end-n_spatial_dims+1:end);
  n_pix = prod(sz_spatial);
  
  %% Calculate the flow
  sz_out = sz_spatial;
  if (length(sz_out) == 1)
    sz_out(2) = 1;
  end
  value_skip = (0:n_values-1)';
  displacement_skip = (0:n_spatial_dims-1)*n_pix;
  memoffset_all = repmat(memoffset,n_values,1) + ...
    repmat(value_skip,1,size(memoffset,2));
  dx = nan([sz_out,n_spatial_dims]);
  for pixelIndex = 1:n_pix
    if (n(pixelIndex) == 0)
      continue
    end
    pixelLookup = (pixelIndex-1)*n_values+1;
    I0 = imdata(pixelLookup+value_skip);
    thisn = n(pixelIndex);
    this_memoffset = memoffset_all(:,1:thisn);
    I = imdata(pixelLookup+this_memoffset);
    %dI = I - repmat(I0,1,size(I,2));
    dI = bsxfun(@minus,I,I0);
    dImean = mean(dI,2);
    % compute the weights as a dot product of each displacement with the
    % mean displacement
    %w = sum(dI .* repmat(dImean,1,size(dI,2)),1);
    w = sum(bsxfun(@times,dI,dImean),1);
    % shift so smallest weight is zero (avoids negative weights)
%     w = w - min(w);
%     % normalize by their sum
%     if (sum(w) > 0)
%       w = w/sum(w);
%     end
    wsum = sum(abs(w));
    if (wsum > 0)
      w = w / wsum;
    end
    % use these weights to calculate the displacement
    %thisdx = sum(repmat(w,n_spatial_dims,1) .* coordoffset(:,1:thisn),2);
    thisdx = sum(bsxfun(@times,w,coordoffset(:,1:thisn)),2);
    % store the result
    dx(pixelIndex+displacement_skip) = thisdx;
  end
  if (nargout > 1)
    map = reshape(1:n_pix,sz_out);
    %maxdx = max(abs(dx),[],ndims(dx));
    maxdx = ones(size(map));
    isnz = maxdx > 0 & n > 0;
    pixel_skip = [1 cumprod(sz_out)];
    colons = repmat({':'},1,n_spatial_dims);
    for dimIndex = 1:n_spatial_dims
      thisdx = dx(colons{:},dimIndex);
      map(isnz) = map(isnz) - round(thisdx(isnz)./maxdx(isnz))*pixel_skip(dimIndex);
    end
  end
  