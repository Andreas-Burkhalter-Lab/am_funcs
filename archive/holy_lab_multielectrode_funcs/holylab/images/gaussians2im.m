function [im,model,im_coords] = gaussians2im(I,mu,C,sz)
% GAUSSIANS2IM: calculate a model image, given gaussian parameters
% Syntax:
%   [im,model,im_coords] = gaussians2im(I,mu,C,sz)
% where
%   I, mu, and C are the outputs of IM2GAUSSIANS;
%   sz is a 1-by-n_dims size vector of the ouptut image
% and
%   im is the output image
%   model is a 1-by-n_objects cell array, containing pixel intensities for
%     the model gaussian in the region of the object; note the pixels can
%     "spill over" the pixels that minimally define the object;
%   im_coords is a 1-by-n_objects cell array, im_coords{i} is a set of
%     indices into the image that correspond to the given object.
%
%  Note that the outputs "model" and "im_coords" are designed to be fed
%  into IM_GAUSSIANS_ERR.
%
% See also: IM2GAUSSIANS, IM_GAUSSIANS_ERR.

% Copyright 2006 by Timothy E. Holy

  im = zeros(sz,'single');
  n_objects = length(I);
  n_dims = size(mu,1);
  
  x = cell(1,n_dims);
  X = cell(1,n_dims);  
  model = cell(1,n_objects);
  im_coords = cell(1,n_objects);
  
  for objIndex = 1:n_objects
    % Create coordinates around the gaussian
    Ctmp = C(:,:,objIndex);
    distmax = 3*sqrt(max(Ctmp(:)));
    for dimIndex = 1:n_dims
      x{dimIndex} = max(1,floor(mu(dimIndex,objIndex) - distmax)) : ...
        min(sz(dimIndex),ceil(mu(dimIndex,objIndex) + distmax));
    end
    X = cell(1,n_dims);
    [X{:}] = ndgrid(x{:});
    % Index these coordinates to the image coordinates
    indexImage = sub2ind(sz,X{:});
    
    % Calculate the gaussians
    Cinv = inv(C(:,:,objIndex));
    for dimIndex = 1:n_dims
      dX{dimIndex} = X{dimIndex} - mu(dimIndex,objIndex);
    end
    szX = size(X{1});
    d2 = zeros(szX,'single');
    for dimIndex1 = 1:n_dims
      for dimIndex2 = 1:n_dims
        % The covariance-weighted distance
        d2 = d2 + dX{dimIndex1} .* (Cinv(dimIndex1,dimIndex2) * ...
          dX{dimIndex2});
      end
    end
    model{objIndex} = (I(objIndex)*sqrt(det(Cinv)) / (2*pi)^(n_dims/2)) * exp(-d2/2);
    im_coords{objIndex} = indexImage;
    im(indexImage) = im(indexImage) + model{objIndex};
  end
  
    