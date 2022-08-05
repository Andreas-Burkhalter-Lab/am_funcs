function I = make_gaussian_image(sz,a,mu,sigma)
% make_gaussian_image: create an image as a sum of Gaussians
%
% This function can be useful in testing certain image processing routines.
% It allows you to place Gaussians at varying points in the image, and even
% permits sub-pixel placement of the individual Gaussians.
%
% Syntax:
%   I = make_gaussian_image(sz,a,mu,sigma)
% where
%   sz is a vector giving the size of the image
%   a is a vector, the ith element giving the amplitude of the ith Gaussian
% and the remaining arguments are as described in make_gaussian_clusters.
%
% I is the output image.

% Copyright 2011 by Timothy E. Holy

  n_dims = length(sz);
  I = zeros(sz);
  I = I(:)';
  x = cell(1,n_dims);
  for i = 1:n_dims
    x{i} = 1:sz(i);
  end
  X = ndgrid_as_matrix(x{:})';
  n_gaussians = length(a);
  % We use a quadratic model to integrate the total intensity over each
  % pixel. The answer turns out to be
  %     N(x) * [1 - Tr(inv(C))/24 + (C\(x-mu))^2/24)]
  % where N(x) is the underlying Gaussian evaluated at the pixel center and
  % C is the covariance matrix.
  for i = 1:n_gaussians
    dX = bsxfun(@minus,X,mu(:,i));
    if isequal(size(sigma),[n_dims n_dims n_gaussians])
      C = sigma(:,:,i);
      C = C*C';
      CidX = C\dX;
      c = a(i)*exp(-sum(dX.*CidX,1)/2)/(2*pi)^(n_dims/2)/sqrt(det(C));
      I = I+c.*(1 - trace(inv(C))/24 + sum(CidX.^2,1)/24);
    else
      if isscalar(sigma)
        s = sigma^2;
      elseif isvector(sigma) && length(sigma) == n_gaussians
        s = sigma(i)^2;
      else
        error('Shape of sigma is not recognized');
      end
      dX2sum = sum(dX.^2,1);
      c = a(i)*exp(-dX2sum/2/s)/(2*pi)^(n_dims/2)/sqrt(s^n_dims);
      I = I + c.*(1 - n_dims/s/24 + dX2sum/s^2/24);
    end
  end
  I = reshape(I,sz);
  
  