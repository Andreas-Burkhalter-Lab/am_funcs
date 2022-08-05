function [z,coef,coefg] = qinterp_coef(x)
% qinterp_coef: quadratic interpolation coefficients
%
%   [z,coef] = qinterp_coef(x)
% Returns the displacements (taking values of -1, 0, or 1), in z, and the
% coefficients, in coef, for interpolation at fractional grid points given
% by x. The integer component of x is discarded. x must be a
% n_pts-by-n_dims matrix. z is an n_dims-by-3^n_dims matrix, giving the
% coordinate displacements to all neighbors of the center point. coef is an
% n_pts-by-3^n_dims matrix, giving the coefficient by which to multiply
% the value of each neighbor. One can calculate the interpolated value as
% the sum across neighbors of these products.
%
%   [z,coef,coefg] = qinterp_coef(x)
% also returns the coefficients for calculating the gradient at x. coefg is
% a cell array, one element for each dimension.
%
% Example:
%   n_dims = 3;
%   % Create a random array
%   r = rand((1:n_dims)+2);
%   % Pick a random point in the interior of the array
%   % (must be more than 0.5 away from the edges)
%   x = (size(r)-2).*rand(1,n_dims) + 1.5;
%   xr = round(x);  % the "center" point
%   [z,coef,coefg] = qinterp_coef(x);
%   zr = bsxfun(@plus,z',xr);  % the coordinates of the neighbors
%   zindx = sub2ind_matrix(size(r),zr)';
%   value = sum(coef.*r(zindx))
%   for i = 1:n_dims
%     deriv = sum(coefg{i}.*r(zindx))
%   end
%
% See also: qinterp_grid, imqinterp.

% Copyright 2011 by Timothy E. Holy

  %% Input parsing
  calc_grad = nargout > 2;
  if (ndims(x) > 2)
    error('The positions must be supplied as a n_pts-by-n_dims matrix');
  end
  [n_pts,n_dims] = size(x);
  
  %% Allocate storage
  n_z = 3^n_dims;
  z = ind2sub_matrix(repmat(3,1,n_dims),1:n_z)'-2;
  coef = zeros([n_pts,n_dims,n_z],class(x));
  if calc_grad
    coefgtmp = coef;
    coefg = cell(1,n_dims);
  end
  
  %% Calculate the 1d coefficients
  xr = round(x);
  alpha = x - xr;
  for thisz = -1:1
    mask = (z == thisz);
    [maski,~] = find(z == thisz);  % the coordinate index
    if (thisz == 0)
      tmp = alpha(:,maski);
      coef(:,mask) = 3/4 - tmp.^2;
      if calc_grad
        coefgtmp(:,mask) = -2*tmp;
      end
    else
      tmp = alpha(:,maski)+thisz/2;
      coef(:,mask) = tmp.^2/2;
      if calc_grad
        coefgtmp(:,mask) = tmp;
      end
    end
  end
  
  %% Compute the products of coefficients
  if calc_grad
    for dimIndex = 1:n_dims
      tmp = coef;
      tmp(:,dimIndex,:) = coefgtmp(:,dimIndex,:);
      coefg{dimIndex} = reshape(prod(tmp,2),[n_pts n_z]);
    end
  end
  coef = reshape(prod(coef,2),[n_pts n_z]);

  