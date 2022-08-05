function coef = qinterp_grid_inverse(A,mode)
% qinterp_grid_inverse: prepare interpolation coefficients that reconstruct grid values
%
% When quadratic interpolation is performed using the grid values as the
% interpolation coefficients, the interpolated values are not equal to the
% original grid values. This function computes a set of coefficients that
% reconstructs the original grid values. This general approach is described
% in
%   P. Thevenaz, T. Blu, M. Unser. "Interpolation revisted." IEEE
%   Transactions on Medical Imaging, 19: 739-758 (2000).
%
% Syntax:
%   coef = qinterp_grid_inverse(A,mode)
% where
%   A is the array to be interpolated (may be of any dimension)
%   mode is a string which affects how values over the edge of the array
%     are handled. The two possible settings are 'nan' and 'reflect'. See
%     qinterp_grid.
%
% See also: qinterp_grid, imqinterp.

% Copyright 2011 by Timothy E. Holy

  n_dims = ndims(A);
  p = [2:n_dims 1];
  coef = A;
  for dimIndex = 1:n_dims
    n = size(coef,1);
    d = repmat(cast(3/4,class(coef)),1,n);
    dm1 = repmat(cast(1/8,class(coef)),1,n-1);
    dp1 = dm1;
    switch lower(mode)
      case 'nan'
        % do nothing
      case 'reflect'
        d([1 end]) = 7/8;
      otherwise
        error('Mode %s not recognized',mode);
    end
    coef = tridiag_inv(dm1,d,dp1,coef);
    coef = permute(coef,p);
  end
  