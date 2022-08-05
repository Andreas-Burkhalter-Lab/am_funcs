function [v,g] = qinterp_grid(x,vgrid,mode)
% qinterp_grid: quadratic interpolation of a regular grid, with optional extrapolation
% This is essentially the equivalent of the MEX file imqinterp, but it
% allows extrapolation.
%
%   v = qinterp_grid(x,vgrid,mode)
% performs quadratic interpolation of the values vgrid, defined on a
% regular grid, to calculate the values at positions specified by x. x is a
% n_pts-by-n_dims matrix, vgrid is an array with at least n_dims dimensions
% (if there are further dimensions, these are treated as additional
% independent grids).
% mode is a string which may assume the following values:
%   'nan': beyond-the-edge grid points are given a value of NaN
%   'reflect': mirror-symmetric boundary conditions on vgrid are used
%
%  [v,g] = qinterp_grid(x,vgrid,mode)
% also evaluates the gradient of the array at each point x. g is a cell
% array, where each element contains the gradient with respect to the
% corresponding spatial coordinate.
%
% Note: if the output points are also on a grid, you may be in a situation
% where you can instead use image_snip_qinterp.
%
% Also note: quadratic interpolation does not preserve on-grid values, and
% (like linear interpolation) tends to smooth the input array. You may want
% to first call qinterp_grid_inverse on the array of coefficients to avoid
% these two "problems."
%
% See also: imqinterp, image_snip_qinterp, qinterp_grid_inverse, qinterp_coef.

% Copyright 2011 by Timothy E. Holy

  %% Input parsing
  calc_grad = nargout > 1;
  if (ndims(x) > 2)
    error('The positions must be supplied as a n_pts-by-n_dims matrix');
  end
  [n_pts,n_dims] = size(x);
  sz = size(vgrid);
  sz(end+1) = 1;
  sz_spatial = sz(1:n_dims);
  checknan = false;
  switch lower(mode)
    case 'nan'
      edgefunc = @(x) edgenan(sz_spatial,x);
      checknan = true;
    case 'reflect'
      edgefunc = @(x) edgereflect(sz_spatial,x);
    otherwise
      error('Mode %s not recognized',mode);
  end
  n_grid_pts = prod(sz_spatial);
  n_values = prod([1 sz(n_dims+1:end)]);
  vgrid = reshape(vgrid,[n_grid_pts n_values]);
  
  %% Allocate storage
  v = zeros(n_pts,n_values,class(vgrid));
  coef = zeros(n_pts,n_values,class(vgrid));
  if calc_grad
    g = cell(1,n_dims);
    for dimIndex = 1:n_dims
      g{dimIndex} = zeros([n_pts n_values],class(vgrid));
    end
  end
  z_all = ind2sub_matrix(repmat(3,1,n_dims),1:3^n_dims)-2;
  
  %% Perform interpolation by looping over displacements from center point
  xr = round(x);
  alpha = x - xr;
  for zIndex = 1:size(z_all,1)
    z = z_all(zIndex,:);
    xz = bsxfun(@plus,xr,z);  % integer position
    xz = edgefunc(xz);        % apply boundary conditions
    xI = sub2ind_matrix(sz_spatial,xz,false);
    % Calculate the interpolation coefficients
    for dimIndex = 1:n_dims
      if (z(dimIndex) == 0)
        coef(:,dimIndex) = 3/4-alpha(:,dimIndex).^2;
      else
        coef(:,dimIndex) = (alpha(:,dimIndex)+z(dimIndex)/2).^2/2;
      end
    end   
    if checknan
      nanflag = isnan(xI);
      xI(nanflag) = 1;
      thisv = vgrid(xI,:);
      thisv(nanflag) = nan;
    else
      thisv = vgrid(xI,:);
    end
    % Add the contribution of this displacement to the total result
    v = v + bsxfun(@times,prod(coef,2),thisv);
    % Calculate the gradient terms, if applicable
    if calc_grad
      for dimIndex = 1:n_dims
        tmpcoef = coef;
        if (z(dimIndex) == 0)
          tmpcoef(:,dimIndex) = -2*alpha(:,dimIndex);
        else
          tmpcoef(:,dimIndex) = alpha(:,dimIndex)+z(dimIndex)/2;
        end
        g{dimIndex} = g{dimIndex} + ...
          bsxfun(@times,prod(tmpcoef,2),thisv);
      end
    end
  end
  % Reshape the value coordinates to original format
  v = reshape(v,[n_pts sz(n_dims+1:end)]);
  if calc_grad
    for dimIndex = 1:n_dims
      g{dimIndex} = reshape(g{dimIndex},[n_pts sz(n_dims+1:end)]);
    end
  end
  
function xo = edgenan(sz,x)
  xo = x;
  xo(x < 1) = NaN;
  xo(bsxfun(@gt,x,sz)) = NaN;

  