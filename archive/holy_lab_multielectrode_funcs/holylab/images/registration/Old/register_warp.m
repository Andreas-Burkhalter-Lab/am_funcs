function [psig,w,sqrtdetJ] = register_warp(psi,g,options)
% REGISTER_WARP: create a warped image
% Syntax:
%   img = register_warp(im,g)
%   [img,w,sqrtdetJ] = register_warp(im,g,options)
% where
%   im is either an image, or the square root of an image (square roots
%     are used during optimization)
%   g is a cell array with the length given by the number of dimensions in
%     im, where g{i} is an array over the grid points of space giving the
%     deformation of the ith coordinate (see REGISTER_G0 for a detailed
%     explanation).
%   options is a structure which may have the following fields:
%    sqrt (default false): if true, the output will be correct assuming im
%      is the square root of an image.
%    covariant (default true): if true, the output intensities will be
%      corrected by the volume factor |detJ|.
%    sqrtdetJ0: if present, this supplies information about a previous
%      deformation with which this one is to be composed.  The
%      returned sqrtdetJ is multiplied by sqrtdetJ0(g(x)).  (Note that,
%      for covariant warping, this factor does not play a role in img,
%      as it is assumed that this previous deformation is incorporated
%      into im.)
% and
%   img is the output warped im;
%   w is a weighting array, of the same size as img, that is 1 for any
%     internal voxel but ranges between 0 and 1 for pixels on the
%     boundary (see IMINTERP);
%   sqrtdetJ is the square root of the determinant of the Jacobian
%     (perhaps multiplied by the interpolated sqrtdetJ0).
%
% See also: IMINTERP,REGISTER_G0.

% Copyright 2006-2007 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'sqrt')
    options.sqrt = false;
  end
  if ~isfield(options,'covariant')
    options.covariant = true;
  end
  if ~isfield(options,'sqrtdetJ0')
    options.sqrtdetJ0 = [];
  end

  if isempty(g{1})
    % Entry with default (identity) deformation. This is a quick exit!
    psig = psi;
    if (nargout > 1)
      w = ones(size(psi),'single')/numel(psi);
    end
    if (nargout > 2)
      if isempty(options.sqrtdetJ0)
        sqrtdetJ = ones(size(psi),'single');
% $$$       % Put NaNs on edges
% $$$       sz = size(sqrtdetJ);
% $$$       n_dims = length(sz);
% $$$       colons = repmat({':'},1,n_dims);
% $$$       for dimIndex = 1:n_dims
% $$$         if (sz(dimIndex) > 1)
% $$$           colons_tmp = colons;
% $$$           colons_tmp{dimIndex} = 1;
% $$$           sqrtdetJ(colons_tmp{:}) = nan;
% $$$           colons_tmp{dimIndex} = sz(dimIndex);
% $$$           sqrtdetJ(colons_tmp{:}) = nan;
% $$$         end
% $$$       end
% $$$       w(isnan(sqrtdetJ)) = 0;
      else
        sqrtdetJ = options.sqrtdetJ0;
      end
    end
    return
  end

  if (nargout > 1)
    [psig,w] = iminterp(psi,g{:});
    w(isnan(psig)) = 0;
    sum_w = sum(w(:));
    if (sum_w == 0)
      error('register:warp','No points to return');
    else
      w = w / sum_w;  % normalize, so that the # of voxels plays no role
    end
  else
    psig = iminterp(psi,g{:});
  end
  
  if (options.covariant || nargout > 2)
    detJ = register_g2detJ(g{:});
    if any(detJ(:) < 0)
      error('register:warp','detJ has negative values');
    end
  end
  
  if (options.sqrt || nargout > 2)
    sqrtdetJ = sqrt(detJ);
  end
  
  if options.covariant
    if options.sqrt
      % We're working with the square-root of an image
      psig = psig .* sqrtdetJ;
    else
      % We're working with a real image
      psig = psig .* detJ;
    end
  end
  
  % "Compose" sqrtdetJ with that from prior deformations
  if (nargout > 2 && ~isempty(options.sqrtdetJ0))
    sqrtdetJ0_g = iminterp(options.sqrtdetJ0,g{:},'extrap');
    sqrtdetJ = sqrtdetJ .* sqrtdetJ0_g;
  end
