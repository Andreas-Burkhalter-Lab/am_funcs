function imMg = register_multigrid_warp(imM,u,rmg_params)
% REGISTER_MULTIGRID_WARP: convenience function for nonrigid image deformations
% Syntax:
%   imMg = register_multigrid_warp(imM,u,rmg_params)
% where
%   imM is the moving image;
%   u is the deformation, defined on a smaller grid (see below);
%   rmg_params defines the grids;
% and
%   imMg is the warped "image".
%
% u is defined on a grid that is rmg_params.gap_data levels below the
% top level (at which imM is defined).
%
% See also: REGISTER_MISMATCH_NONCOV.
  
% Copyright 2010 by Timothy E. Holy
  
  g0 = [];
  if (nargin > 2)
    if ~isempty(u)
      % Prolong u up to full-size
      u = register_u_prolong(u,rmg_params);
    end
    if isfield(rmg_params.image_grid,'g0')
      g0 = rmg_params.image_grid(1).g0;
    end
  end
  
  if isempty(g0)
    if (~isempty(u) || isfield(rmg_params,'shift'))
      g0c = register_g0(size(imM),class(imM));
      g0 = cat(ndims(imM)+1,g0c{:});
    end
  end
  
  if (nargin > 2)
    if isfield(rmg_params,'shift')
      colons = repmat({':'},1,rmg_params.n_dims);
      for dimIndex = 1:rmg_params.n_dims
        g0(colons{:},dimIndex) = g0(colons{:},dimIndex)+rmg_params.shift(dimIndex);
      end
    end
  end
    
  if ~isempty(u)
    g = g0+u;
    imMg = imqinterp(g,imM);
  else
    imMg = imM;
  end
