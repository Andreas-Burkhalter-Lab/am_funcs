function uH = register_rigid2nonrigid(dx,rmg_params)
% REGISTER_RIGID2NONRIGID: convert translation into a multigrid-able deformation
% Syntax:
%   uH = register_rigid2nonrigid(dx,rmg_params)
% where
%   dx is the shift in pixels, e.g., the output of register_rigid;
%   rmg_params is the registration structure computed by
%     register_multigrid_options
% and
%   uH is an appropriate initial guess for the deformation in routines like
%     register_multigrid_vcycle.
%
% See also: REGISTER_RIGID, REGISTER_MULTIGRID_VCYCLE.

% Copyright 2010 by Timothy E. Holy

  n_dims = rmg_params.n_dims;
  n_levels = length(rmg_params.image_grid);
  restrictFlag = cat(1,rmg_params.image_grid.restrict);
  cdF = cumsum(restrictFlag,1);
  level = rmg_params.gap_data+1;
  scale = 2.^cdF(level,:);
  sz = rmg_params.image_grid(level).sz;
  uH = zeros([sz n_dims],class(rmg_params.image_grid(1).imFixed));
  colons = repmat({':'},1,n_dims);
  for i = 1:n_dims
    uH(colons{:},i) = dx(i)/scale(i);
  end
end