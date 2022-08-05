function [u,rmg_params] = register_change_level_gap(u,rmg_params,direction)
% REGISTER_CHANGE_LEVEL_GAP: change scale of deformation grid in registration
% 
% If you want to decrease/increase g_level_gap (see
% REGISTER_MULTIGRID_OPTIONS) for a registration that you've already made
% progress on, use this function.
%
% Syntax:
%   [u,rmg_params] = register_refine_u(u,rmg_params,direction)
% where
%   u is your existing estimate of the deformation
%   rmg_params contains your current settings (see
%     REGISTER_MULTIGRID_OPTIONS), 
%   direction is a character, 'f' to make u finer and 'c' to make u
%     coarser.
%
% See also: REGISTER_MULTIGRID_OPTIONS.

% Copyright 2010 by Timothy E. Holy

gap = rmg_params.gap_data;
sz = size(u);
n_dims = sz(end);
colons = repmat({':'},1,n_dims);
switch(lower(direction))
  case 'f'
    u = register_u_prolong(u,rmg_params,gap);
    rmg_params.gap_data = gap - 1;
  case 'c'
    restrictFlag = rmg_params.image_grid(gap+1).restrict;
    for dimIndex = 1:n_dims
      u(colons{:},dimIndex) = (1.0/(1+restrictFlag(dimIndex)))*array_restrict(u(colons{:},dimIndex),restrictFlag);
    end
    rmg_params.gap_data = gap + 1;
  otherwise
    error(['Direction ' direction ' not recognized'])
end
