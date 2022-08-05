% function [ucompc,detJout,val,grad] = register_logdetpenalty_composition(uoldc,unewc,c)
% 
% register_logdetpenalty_composition: explicit composition with a log(det(J)) penalty
%
% This function composes two deformations and computes a penalty based on
% the local volume change. Optionally, it also computes the derivative of
% the penalty with respect to the "new" deformation.
%
% The full deformation g(x) is defined as
%    g(x) = x + u(x)
% so that u(x) is the displacement from the identity. The composition is
%    gcomp(x) = gold(gnew(x))
% for two deformations, the "old" (original) and "new" (typically, a
% refinement of the original).
%
% The penalty is the mean (log(det(J)/c))^2, where J is the Jacobian of the
% deformation in an individual "cell" (i.e., the set of grid points
% surrounding the half-grid centers).
%
%   [ucompc,detJ,val,grad] = register_logdetpenalty_composition(uoldccoef,unewc)
%   [ucompc,detJ,val,grad] = register_logdetpenalty_composition(uoldccoef,unewc,c)
% This calculates the composed deformation (ucompc), the determinant of J
% in each "cell", the value of the regularization penalty (val), and its
% gradient with respect to unew. uoldc and unewc are cell arrays of length
% n_dims, giving the displacement at each grid point (each element is the
% size of the grid).  Important note: uoldccoef are the coefficients for
% quadratic interpolation of uold, as produced by qinterp_grid_inverse.
%
%   [ucompc,detJ,val,grad] = register_logdetpenalty_composition([],unewc,c)
% uses the identity map as the "old" deformation.
%
% See also: qinterp_grid_inverse, register_block_penalty.

% Copyright 2011 by Timothy E. Holy
  
