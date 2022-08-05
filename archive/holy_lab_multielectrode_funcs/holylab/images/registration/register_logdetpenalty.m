function [detJnew,val,grad,hess] = register_logdetpenalty(u,c,detJprev)
% REGISTER_LOGDETPENALTY: penalize deformations based on log(det(Jacobian))
% Syntax:
%   [detJnew,val,grad] = register_logdetpenalty(u)
%   [detJnew,val,grad] = register_logdetpenalty(u,c)
%   [detJnew,val,grad] = register_logdetpenalty(u,c,detJprev)
%   [detJnew,val,grad,hess] = register_logdetpenalty(...)
% where
%   u is the deformation (an array of size [gridsize n_dims])
%   pixel_spacing is the grid spacing along each axis
%   c (default 1) is the desired value for detJ (1 = no volume change)
%   detJprev, if supplied, allows one to compose deformations together,
%     phitot(x) = phi1(phi2(x)).  Here detJprev is det(J(phi1(x))), and u
%     corresponds to phi2.  Omit this argument if there is no "earlier"
%     deformation.
% and
%   detJnew is the array of jacobian determinants at half-grid points;
%   val = sum((log(det(Jtot))).^2)
%   grad is the gradient of val with respect to u
%   hess is the Hessian with respect to u.
  
% Copyright 2010 by Timothy E. Holy. Implemented as a MEX file.
