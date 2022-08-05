% REGISTER_G2DETJ: calculate the determinant of the Jacobian without storing J
%
% The Jacobian takes a lot of memory: d^2 * the number of pixels in the
% stack.  However, for computing the warped image we only need the
% determinant of the Jacobian.  This function calculates the determinant
% without bothering to store the Jacobian itself, thus saving a lot of
% memory.
%
% Syntax:
%   detJ = register_g2detJ(g1,g2,...)
%   detJ = register_g2detJ(g1,g2,...,handle_edges)
% where
%   g1, g2 are the deformation arrays, for example from g{:} (see
%     REGISTER_G0 for a detailed explanation).
%   handle_edges is a boolean (default true) that, if true, results
%     in finite values being produced on edges (if false, NaNs);
% and
%   detJ is the determinant of the Jacobian at each spatial location.
%
% See also: REGISTER_G0, REGISTER_JACOBIAN, REGISTER_DETJ.

% Implemented as a MEX file.
% Copyright 2006 by Timothy E. Holy
