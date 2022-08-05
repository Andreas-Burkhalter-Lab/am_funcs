function M = raytrace_matrix
% raytrace_matrix: a collection of utilities for matrix optics
%
% Paraxial rays can be traced using a matrix formulation, where the
% 2-vector
%      v = [n*u; y]
% encodes the ray height (y), refractive index (n), and ray angle (u).
%
% The effect of a surface can be described by a 2-by-2 matrix
%   Msurface = = eye(2,2) + C*N, where N = [0 -(np-n)C; 0 0] 
% where n is the refractive index on the input side, np is the refractive
% input on the output side, and C=1/R where R is the radius of curvature
% (sign convention as in Smith, i.e., C > 0 if the center of curvature is
% to the right).
%
% The effect of propagation a distance L in a medium of refractive index n
% can be described by a 2-by-2 matrix
%   Mpropagation = eye(2,2) + L*P, where P = [0 0; 1/n 0]
%
% You obtain access to these matrices in the following way:
%   M = raytrace_matrix;
% after which
%   M.surface(n,np,C) generates the matrix Msurface
%   M.N(n,np) generates N
%   M.propagation(n,L) generates Mpropagation
%   M.P(n) generates P
%   M.thinlens(f) generates a "thin lens" matrix [1 -1/f; 0 1]
%
% You also have access to the following utilities:
%   M.plot(v,n,x1,x2) draws a line corresponding to the ray v starting at
%     position x1 along the optic axis and propagating to position x2 along
%     the optic axis.

M.surface = @(n,np,C) [1 -(np-n)*C; 0 1];
M.N = @(n,np) [0 -(np-n); 0 0];
M.propagation = @(n,L) [1 0; L/n 1];
M.P = @(n) [0 0; 1/n 0];
M.thinlens = @(f) [1 -1/f; 0 1];
M.plot = @(v,n,x1,x2) line([x1 x2],[0 (x2-x1)/n*v(1)]+v(2));
