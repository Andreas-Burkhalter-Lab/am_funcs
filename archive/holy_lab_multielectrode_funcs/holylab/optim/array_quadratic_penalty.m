% array_quadratic_penalty: compute value and gradient of quadratic function on unit cells of an array
% Penalty functionals such as the gradient penalty,
%      E = \int dx (grad f(x))^2,
% are quadratic on "unit cells" of an array A that is a discrete
% representation of the function f(x). This function makes it easy to
% compute the value and gradient of such penalty functions.
%
% Syntax:
%   [v,g] = array_quadratic_penalty(A,vertex,Q,boundary_conditions)
% where
%   A is a d-dimensional array;
%   vertex is a 2^d-by-d matrix of 0s and 1s, giving the displacement
%     of each "corner" of a unit cell (this matrix defines how Q is
%     interpreted);
%   Q is a 2^d-by-2^d matrix, specifying the penalty on a single unit cell:
%        Pcell = (1/2) * a^T * Q * a,
%     where a would be a vector of length 2^d representing the values
%     of A at the vertices of a single unit cell, in the order
%     specified by the variable "vertex";
%   boundary_conditions (optional) is a string, one of 'none'
%     (default) or 'circular', specifying how boundaries are supposed
%     to be treated.  'none' means that the penalty is applied only
%     to cells that lie fully within the bounds of the array; 'zero'
%     (not yet implemented) assumes that the function values are zero
%     beyond the edge of the array; 'circular' uses periodic boundary
%     conditions.
% On output,
%   v is the value of the penalty, summed over all unit cells (a scalar);
%   g is the gradient of v with respect to the coordinates of A, and
%     therefore is an array of the size of A.
%
% For example, in two dimensions, if we wanted to use a penalty function
%     Pcell = (1/2) (grad A)^2,
% then
%   vertex = [0 0;
%             0 1;
%             1 0;
%             1 1];
% which indicates a unit cell vertex numbering scheme like this:
%
%        3     4
%
%        1     2
%
% (In this coordinate system, we'll call the second coordinate "x",
% which increases to the right; y increases upward.)
% The x-gradient is
%       (a(2)+a(4)-a(1)-a(3))/2
% (the average value on the two right vertices, minus the average
% value on the two left vertices), and the y-gradient is
%       (a(3)+a(4)-a(1)-a(2))/2.
% Let gx = [-1/2, 1/2, -1/2, 1/2]^T, gy = [-1/2,-1/2,1/2,1/2]^T, in
% terms of which the gradients are expressed as a dot product.  We
% have
%       Q = gx gx^T + gy gy^T,
% and hence
%            [ 1/2   0    0  -1/2 ] 
%      Q =   [  0   1/2 -1/2   0  ]
%            [  0  -1/2  1/2   0  ]
%            [-1/2   0    0   1/2 ]
%
% See also: gradient_penalty.

% MEX file. Copyright 2011 by Timothy E. Holy
