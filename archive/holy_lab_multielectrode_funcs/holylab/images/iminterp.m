% iminterp: interpolate an image quickly
% Syntax:
%   imout = iminterp(im,X1,X2,...)
%   [imout,w] = iminterp(im,X1,X2,...)
%   imout = iminterp(im,X1,X2,...,'extrap')
% where
%   im is the input image (may be multidimensional)
%   X1, X2, ... are arrays of the same number of dimensions as the image
%     (but not necessarily the same size), specifying a set of evaluation
%     locations.  These arrays specify a deformation X1 = g_1(x) for x over
%     the points on a grid.
% and
%   imout is the estimated image value at the grid points specified by
%     g.  Points lying outside the image region are returned as NaN, with
%     a limited exception occuring when you request w;
%   w is a weight array with the same size as the image.  w is 1 in
%     the interior but between 0 and 1 on the boundary.  For example,
%     in d=1, a pixel evaluated at x = 0.8 is just beyond the left edge
%     of the image.  Without the w output, the interpolated value would
%     be NaN.  With the w output, the pixel would have the value of the
%     pixel at x=1 but would have a weight of 0.8.  This is useful in
%     ensuring continuity of functions of interpolated images (e.g., in
%     image registration).
%
% If you want to have the image value linearly extrapolated beyond
% the edge of the image, supply the extra string argument
% 'extrap'. This mode is incompatible with the "w" weight array
% output.
%
% Note: both im and the Xi must be of the same data type (e.g., single or
% double).
%
% See also: INTERP1, INTERPN.

% Copyright 2006-2007 by Timothy E. Holy

% Implemented as a MEX function
