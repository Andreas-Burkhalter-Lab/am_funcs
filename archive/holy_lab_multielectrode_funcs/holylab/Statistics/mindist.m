% MINDIST: Find nearest neighbor pairs across sets of points
% Usage:
%  [dist,index] = mindist(x,y)
% where
%   x is a d-by-N matrix of points in d-dimensional space;
%   y is a d-by-q matrix of points in d-dimensional space;
% and
%   dist is a 1-by-N vector giving, for each point in x, the square distance
%     to the closest point in y;
%   index is a 1-by-N vector indicating, for each point in x, the index
%     of the closest point in y.
%
% See also: SQRDIST.

% Copyright 2005 by Timothy E. Holy
