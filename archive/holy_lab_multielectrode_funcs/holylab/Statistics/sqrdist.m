% SQRDIST: Find Euclidean square distances between points
% Usage:
%  sd = sqrdist(x,y)
% where
%   x is a d-by-N matrix of points in d-dimensional space;
%   y is a d-by-q matrix of points in d-dimensional space;
% and
%   sd is a N-by-q matrix containing the square distances between each pair
%     of x,y, i.e., sqrt(d(i,j)) is the distance between x(:,i) and y(:,j).
%
% This function is multithreaded, with different points in y handled by
% different CPUs. You can thus decrease execution time by making sure that
% y contains a reasonable number of points; it is inefficient to call this
% with N large and q=1 (better N=1 and q large).
%
% NOTE: If you are using old code and it returns an error on this function
% it may have been written using an older version of sqrdist (which was
% written in matlab code and only accepted one input); that version can
% still be found off the matlab path in /usr/lab/etc/old_code.
%
%  See also: MINDIST.