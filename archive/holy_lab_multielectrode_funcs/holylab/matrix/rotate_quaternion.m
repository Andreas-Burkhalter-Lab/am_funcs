function Xr = rotate_quaternion(X,q)
% rotate_quaternion: perform 3d rotations around origin, quickly
%
% Quaternions provide an efficient mechanism for performing 3d
% rotations.  See
%  http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
% 
% Syntax:
%   Xr = rotate_quaternion(X,q)
% where
%   X is a 3-by-N matrix of points in three-dimensional space
%   q is a 4-vector containing the coordinates of the quaternion;
%     if xhat is the unit vector of the axis of rotation, and alpha the
%     angle of rotation, then
%       q = [cos(alpha/2) sin(alpha/2)*xhat]
% and
%   Xr contains the rotated points.
  
% Copyright 2010 by Timothy E. Holy
  