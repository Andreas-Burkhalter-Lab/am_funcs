% adaptive_meanshift_step: take one step of an adaptive meanshift
%
% Syntax:
%   [yout,R2,n] = adaptive_meanshift_step(x,y,factor,minN)
% where
%   x is a d-by-N matrix of data points in d-dimensional space (must be of
%     class double);
%   y is a d-by-q matrix of landmarks that will move up the density
%      defined by x (must be of class double);
%   factor (default 3) is the number of times that the step size needs
%      to exceed the standard error by to be considered significant;
%   minN (default 10) is the minimum number of neighbors to use to
%      initially judge whether a step should be taken.
% and
%   yout is the new set of landmark positions;
%   R2 is a 1-by-q vector containing the set of square radii used
%     for each landmark;
%   n is a 1-by-q vector containing the number of data points (x)
%     within the radius at each landmark.
%
% See also: ADAPTIVE_MEANSHIFT.

% Implemented as a MEX file for speed.

% Copyright 2006 by Timothy E. Holy
