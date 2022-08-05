% cn_mclinesearch: Monte Carlo T^2 simulation
% This function calculates the value of Hotelling's T^2 for
% randomly-generated point neighborhoods of different sizes.  Two different
% "flavors" of T^2 are computed, one from a fixed basepoint and the other at
% the maximum along a linesearch.
% 
% Syntax:
%   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim)
%   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim,npoints)
%   [T2,T2ls] = cn_mclinesearch(d,nlist,covarianceModel,nsim,npoints,nthreads)
% where
%   d is the dimensionality (an integer)
%   nlist is a vector of neighborhood sizes
%   covarianceModel is 'isotropic', 'diagonal', or 'full'
%   nsim is the number of simulations to run
%   npoints is the total number of random points generated
%     (default 2*d*max(nlist))
%   nthreads is the number of computational threads to use in the MEX
%     file (default ncores-1, up to a maximum of 8)
% and
%   T2 is a matrix of size length(nlist)-by-nsim, containing the value of
%     T^2 using a fixed basepoint for the different-sized neighborhoods
%   T2ls is similar, except that it represents the maximum T^2 for
%     displacements along the line from the origin to the neighborhood
%     mean subject to the constraint that the n points in the
%     neighborhood are the closest n points.
%
% This function can take a while to execute, and because it is a MEX file
% it cannot be interrupted.  To get a feeling for how long a full run
% takes, first time the execution with a relatively small value of nsim
% (say, under 1000) and then extrapolate linearly.

% Copyright 2011 by Timothy E. Holy

% Implemented as a MEX file for speed