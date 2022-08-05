% MSAMS_CONVERGE_MEX: move point(s) by adaptive meanshift, using landmarking
%
% This is the core algorithmic component of MSAMS.  This implementation
% is a much faster version of MSAMS_CONVERGE, with some differences
% in call syntax.
%
% Syntax:
%
%   outdata = msams_converge_mex(x,lminfo,y)
%   outdata = msams_converge_mex(x,lminfo,y,options)
%
% where
%   x is a d-by-N matrix of data points in d-dimensional space;
%   lminfo is a landmark structure of the type returned by
%     CHOOSE_LANDMARKS;
%   y is a d-by-q matrix of probe points in d-dimensional space;
%   options is a structure which may have the following fields:
%     min_to_check (default 3): the minimum # of points needed to satisfy the
%       MSAMS criterion
%     factor (default 3): the factor by which the mean shift has to
%       exceed the standard error to be considered significant
%     leapfrog (default true): use this to greatly increase the speed
%       (but not if you want to study the trajectories that the algorithm
%       produces)
%     backtrack (default true): if true, once the MSAMS criterion is
%       satisfied, the accepted movement (& diagnostic data) are taken
%       from the first time the mean shift exceeds the standard error
%       (i.e., as if factor = 1 had been supplied).  Making this false
%       will tend to over-aggregate clusters.
%     any_coordinate (default false): if true, the MSAMS criterion is
%       applied to each coordinate independently
%     convergence_thresh (default 1e-12): distance scale for
%       detection of cycling
%     max_iter (default 1000): maximum # of MSAMS steps that a given
%       probe point can undergo
%     n_threads (default # of cores): the # of threads to use; each
%       thread is assigned a different set of probe points.  It's not
%       useful to use more threads than there are processors on your
%       computer.  If you have a lot of processors, it might not even be
%       useful to use all of them, as bus contention for memory will
%       probably prevent this from scaling to really large numbers of
%       processors.
% and outdata is a structure which has the following fields:
%   yf is d-by-q matrix containing the "final" position of the
%     input probe points (before leapfrog convergence set in)
%   closestLandmark is the 1-by-q index of the landmark closest to each
%     output y
%   n is the 1-by-q vector containing the number of points
%     contributing to each probe point at their final positions
%   R2 is the squared radius for each final probe point
%   n_iter is the number of iterations required for convergence for
%     each probe point.  Note that because of cycling this may not be
%     the iteration that led to the final position.
%   convergedFlag (1-by-q) is a flag indicating whether each point
%     converged (0 = no convergence, 1 = converged, 2 = "converged" by
%     cycling)
%   settings is a copy of the options structure with defaults filled
%     in (checking this is the most reliable way to learn what the
%     defaults really are!)
%   ntraj (present only for q = 1) is the sequence of # of contributing
%     neighbors on each mean shift cycle
%   ytraj (present only for q = 1) is a d-by-niterations matrix that
%     stores the history of positions visited
%
% See also: MSAMS, CHOOSE_LANDMARKS, MSAMS_CONVERGE.

% Copyright 2007 by Timothy E. Holy

% Implemented as a MEX file for speed.
