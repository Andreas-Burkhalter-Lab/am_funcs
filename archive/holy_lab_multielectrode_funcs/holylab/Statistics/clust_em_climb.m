function landmarkClust = clust_em_climb(x,y,Rfactor,options)
% CLUST_EM_CLIMB: cluster data based on an EM estimate of the density,
% and the climbing of landmarks up to peaks of the density
% Syntax:
%   landmarkClust = clust_em_climb(x,y,Rfactor,options)
% where
%   x is a d-by-N data matrix;
%   y is a d-by-q matrix containing landmark positions;
%   Rfactor is the scaling factor applied to the radius of each
%     gaussian in climbing uphill;
%   options is a structure which may have the following fields:
%     ploteach (default false): if true, the data & centers are plotted
%       on each iteration;
%     plotpause (default false): if true, the algorithm pauses for
%       keyboard input after plotting;
%     plotproject: if supplied, projects the data & centers down to two
%       dimensions before plotting (a 2-by-d matrix) (default: use first
%       two dimensions of the input);
% and
%   landmarkClust is a cluster number assigned to each landmark.
%
% See also: MEANSHIFT.
  
% Copyright 2005 by Timothy E. Holy

  % Save variables from last call, in case we can skip the step in which
  % we build a gaussian model
  persistent lastx lasty alpha mu sigma2

  % Do we also want to re-compute landmarkIndex in terms of maximum
  % likelihood rather than nearest point?
  if (nargin < 4)
    options = struct;
  end
  options.hard = 0;
  options.movemu = 0;
  if (~isequal(x,lastx) || ~isequal(y,lasty))
    % Build a gaussian model of the data
    [alpha,mu,sigma2] = em_gauss(x,y,options);
  end
  lastx = x;
  lasty = y;
  % Now aggregate the landmarks by flowing them uphill
  landmarkClust = climb_gaussian(y,alpha,mu,Rfactor^2*sigma2,options);
  %landmarkClust = nzIndex(landmarkClust); % in case any landmarks disappeared
  