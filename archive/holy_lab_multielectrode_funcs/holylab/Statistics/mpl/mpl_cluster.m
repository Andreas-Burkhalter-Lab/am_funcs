function index_clust = mpl_cluster(xdens,options)
% MPL_CLUSTER: cluster data by flowing up the density
% Syntax:
%   index_clust = mpl_cluster(xdens,options)
% where
%   xdens is a d-by-N matrix containing one data point per column;
% and
%   index_clust is a 1-by-N vector containing the cluster number assigned
%     to each data point.
% You can control the behavior through the fields of options:
%   num_landmarks: This will determine how many points are used to
%     estimate the density.  (All of the points will contribute in some
%     fashion through the multiplicity assigned to the landmarks,
%     however.) See LANDMARK. For d > 1, the default value is 500; for 
%     d = 1, the default is Inf. Set to Inf to disable landmarking.
%   mode (default 'isotropic'): determines whether W is 'full', 'diag',
%     or 'isotropic'.
%   Kneighbors (default 2*d+1): sets the number of neighbors examined in
%     finding the most uphill point.
%   oversmooth: if present, is used as a scale factor on the optimal W to
%     force the density to be more smooth (if oversmooth > 1) than that
%     determined by cross-validation.
%
% You can also pass on any arguments to mpl_optw and mpl.  In particular,
% 'kernel_only' = 1 can (barely) speed up length calculations.
%
% See also: LANDMARK, MPL_OPTW, MPL.
  
% Copyright 2006 by Timothy E. Holy
  
  [d,N] = size(xdens);
  if (nargin < 2)
    options = struct;
  end
  if (d > 1 && N > 1000 && ~isfield(options,'num_landmarks'))
    options.num_landmarks = 500;
    warning('Landmarking. You can set options.num_landmarks = Inf to disable');
  end
  if ~isfield(options,'num_landmarks')
    options.num_landmarks = Inf;
  end
  if ~isfield(options,'mode')
    options.mode = 'isotropic';
  end
  if ~isfield(options,'Kneighbors')
    options.Kneighbors = 2*d+1;
  end
  
  % Chop down the size of the data set (for d > 1, this is O(N^2))
  [xlandmark,n_landmark,landmark_assignment] = ...
      landmark(xdens,options.num_landmarks);
  % Calculate the density at each landmark point
  [p,W] = mpl_optw(xlandmark,n_landmark,options);
  if isfield(options,'oversmooth')
    p = mpl(xlandmark,n_landmark,W*options.oversmooth,options);
  end
  % Find the most-uphill point, among nearest neighbors
  sd = sqrdist(xlandmark,xlandmark);
  [ssd,nbrrank] = sort(sd);
  [maxp,maxi] = max(p(nbrrank(1:options.Kneighbors,:)));
  clustmap = zeros(size(maxi));
  for i = 1:length(maxi)
    clustmap(i) = nbrrank(maxi(i));
  end
  % Flow each (landmark) point uphill, by iterating the cluster map
  clustmap0 = clustmap;
  clustmapc = clustmap;
  clustmap = clustmap(clustmap);
  while any(clustmap ~= clustmapc)
    clustmapc = clustmap;
    clustmap = clustmap(clustmap);
  end
  % Assign the cluster #s to the landmarks
  [cl,indx_tmp,index_clust] = unique(clustmap);
  % Assign the cluster #s to the original points
  index_clust = index_clust(landmark_assignment);

