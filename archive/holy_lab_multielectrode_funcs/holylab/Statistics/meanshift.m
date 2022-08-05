function clust = meanshift(x,y,options)
% MEANSHIFT: Use mean-shift algorithm to move landmarks
% This is a modified mean-shift algorithm, in that it moves a
% set of q landmark points "uphill," up the density of the supplied
% points.  With each iteration of the algorithm, a fraction of the input
% points are discarded.  This procedure seems to "jostle" the
% landmarks out of small local maxima in the density.
%
% Syntax:
%   clust = meanshift(x,y,R)
%   clust = meanshift(x,y,options) % here, R must be set as options.R
% where
%   x is a d-by-N matrix of data points in d-dimensional space;
%   y is a d-by-nlandmarks matrix of landmark points in d-dimensional
%     space;
%   R is the mean-shift radius, either a scalar or vector (one for each
%     landmark or data point, if length(R) = N, then they are attached to
%     the data points);
%   options is an options structure with the following fields:
%     shrinkfrac (default 0.3): percentage of points in x to throw out
%       with each iteration
%     usey (default true): if true, y_i is thrown into the crossmoveR step.
%     ploteach (optional): if true, will show the status of the
%       clustering at each iteration, using the first two dimensions of
%       x, or the first dimension if x is one-dimensional.
%       (but see "plotproject") If it's an axis handle, it will use this
%       axis to do the plot.
%     plotproject (optional): if supplied, will do visualization in
%       terms of x coordinates projected on these vectors. Should be a
%       2-by-d matrix.
%     plotpause (optional): if true, will pause after each iteration and
%       wait for the user to hit a key. This only has effect if ploteach is
%       true. Default: false.
%  and
%    clust is a 1-by-nlandmarks vector, giving the cluster number
%      assigned to each landmark.
%
% HISTORY:
%
%   2006-09-20 (RCH)    Changed calling syntax to be consistent with syntax
%                       used for adaptive meanshift syntax, so that calls
%                       that occur when cluster_func = @meanshift still
%                       work
%
% See also: RMEANS, CLUST_EM_CLIMB.
  
% Copyright 2004 Timothy E. Holy
  
  if isnumeric(options)
    R = options;
    options = struct;
  elseif isstruct(options)
    if isfield(options,'R')
      R = options.R;
    else
      R = options.Rfactor; % TEH: why is this Rfactor?? That's a scale factor, not an absolute distance
    end
  end
  options = meanshiftoptions(options);
  [d,N] = size(x);
  if ~isfield(options,'plotproject')
    dlim = min(2,d);   % Plotting in 1 or 2 dimensions
    options.plotproject = zeros(dlim,d);
    options.plotproject(1:dlim,1:dlim) = eye(dlim);
  end
  if (length(R) == N)
    options.R_on_data = 1;
  end
  if (options.R_on_data && options.overlap)
    error('Can''t do overlap mode when R is assigned to data points');
  end

  q = size(y,2);
  if (options.ploteach)
    if ishandle(options.ploteach)
      axes(options.ploteach)
    end
    meanshiftplot(x,y,R,options)
  end
  yold = y;
  if options.R_on_data
    y = crossmoveRx(x,y,R);
  else
    y = crossmoveR(x,y,R,options.usey);
  end
  isdone = isequal(y,yold);
  if options.shrinkfrac
    isdone = size(x,2) < q;
  end
  while ~isdone
    % Plot
    if options.ploteach
      if ishandle(options.ploteach)
        axes(options.ploteach)
      end
      meanshiftplot(x,y,R,options)
    end
    % Toss points
    npts = size(x,2);
    ntoss = round(options.shrinkfrac * npts);
    tossindx = round((npts-1)*rand(1,ntoss))+1;
    x(:,tossindx) = []; % because of duplicates, will be less than ntoss, but who cares?
    % Update y
    yold = y;
    if options.R_on_data
      y = crossmoveRx(x,y,R);
    else
      y = crossmoveR(x,y,R,options.usey);
    end
    %isdone = isequal(y,yold);
    dy = y-yold;
    isdone = all(sum(dy.^2,1) < 1e-8*R.^2);
    if options.shrinkfrac
      isdone = size(x,2) < q;
    end
  end
  % Now do a single pass where we throw out all the remaining x points, just
  % using the ys, to make very close points aggregate
  if ~options.overlap
    y = crossmoveR(y,y,median(R),0);  % use scalar R and ~usey, so that all
                                   % y will have the same view of the
                                   % world; if y_i can see y_j, so too can
                                   % y_j see y_i
    % Identify the points belonging to a common cluster
    [yu,yind,clust] = unique(y','rows');
    clust = clust';
  else
    % Determine whether circles overlap; any that overlap are in the same
    % cluster
    if isscalar(R)
      R = repmat(R,1,q);
    end
    ydist = sqrt(sqrdist(y'));
    sumR = repmat(R,q,1) + repmat(R',1,q);
    isconnected = ydist < sumR;
    comp = connected_components(isconnected);
    ncomp = length(comp);
    clust = zeros(1,q);
    for i = 1:ncomp
      clust(comp{i}) = i;
    end
  end
    
  
  
function op = meanshiftoptions(op)
  if ~isfield(op,'ploteach')
    op.ploteach = 0;
  end
  if ~isfield(op,'shrinkfrac')
    op.shrinkfrac = 0.2;
  end
  if ~isfield(op,'plotpause')
    op.plotpause = 0;
  end
  if ~isfield(op,'usey')
    op.usey = 1;
  end
  if ~isfield(op,'R_on_data')
    op.R_on_data = 0;
  end
  if ~isfield(op,'overlap')
    op.overlap = 0;
  end
  


  
function meanshiftplot(x,y,R,options)
  if (isscalar(R) & size(y,2) > 1)
    R = repmat(R,1,size(y,2));
  end
  xp = options.plotproject*x;
  yp = options.plotproject*y;
  dim = size(xp,1);  % # of plot dimensions
  odim = size(x,1);  % # of original dimensions
  if (dim > 1)
    plot(xp(1,:),xp(2,:),'.')
    hold on
    if options.R_on_data
      scatter(yp(1,:),yp(2,:),'ro')
    else
      angle = linspace(0,2*pi,200);
      for i = 1:size(yp,2)
        xc = yp(1,i) + R(i)*cos(angle)*sqrt(dim/odim);
        yc = yp(2,i) + R(i)*sin(angle)*sqrt(dim/odim);
        plot(xc,yc,'r')
      end
    end
    hold off
    axis equal
  else
    hist(xp,ceil(sqrt(size(x,2))));
    hold on
    scatter(yp,zeros(1,size(y,2)),'ro');
    hold off
  end
  if isfield(options,'plottitle')
    title(options.plottitle);
  end
  if options.plotpause
    pause
  else
    drawnow
  end

  
  