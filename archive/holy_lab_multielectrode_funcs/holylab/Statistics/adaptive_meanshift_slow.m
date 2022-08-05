function [comp,y] = adaptive_meanshift(x,y,options)
% ADAPTIVE_MEANSHIFT: perform locally-adaptive meanshift (LAMS)
% Syntax:
%   clust = adaptive_meanshift(x,y)
%   [clust,yf] = adaptive_meanshift(x,y,options)
% where
%   x is a d-by-N matrix of data points in d-dimensional space;
%   y is a d-by-nlandmarks matrix of landmark points in d-dimensional
%     space;
% and
%   comp is a cell array of length nclusters, each entry containing the
%     indices of the landmarks belonging to a particular cluster;
%   yf is a d-by-nlandmarks matrix, with each column containing the
%     final resting point of each landmark.
%
% See also: MEANSHIFT, RMEANS, ADAPTIVE_MEANSHIFT1.
  
% Old output:
%   clust is a 1-by-nlandmarks vector, giving the cluster number
%     assigned to each landmark.

% Copyright 2006 by Timothy E. Holy
  
  [d,N] = size(x);
  [dy,q] = size(y);
  if (d ~= dy)
    error(['The dimensionality of the data and the landmarks are not ' ...
           'consistent']);
  end
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'thresh')
    options.thresh = 1e-6;
  end
  if ~isfield(options,'mergetol')
    options.mergetol = 0.1;
  end
  if ~isfield(options,'maxiter')
    options.maxiter = 100;
  end
  if ~isfield(options,'showprogress')
    options.showprogress = false;
  end
  if ~isfield(options,'ploteach')
    options.ploteach = false;
  end
  if ~isfield(options,'plotproject')
    dlim = min(2,d);   % Plotting in 1 or 2 dimensions
    options.plotproject = zeros(dlim,d);
    options.plotproject(1:dlim,1:dlim) = eye(dlim);
  end
  if ~isfield(options,'plotpause')
    options.plotpause = false;
  end
  
  ismoving = true(1,q);
  R2 = zeros(1,q);
  yOld = y;
  iter = 0;
  while (any(ismoving) && iter < options.maxiter)
    %for i = find(ismoving)
    %  [y(:,i),R2(i)] = adaptive_meanshift1(x,y(:,i),options);
    %end
    % NOTE: to use adaptive_meanshift_step(), comment out the 3 lines above
    %       and use next line:
    [y(:, ismoving), R2(ismoving)]=adaptive_meanshift_step(x,y(:,ismoving));
    
    % See whether the current position has been visited before (i.e., we're
    % on a cycle).
    ydiff = repmat(y,[1 1 size(yOld,3)]) - yOld;
    sd = sum(ydiff.^2,1);
    minsd = squeeze(min(sd,[],3));  % Distance to closest predecessor
    ismoving = minsd > options.thresh*R2;
    yOld(:,:,end+1) = y;
    iter = iter+1;
    if options.showprogress
      fprintf('  %d\n',sum(ismoving));
    end
    if options.ploteach
      meanshiftplot(x,y,sqrt(R2),options);
    end
  end
  % Now combine close landmarks---those that include each other's centers
  % within their radius of influence
  connectivity_graph = false(q,q);
  sd = sqrdist(y,y)/options.mergetol^2;
  for i = 1:q
    connectivity_graph(i,i) = true;
    for j = i+1:q
      %if (sd(i,j) < R2(i) && sd(i,j) < R2(j))
      if (sd(i,j) < R2(i) || sd(i,j) < R2(j))
        connectivity_graph(i,j) = true;
        connectivity_graph(j,i) = true;
      end
    end
  end
  comp = connected_components(connectivity_graph);
  %yclump = crossmoveR(y,y,0.5*sqrt(R2),0);
  %[yu,yind,clust] = unique(yclump','rows');
  %clust = clust';
  
  
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
    angle = linspace(0,2*pi,200);
    for i = 1:size(yp,2)
      xc = yp(1,i) + R(i)*cos(angle)*sqrt(dim/odim);
      yc = yp(2,i) + R(i)*sin(angle)*sqrt(dim/odim);
      plot(xc,yc,'r')
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
  