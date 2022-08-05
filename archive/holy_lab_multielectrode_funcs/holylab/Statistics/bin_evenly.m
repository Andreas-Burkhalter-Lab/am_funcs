function [index,c] = bin_evenly(x,varargin)
  % bin_evenly: split points into bins with similar numbers of points
  % Syntax:
  %    [index,c] = bin_evenly(x,k)
  % Splits points x, a d-by-N matrix with each point corresponding to a
  % column, into k bins. index is the assignment of each point, and c is a
  % d-by-k matrix containing the centers. The final assignment is always to
  % the closest center, even though the algorithm internally uses a
  % different intermediate strategy. The initial center positions are drawn
  % randomly from x.
  %
  %    [index,c] = bin_evenly(x,c0)
  % Specify an initial guess for c, rather than just the number k of centers.
  %
  %    [index,c] = bin_evenly(x,w,...)
  % Allow each point to have a "weight" (default 1). In this case, the
  % total weight per bin is approximately even.
  %
  %    [index,c] = bin_evenly(...,options)
  % Manually specify certain parameters as fields of the structure
  % "options":
  %    thresh: the threshold for the minimum weight to qualify as
  %      "even binning." The default is Wtot/k - sqrt(Wtot/k), where Wtot
  %      is the total weight (or total # of points when all weights are 1).
  %      This choice effectively corresponds to "equal within sampling
  %      variability."
  %    iter_max (default inf): constrains the maximum number of iterations
  %      of the algorithm. When set to inf, the algorithm proceeds until
  %      the total cost stops decreasing.
  % 
  % See also: kmeans_hard, likelihood_ratio_test.
  
  % Copyright 2012 by Timothy E. Holy
  
  N = size(x,2);
  carg = 1;
  if isequal(size(varargin{carg}),[1,N])
    w = varargin{carg};
    carg = carg+1;
  else
    w = ones(1,N);
  end
  if isscalar(varargin{carg})
    % k syntax
    k = varargin{carg};
    rp = randperm(N);
    c = x(:,rp(1:k));
  else
    % "improve" syntax
    c = varargin{carg};
    k = size(c,2);
  end
  Wtot = sum(w);
  carg = carg+1;
  options = struct;
  if (length(varargin) >= carg)
    options = varargin{carg};
  end
  options = default(options,'thresh',Wtot/k-sqrt(Wtot/k),'iter_max',inf);
  
  if (options.thresh > Wtot/k)
    error('No solution possible; your threshold is too high.')
  end
 
  % Iterate until the cost stops decreasing or we exceed the iteration
  % count
  iter = 0;
  costOld = inf;
  while iter < options.iter_max
    d2 = sqrdist(c,x);
    % For each center, assign points by increasing distance until the
    % weight threshold is crossed
    [~,sortOrder] = sort(d2,2);
    assignment = zeros(k,N);
    curw = zeros(k,1);  % The weight assigned to each centroid so far
    wgain = w; % the potential increase in weight for assigning each data point
    while any(curw < options.thresh)
      cumw = cumsum(wgain(sortOrder),2);  % the cumulative weight as each neighborhood grows
      for i = 1:k
        if (curw(i) < options.thresh)
          % Add more points to this center
          indx = find(cumw(i,:)>=options.thresh,1,'first');
          assignment(i,sortOrder(i,1:indx)) = 1;
        end
      end
      sa = sum(assignment,1);
      wnorm = zeros(1,N);
      wnorm(sa>0) = w(sa>0) ./ sa(sa>0);
      wgain = w ./ (sa+1);
      % Of the next two lines:
      % The first forces the normalized weight to exceed threshold;
      % The second avoids this requirement, only the "naive" weight
      %   must exceed threshold.
      % The second is recommended, both because it is much faster and
      % because it seems to progress more smoothly to a near-minimum.
      curw = sum(bsxfun(@times,wnorm,assignment),2); %#ok<NASGU>
      curw = inf;  % this replaces the previous value and will terminate the while loop
    end
    % All unclaimed points are assigned to their closest center
    flag = sa == 0;
    [~,indexClosest] = min(d2(:,flag),[],1);
    indx = sub2ind(size(assignment),indexClosest,find(flag));
    assignment(indx) = 1;
    assignment = bsxfun(@rdivide,assignment,sum(assignment,1));
    % Compute the cost and check for termination
    wassignment = bsxfun(@times,w,assignment);
    cost = sum(wassignment(:).*d2(:));
    fprintf('cost %g\n',cost);
    if (cost >= costOld)
      break
    end
    costOld = cost;
    % Shift the centers to the weighted mean
    wnorm = sum(wassignment,2);
    for i = 1:k
      c(:,i) = sum(bsxfun(@times,wassignment(i,:),x),2)/wnorm(i);
    end
    iter = iter+1;
  end
  % Finally, assign each point to its closest center
  [~,index] = min(d2,[],1);
end
