function [map,leapfrog_info] = leapfrog(map,n)
% LEAPFROG: rapid convergence for mean shift algorithms by graph iteration
% Syntax:
%   mapout = leapfrog(mapin)
% where map defines a directed graph, i->map(i).
%   mapin is the input map;
%   mapout is the output map, where each point maps to its "final"
%     location.
% The notion of a "final" location can run into trouble when the input
% graph contains cycles.  In that case, cycles need to be broken, so that
% all connected points map to the same location.  The syntax
%
%   mapout = leapfrog(mapin,priority)
%
% with the optional input "priority", a vector with one entry per map
% point, allows one to control how cycles are broken: points with the
% highest priority will win.  In cases of a tie, the point with the
% smallest index wins.  The default priority is 1 for all points.
%
%
%  [map,leapfrog_info] = leapfrog(mapin,...)
%
% allows you to obtain some information about the performance of leapfrog.
%
% See also: MSAMS, MEANSHIFT.

% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 2)
    % n can be used so that consolidation occurs at the peak. But if it's
    % not available, set a default so that it will simply go to the
    % lowest index.
    n = ones(size(map));
  end
  if (nargout > 1)
    leapfrog_info = struct('n_leapfrog_steps',0,'n_repair_steps',0);
  end
  sz = size(map);
  map = map(:)';
  n = n(:)';

  % Every landmark "maps" to some landmark, either itself (usually
  % rare) or to another landmark.
  % We can accelerate convergence by iterating this map, so that points
  % leapfrog uphill.
  map0 = map;
  umap = unique(map);
  umapOld = [];
  while ~isequal(umapOld,umap)
    umapOld = umap;
    map = map(map);
    umap = unique(map);
    if (nargout > 1)
      leapfrog_info.n_leapfrog_steps = leapfrog_info.n_leapfrog_steps+1;
    end
  end
  % The leapfrogging process can cut off some points: consider a map of
  % only three points, where the mapping is
  %    2 3 2
  % i.e., 1->2, 2->3, 3->2. Note there is a cycle. If we iterate this map
  % we get
  %    3 2 3
  % which is stable under further iterations with itself. However, it looks
  % like there are two "clusters," [1 3] and [2]. But in fact they were all
  % connected originally. Consolidate these splits due to cycles by
  % checking the original map again, and moving each point to our best
  % guess for the most uphill point of two possible mappings. (We could
  % just choose, say, the lower index of two possible mappings, but in fact
  % the map is interesting in its own right, and we might as well try to
  % preserve as much of its character as possible.) 
  % This usually converges very quickly, although there may not be
  % guarantees that it will do so.  This could be replaced by a 'real'
  % (well-optimized) graph theory approach, computing the reduced graph of
  % vertices and then calculating the number of components by standard
  % algorithms. But this converges so fast, by comparison with the other
  % steps, that it seems that a low-tech approach is fine.
  is_changed = true;
  iter = 0;
  maxiter = length(umap);
  while (any(is_changed) && iter < maxiter)
    map_compare = [map; map(map0)];
    is_changed = (map_compare(1,:) ~= map_compare(2,:));
    map_compare = map_compare(:,is_changed);
    n1 = n(map_compare(1,:));  % estimate max density from landmark that ...
    n2 = n(map_compare(2,:));  % ... averages over the most data points
    [tmp,index_max] = max([n1; n2]);
    is_tie = (n1 == n2);
    if any(is_tie)
      % Ties can result it an infinte loop. Break them by choosing the one
      % with the lower index
      [tmp,index_max(is_tie)] = min(map_compare(:,is_tie));
    end
    map_uphill_index = sub2ind([2 length(n1)],index_max,1:length(n1));
    map(is_changed) = map_compare(map_uphill_index);
    iter = iter+1;
  end
  if (iter > maxiter)
    error('Map repair failed to be fast')
  end
  if (nargout > 1)
    leapfrog_info.n_repair_steps = iter;
  end
  map = reshape(map,sz);
  