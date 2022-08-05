function npb = binpairs(t,trange)
% BINPAIRS: Calculate number of points falling between pairs of values
% Syntax:
%   npb = binpairs(t,trange)
% where
%   t is a vector containing the data points (does not have to be sorted);
%   trange is an n-by-2 matrix, where each row contains the start and
%     stop for each bin;
% and
%   npb is the number of points in each bin.
% 
% See also: HISTSPLIT.
  
% Copyright 2005 by Timothy E. Holy
  lt = length(t);
  np = size(trange,1);
  % The algorithm: calculate where each boundary would fall in the sorted
  % data vector. Then, compute the number of points between pairs as the
  % difference between the positions of the end boundaries and beginning
  % boundaries.
  for i = 1:2
    % Intercalate boundaries into points
    [ts,ii] = sort([t(:)',trange(:,i)']);
    % Find the location of the boundaries in sorted vector
    indx = find(ii > lt);
    % Calculate their location without the boundary points themselves
    tmp = indx - (0:np-1);
    % Calculate the order in which they intercalated
    p = ii(indx)-lt;
    % Invert this permutation
    [ps,ip] = sort(p);
    % Put them in the right order
    pos(i,:) = tmp(ip);
  end
  % Perhaps the boundaries are not in increasing order, so we have to
  % throw out negative values
  npb = max([diff(pos,1,1); zeros(1,size(pos,2))]);
