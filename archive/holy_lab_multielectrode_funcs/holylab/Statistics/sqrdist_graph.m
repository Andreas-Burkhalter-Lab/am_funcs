function d2 = sqrdist_graph(D2,nodeList,w)
% SQRDIST_GRAPH: square-distances to an off-graph point
% 
% In a graph, one way to define an off-graph point is as a weighted sum of
% a set of on-graph points.  This function computes the square-distance of
% all nodes to a single off-graph point.
%
% Syntax:
%   d2 = sqrdist_graph(D2,nodeList)
%   d2 = sqrdist_graph(D2,nodeList,w)
% where
%   D2 is the N-by-N matrix of distances
%   nodeList is the list of n nodes making up the "neighborhood" defining
%     the off-graph point
%   w is the list of n weights, in the order defined by nodeList. If
%     absent, all points in the list are weighted equally. The weights
%     should typically sum to 1, but this is not checked.
% and
%   d2 is an N-by-1 vector containing the effective square-distances to
%     the point. 
  
% Copyright 2010 by Timothy E. Holy
      
  m = length(nodeList);
  if (m == 1)
    d2 = D2(:,nodeList);
  else
    D2nbrhd = D2(nodeList,nodeList);
    if (nargin < 3 || isempty(w))
      % Unweighted case
      D2nbrhd_sum = sum(D2nbrhd(:));
      d2 = sum(D2(:,nodeList),2)/m - D2nbrhd_sum/(2*m^2);
    else
      % Weighted case
      wwD2n = bsxfun(@times,w(:),bsxfun(@times,w(:)',D2nbrhd));
      D2nbrhd_sum = sum(wwD2n(:));
      d2 = sum(bsxfun(@times,w(:),D2(:,nodeList)),2) - D2nbrhd_sum/2;
    end
  end
