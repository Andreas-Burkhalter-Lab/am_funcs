function [Xu,nu] = unique_data(X,n)
% UNIQUE_DATA: reduce to a set of unique data points & multiplicities
% Syntax:
%   [Xu,nu] = unique_data(X)
%   [Xu,nu] = unique_data(X,n)
% where
%   X is a N-by-d matrix, with each data point on one ROW of X;
%   n (optional) is a 1-by-N vector containing the multiplicity
%     associated with each input data point (default: all 1)
% and
%   Xu is the set of unique rows of X;
%   nu is a vector, where nu(i) is the number of times Xu(i,:) appears.

% Copyright ?-2010 by Timothy E. Holy
  
  N = size(X,1);
  if (nargin < 2 || isempty(n))
    n = ones(1,N);
  end
  [Xu,~,Xindex] = unique(X,'rows');
  M = size(Xu,1);
  nu = zeros(1,M);
  for i = 1:N
    nu(Xindex(i)) = nu(Xindex(i))+n(i);
  end
  