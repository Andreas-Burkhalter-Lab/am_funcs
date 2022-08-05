function out = cluster_compare(clust1,clust2,method)
% CLUSTER_COMPARE: compare two clusterings
%
% Given a data set with N points, a clustering corresponds to an
% assignment of a cluster index to each point.  Two different cluster
% assignments might permute the cluster labels but be otherwise
% identical, or there might be some disagreements.  
%  
% This function allows you to use any of several different algorithms (or
% all of them, if you don't specify a method).  The methods are:
%    'mi': mutual information
%       Nonnegative. 0 means no relationship between cluster
%       assignments, but there can be sample-size effects.
%    'mis','misN': mutual information with shuffle correction. The shuffle
%      correction accounts for sample-size effects given the number of
%      observations and points.  'mis' by itself does a single shuffle;
%      'misN', where N is an integer, does N shuffles.
%    'crand': corrected Rand (Hubbert, L. J., Arabie, P. (1985),
%       Comparing partitions, Journal of Classiffication, 2:63-76.)
%       It varies between -1 and 1, with values above 0 indicating
%       increasing agreement (and 1 indicating perfect agreement).
%
% If you want to rely on a single measure, 'mis' is recommended.
%
% Syntax:
%   indx = cluster_compare(clust1,clust2,method)
%   results = cluster_compare(clust1,clust2)
% Inputs:
%   clust1 and clust2 are the two vectors of cluster assignments;
%   method is a string containing the name of the method, or a cell array
%     containing a list of methods. If missing, all methods will be used.
% Output:
%   indx is the selected measure of similarity if only one was selected, OR
%   results is a structure containing the results of multiple measures.
  
% Copyright 2008 by Timothy E. Holy

  clust1 = clust1(:);
  clust2 = clust2(:);
% Conceptually we want to do the following:
%   for i = 1:n1
%     for j = 1:n2
%       count(i,j) = sum((clust1 == i) & (clust2 == j));
%     end
%   end
% but:
%   1. We don't want to assume that the cluster labels go from 1...n
%   2. If the number of clusters is big, we want to do this in an efficient
%      way.
% So, instead do the following:
  [urows,tmp,pairIndex] = unique([clust1 clust2],'rows');
  [u1,tmp,indx1] = unique(urows(:,1));
  [u2,tmp,indx2] = unique(urows(:,2));
  n1 = length(u1);
  n2 = length(u2);
  count = zeros(n1,n2);
  [pairIndexC,pairIndexN] = agglabel(pairIndex);
  
  countIndex = sub2ind([n1 n2],indx1,indx2);
  count(countIndex) = pairIndexN;  % there!
  
  % Do the marginals
  count1 = sum(count,2);
  count2 = sum(count,1);
  N = sum(count1);
  
  pjoint = count / sum(count(:));
  p1 = count1/N;
  p2 = count2/N;
  
  methods = {'crand','mi','mis'};
  if (nargin < 3)
    method = methods;
  end
  if ischar(method)
    method = {method};
  end
%   if ~isempty(setdiff(method,methods))
%     error('Some methods not recognized');
%   end
  
  if ~isempty(strmatch('mi',method,'exact'))
    result.mi = entropy(p1) + entropy(p2) - entropy(pairIndexN/N);
  end
  misIndex = strmatch('mis',method);
  if ~isempty(misIndex)
    Nshuf = sscanf(method{misIndex},'mis%d',1);
    if isempty(Nshuf)
      Nshuf = 1;
    end
    Hj = zeros(1,Nshuf);
    for i = 1:Nshuf
      [urowsS,tmp,pairIndexS] = unique([clust1(randperm(N)) clust2],'rows');
      [tmp,pairIndexNS] = agglabel(pairIndexS);
      Hj(i) = entropy(pairIndexNS/N);
    end
    result.(method{misIndex}) = mean(Hj) - entropy(pairIndexN/N);% + entropy(p1) + entropy(p2);
  end
  if ~isempty(strmatch('crand',method,'exact'))
    count_choose2 = count .* (count-1)/2;
    scount1_choose2 = sum(count1 .* (count1-1)/2);
    scount2_choose2 = sum(count2 .* (count2-1)/2);
    N_choose = N*(N-1)/2;
    corrfac = scount1_choose2*scount2_choose2/N_choose;
    
    result.crand = (sum(count_choose2(:)) - corrfac) / ...
	   ((scount1_choose2 + scount2_choose2)/2 - corrfac);
  end
  
  if (length(method) == 1)
    out = result.(method{1});
  else
    out = result;
  end
  
  
  function H = entropy(p)
    p = p(p > 0);
    H = - sum(p .* log2(p));
    