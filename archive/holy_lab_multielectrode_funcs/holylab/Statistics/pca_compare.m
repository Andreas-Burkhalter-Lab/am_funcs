function sv = pca_compare(components,options)
% PCA_COMPARE: calculate the degree of overlap between projections
% 
% Syntax:
%   sv = pca_compare(components)
% where
%   components is a d-by-p-by-npca matrices: components(:,:,i) is the ith
%     set of p principal components, for data in d dimensions.
%
% This function calculates the redundancy among component sets in the
% following way: for each pair of sets, the union of both sets is
% formed.  An SVD is performed to determine the singular values; on
% output, sv is a vector (of length 2p) containing the singular values,
% averaged across all pairs.
%
% You can estimate the number of consistent dimensions as the number of
% "small" elements of sv, or more robustly as 2p - sum(sv.^2).
%
% This naturally assumes that d > 2p.
%
% See also: ?.
  
% Copyright 2005 by Timothy E. Holy
  
  [d,p,npca] = size(components);
  dotprod = zeros(1,p);
  for i = 1:npca
    for j = i+1:npca
      dotprod = dotprod + squeeze(abs(sum(components(:,:,i).*components(:,:,j),1)));
    end
  end
  sv = dotprod/(npca*(npca-1)/2);
  return;
  
  [d,p,npca] = size(components);
  if (d < 2*p)
    error('This makes no sense with d < 2p');
  end
  sv = zeros(1,2*p);
  for i = 1:npca
    for j = i+1:npca
      tComp = [components(:,:,i) components(:,:,j)];
      [U,S,V] = svd(tComp,0);
      sv = sv + diag(S)';
    end
  end
  sv = sv/(npca*(npca-1)/2);
