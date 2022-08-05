function Y = maps2sim(maps)
% MAPS2SIM: convert a set of MSAMS-maps to a similarity matrix
% Syntax:
%   Y = maps2sim(maps)
% where
%   maps is a n-by-q array of maps (e.g., from msams_resample)
% and
%   Y is a similarity matrix. Y(i,j) describes the similarity between
%   probe points i and j, where "1" indicates that two landmarks are
%   always linked, "0" indicates that they are never linked, and
%   intermediate values are possible.
%
% 1-Y is a good candidate for hierarchical clustering.
%
% See also: MSAMS_RESAMPLE, LINKAGE, SIM2CLUST.

% Copyright 2008 by Timothy E. Holy

[n,n_lm] = size(maps);

Y = zeros(n_lm,n_lm);
for i = 1:n
  [ul,tmp,clust] = unique(maps(i,:));
  clustc = agglabel(clust);
  for j = 1:length(clustc)
    Y(clustc{j},clustc{j}) = Y(clustc{j},clustc{j}) + 1;
  end
end
Y = Y/n;