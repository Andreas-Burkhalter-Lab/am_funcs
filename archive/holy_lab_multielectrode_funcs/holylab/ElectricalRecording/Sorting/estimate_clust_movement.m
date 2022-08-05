function ell = estimate_clust_movement(clust,clustsize,nsnipsperblock)
% ESTIMATE_CLUST_MOVEMENT: estimate length scale of cluster movement
% between blocks of spikes
%
% The length scale of cluster movement between two adjacent blocks is
% estimated as
%     mean(r_i * |f_i,1 - f_i,2|/(f_i,1 + f_i,2))
%          (mean is across i = landmark index)
% where f_i,1 is the fraction of total spikes associated with landmark i
% in block 1, similarly for f_i,2.
% 
% Syntax:
%   ell = estimate_clust_movement(clust,clustsize,nsnipsperblock)
% where
%   clust is a vector of "cluster" numbers (e.g., assigned by kmeans)
%     assigned to each snippet;
%   clustsize is a vector of length nclusters, giving a measure of the size
%     of each cluster (e.g., see the 3rd output of kmeans_hard);
%   nsnipsperblock is the number of snippets in each "block" to be
%     compared; sum(nsnipsperblock) = length(clust)
% and
%   ell is a measure of the distance moved between adjacent blocks.
%
% See also: AUTOSORT_CALC_T2V.

  nblocks = length(nsnipsperblock);
  nsnipscum = [0 cumsum(nsnipsperblock)];  % Lookup table for breaks between blocks
  %[clust,c,rmsd] = kmeans_hard(cat(2,snip{:}),options.k);
  % Find number in each cluster, broken down by blocks
  uclust = unique(clust);
  nclusts = length(uclust);
  nperclust = zeros(nblocks,nclusts);
  for i = 1:nblocks
    clabel = agglabel(clust(nsnipscum(i)+1:nsnipscum(i+1)));
    for j = 1:length(clabel)
      nperclust(i,j) = length(clabel{j});
    end
  end
  % Length scale of movement between two blocks:
  %     mean(r_i * |f_i,1 - f_i,2|/(f_i,1 + f_i,2))
  % where f_i,1 is the fraction of total spikes in pseudo-clust i at time 1
  for i = 1:nblocks-1
    tnpc = nperclust([i i+1],:);
    % Normalize
    tfpc = tnpc ./ repmat(sum(tnpc,2),1,size(tnpc,2));
    nzindx = find(sum(tfpc) > 0);
    tfpc = tfpc(:,nzindx);
    ell(i) = mean(clustsize(nzindx).*abs(diff(tfpc))./sum(tfpc));
  end
  