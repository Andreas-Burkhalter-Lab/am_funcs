function [components,options] = lda_sort_fixed(snip,options)
% LDA_SORT_FIXED: choose components from consistency across fixeds
%
% Fixed replicates are generated from an input set of spike
% waveforms.  LDA is performed on each replicate.  Two replicates are
% compared by computing the dot product between all principal
% components.  The absolute value of the dot product is averaged across
% all pairs of replicates.  The number of components is chosen by a
% threshold on the average dot product.
%
% Syntax:
%   [components,options] = lda_sort_fixed(snip,options)
% where
%   snip is a snipLength-by-nSnips matrix of spike snippets;
%   options is a structure array with the following possible fields:
%     lda_n_components (default 3): picks the number of components;
%     lda_n_clusters (default 2*n_components): the number of clusters
%       (picked by k-means) used to define the LDA directions
%
% On output,
%   components is a snipLength-by-nComponents matrix;
%   options is an output copy of the input structure options, with values
%     for all defaults supplied.
%
% See also: LDA, LDA_COMPARE.
  
  if ~isfield(options,'lda_n_components')
    options.lda_n_components = 3;
  end
  if ~isfield(options,'lda_n_clusters')
    options.lda_n_clusters = 2*options.lda_n_components;
  end
  
  indx = kmeans_hard(snip,options.lda_n_clusters);
  components = lda(snip,indx);

  % Keep only the number requested by the user
  n_to_keep = min(options.lda_n_components,size(components,2));
  components = components(:,1:n_to_keep);
