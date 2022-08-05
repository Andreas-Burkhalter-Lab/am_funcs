function [components,options] = pca_sort_fixed(snip,options)
% PCA_SORT_FIXED: choose components from consistency across fixeds
%
% Fixed replicates are generated from an input set of spike
% waveforms.  PCA is performed on each replicate.  Two replicates are
% compared by computing the dot product between all principal
% components.  The absolute value of the dot product is averaged across
% all pairs of replicates.  The number of components is chosen by a
% threshold on the average dot product.
%
% Syntax:
%   [components,options] = pca_sort_fixed(snip,options)
% where
%   snip is a snipLength-by-nSnips matrix of spike snippets;
%   options is a structure array with the following possible fields:
%     pca_nu (default 0.5): controls the robustness of the PCA, see PCA;
%     pca_n_components (default 3): picks the number of components;
%
% On output,
%   components is a snipLength-by-nComponents matrix;
%   options is an output copy of the input structure options, with values
%     for all defaults supplied.
%
% See also: PCA, PCA_COMPARE.
  
  if ~isfield(options,'pca_nu')
    options.pca_nu = 0.5;  % Use robust form of PCA
  end
  if ~isfield(options,'pca_n_components')
    options.pca_n_components = 2;
  end
  
  % Power-law scale the data for more robust PCA
  Y = plsdata(snip',options.pca_nu);
  components = pca(Y,struct('nocenter',1));
  n_to_keep = min(options.pca_n_components,size(components,2));
  components = components(:,1:n_to_keep);
