function [components,options] = pca_sort_bootstrap_AM(snip,options,progressbar_on)
% PCA_SORT_BOOTSTRAP: choose components from consistency across bootstraps      
    %%% AM edited 8/10/15 on vivid
    %%% -added option to suppress progressbar
% Bootstrap replicates are generated from an input set of spike
% waveforms.  PCA is performed on each replicate.  Two replicates are
% compared by computing the dot product between all principal
% components.  The absolute value of the dot product is averaged across
% all pairs of replicates.  The number of components is chosen by a
% threshold on the average dot product.
%
% Syntax:
%   [components,options] = pca_sort_bootstrap(snip,options)
% where
%   snip is a snipLength-by-nSnips matrix of spike snippets;
%   options is a structure array with the following possible fields:
%     pca_nu (default 0.5): controls the robustness of the PCA, see PCA;
%     pca_n_bootstraps (default 10): the number of replicate resamplings
%       of the data;
%     pca_dotprod_thresh (default 0.9): the threshold on the average
%       |dotproduct| of components in pairwise comparisons between
%       replciates;
%     pca_min_components (default 3): the minimum number of components to
%       choose;
%     pca_bootstrap_plot (default false): if true, causes the dot product
%       to be plotted.
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
  if ~isfield(options,'pca_n_bootstraps')
    options.pca_n_bootstraps = 10;
  end
  if ~isfield(options,'pca_dotprod_thresh')
    options.pca_dotprod_thresh = 0.9;
  end
  if ~isfield(options,'pca_min_components')
    options.pca_min_components = 3;
  end
  toplot = 0;
  if isfield(options,'pca_bootstrap_plot')
    toplot = options.pca_bootstrap_plot;
  end
  progress_options = struct('max',options.pca_n_bootstraps+1,...
    'what','Finding components by PCA bootstrap');
  progress_options.progress = 0;
  
        %% AM 8/10/15 added conditional
    if exist('progressbar_on','var') && progressbar_on  
        progress_options = progress(progress_options);
    end
  
  % Power-law scale the data for more robust PCA
  Y = plsdata(snip',options.pca_nu);
  % Now do the replicates of PCA
  N = size(Y,1);
  for i = 1:options.pca_n_bootstraps
    bsIndex = round(N*rand(1,N)+0.5);
    pcaboot(:,:,i) = pca(Y(bsIndex,:),struct('nocenter',1));

        %% AM 8/10/15 added conditional
    if exist('progressbar_on','var') && progressbar_on
        progress_options.progress = i+1;
        progress_options = progress(progress_options);
    end
    %%
    
  end
  % Compare the components from the different replicates
  dp = pca_compare(pcaboot);
  % Find the first component with less than threshold consistency
  n_to_keep = find(dp < options.pca_dotprod_thresh,1,'first')-1;
  if (isempty(n_to_keep) || n_to_keep < options.pca_min_components)
    n_to_keep = options.pca_min_components;
  end
  % Check to make sure there aren't any later components with larger
  % consistency
  if (find(dp > options.pca_dotprod_thresh,1,'last') > n_to_keep)
    warning('Inconsistency in dot product thresholds')
    %toplot = 1;
  end
  n_to_keep = max(n_to_keep,options.pca_min_components);
  % Do one last PCA on the whole data set & keep only selected components
  components = pca(Y);
  components = components(:,1:n_to_keep);
  % Plot, if indicated
  if toplot
    figure
    semilogx(dp,'b.')
    line(1:n_to_keep,dp(1:n_to_keep),...
         'Color','r','LineStyle','none','Marker','.');
    line([1 length(dp)],[1 1]*options.pca_dotprod_thresh,...
         'Color','k','LineStyle','--');
  end
  
          %% AM 8/10/15 added conditional
    if exist('progressbar_on','var') && progressbar_on
        close(progress_options.handle)
    end
    %%