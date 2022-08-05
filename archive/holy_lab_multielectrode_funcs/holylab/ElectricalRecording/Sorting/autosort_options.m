function options = autosort_options(options)
% AUTOSORT_OPTIONS: options structure defaults for autosorting
% Syntax:
%   options = autosort_options(options)
% You supply an input set of options, and this function fills in the
% defaults for any unsupplied options.
%
% Here are the recognized fields and their default values:
%
% Fields dealing with resource consumption:
%   max_snip_memsize (default 5e6): the number of values to allocate to
%     snippets.  This is the product of the sniplength and the number of
%     snippets.  Note this is the number of doubles (or whatever format
%     they're loaded in), not in bytes: 3e6 corresponds to 24Mb in
%     doubles.  Note that this refers only to snippet waveforms, other
%     variables will consume additional memory.
%   min_snips_to_cluster (default 100): if there are fewer snippets than
%     this, don't even bother with this channel.
%   max_snips_to_cluster (default 20000): the name says it all...almost,
%     that is.  AUTOSNIP will not use more memory for holding snippets
%     than specified by max_snip_memsize.  This parameter allows you to
%     restrict the total number of snippets directly, if desired. Set to
%     "Inf" if you don't want to restrict the number by anything other than
%     memory limitations.
%   Also note that "n_landmarks" and "n_snips_projection" (described
%   below) have impact on resource consumption.
%
% Fields dealing with dimensionality reduction:
%   use_projection (default true): if true, reduce dimensionality by
%     projecting the spike waveforms onto a set of filters.
%   n_snips_projection (default 2000): the maximum number of waveforms to
%     use in calculating the filters.
%   projection_func (default @pca_sort_bootstrap): a function handle
%     specifying the function used to calculate the projection directions
%     (filters) from the spike waveforms.  This function must have the
%     syntax:
%         filters = projection_func(waveforms)
%     where spike waveforms are the columns of 'waveforms', and each column of
%     'filters' is a projection direction.
%     Good choices are @pca_sort_fixed and @pca_sort_bootstrap.  There
%     are additional possible options depending on which projection
%     function you choose; see the relevant help for each choice.
%
% Fields dealing with the tradeoff between time and voltage in clustering:
%   t2V (by default, not set): the constant of proportionality converting
%     differences in spike time to differences in voltage.  Spike time
%     can be used as an additional coordinate in sorting, to give a very
%     fluid mechanism for tracking changes in spike shape/amplitude over
%     time.  You can supply this value directly, or you can calculate it
%     from the data (see next field). You can also set it interactive;
%     see the option below for 'interactive_landmarking'.
%   t2V_func (default @autosort_calc_t2V): a function handle to use in
%     calculating t2V from the data set.  See autosort_calc_t2V for the
%     options that this function takes; or write your own.  If you write
%     your own, it must have the syntax described in the help for
%     AUTOSORT_CALC_T2V.
%   Note also that "interactive_landmarking" (described below) has an
%   influence on t2V as well.
%
% Fields dealing with landmark choice:
%   landmark_func (default @sort_choose_landmarks): the function to call
%     for assigning landmarks.  Two good choices are
%     @sort_choose_landmarks and @kmeans_hard.  If you're clustering with
%     ADAPTIVE_MEANSHIFT, you probably want the former.  See
%     SORT_CHOOSE_LANDMARKS for the required interface to the
%     landmark_func.
%   n_landmarks (default 1000): the number of landmarks to use.  This is
%     a sensible choice when the clustering function is ADAPTIVE_MEANSHIFT.
%     For other choices of clustering function, a warning will be
%     generated if n_landmarks > 200 due to probable long computation times.
%   interactive_landmarking (default false): if true, a GUI will pop up asking
%     you to approve the landmarks and the value of t2V; you can edit t2V
%     if desired.  Note an older name, "interactive_kmeans," is also
%     supported.
%   n_replicates (default 1): landmarks will be assigned, e.g., by K-means,
%     which introduces an element of randomness.  You can see how
%     consistent the clustering is across sets of landmarks by doing the
%     sorting multiple times.  The sorting GUI, CASS, allows you to
%     select among replicates.  Note: with larger numbers of landmarks,
%     the need for replicates is reduced.
%
% Fields dealing with clustering algorithm:
%   cluster_func (default @adaptive_meanshift): function to use in aggregating
%     landmarks for doing the actual clustering.  This function must have
%     the syntax:
%        clustIndex = cluster_func(snip,landmark,Rfactor,options)
%     Two good choices are @adaptive_meanshift and @clust_em_climb.
%   Rfactor (default 2): a scalar or vector, giving a multiplier(s) used
%     in the amount of smoothing performed in flowing points up the
%     density gradient.  If you supply a vector, the clustering will be
%     performed for all the choices, which you can browse later.
%     This parameter has no influence when adaptive_meanshift clustering
%     is used. (??: interpret it as the "factor" parameter?)
%   ploteach, plotclust, plotpause, plotclustpause: these all affect the
%     amount of graphical feedback given during clustering.  See
%     MEANSHIFT and RMEANS for a description.
%
% See also: AUTOSNIP, AUTOSORT_CALC_T2V, SORT_CHOOSE_LANDMARKS,
% KMEANS_HARD, PCA_SORT_FIXED, PCA_SORT_BOOTSTRAP,
% CLUST_EM_CLIMB, ADAPTIVE_MEANSHIFT, MEANSHIFT.

% Copyright 2005-2006 by Timothy E. Holy
  
  options = autosort_meanshift_options(options);
  options = autosort_landmark_options(options);
  if ~isfield(options,'max_snip_memsize')
    options.max_snip_memsize = 5e6;  % in # of doubles
  end
  if ~isfield(options,'min_snips_to_cluster')
    options.min_snips_to_cluster = 100;
  end
  if ~isfield(options,'max_snips_to_cluster')
    options.max_snips_to_cluster = 20000;
  end
  if ~isfield(options,'n_replicates')
    options.n_replicates = 1;
  end
  if ~isfield(options,'t2V_func')
    options.t2V_func = @autosort_calc_t2V;
  end
  if ~isfield(options,'use_projection')
    options.use_projection = 1;
  end
  options = autosort_projection_options(options);
  
  
function options = autosort_projection_options(options)
  if ~isfield(options,'projection_func')
    options.projection_func = @pca_sort_bootstrap;
  end
  if ~isfield(options,'n_snips_projection')
    options.n_snips_projection = 2000;
  end
  if isequal(options.projection_func,@pca_sort_bootstrap)
    if ~isfield(options,'pca_bootstrap_plot')
      options.pca_bootstrap_plot = 0;
    end
  end


function options = autosort_meanshift_options(options)
  if ~isfield(options,'cluster_func')
    %options.cluster_func = @clust_em_climb;
    options.cluster_func = @adaptive_meanshift;
  end
  if ~isfield(options,'Rfactor')
    % Unused?
    options.Rfactor = 2;
  else
    options.Rfactor = sort(options.Rfactor);
  end
  if ~isfield(options,'ploteach')
    options.ploteach = 0;
  end
  if ~isfield(options,'plotclust')
    options.plotclust = 0;
  end
  if ~isfield(options,'plotpause')
    options.plotpause = 0;
  end
  if ~isfield(options,'plotclustpause')
    options.plotclustpause = 0;
  end

function options = autosort_landmark_options(options)
  if ~isfield(options,'n_landmarks')
    if isequal(options.cluster_func,@adaptive_meanshift)
      options.n_landmarks = 1000;
    else
      options.n_landmarks = 40;
    end
  end
  if ~isfield(options,'landmark_func')
    options.landmark_func = @sort_choose_landmarks;
    % another choice is kmeans_hard
  end
  % Copy the old name "interactive_kmeans" to "interactive_landmarking"
  if isfield(options,'interactive_kmeans')
    options.interactive_landmarking = options.interactive_kmeans;
  end
  if ~isfield(options,'interactive_landmarking')
    options.interactive_landmarking = false;
  end