function options = cass_options
  % Note any of these will be overwritten by param/value pairs
  options.max_snip_memsize = 5e6;  % in # of doubles
  options.max_snips_to_cluster = 1e6;
  options.dirname = '';
  options.channel = 1;  % the channel to start with
  options.RfactorDefault = 1.5;
  options.corr_tmax = 3;  % 3 second max correlation time
  options.n_log_bins = 30;  % 30 bins in auto- & cross-correlation
  %options.cluster_func = @meanshift;
  options.cluster_func = @clust_em_climb;
  options.force_autosort_info = 0;  % 1 causes autosort_info to be chosen
  options.fig_positions = NaN;
  options.set_correlation_yaxes_without_fastest = 0;