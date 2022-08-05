function cluster_fits(fitFiles,sortfile,options)
% CLUSTER_FITS: cluster the amplitudes from spike waveform decomposition
% This function performs automated clustering (MSAMS) on spike
% amplitudes.
% Syntax:
%   cluster_fits(fitFiles,outputfile,options)
% where
%   fitfiles is a cell array of ".fitcomp" files
%   outputfile is a string containing the name to use for the results
%     file
%   options is a structure which may have the following fields:
%     mode (default 'minmax'): 'minmax' (performs clustering on each
%       spike's min/max voltage on each channel) or 'amplitudes' (performs
%       clustering on the projections onto the component waveforms).
%     t2V (default 0): if nonzero, the spike time is appended as a
%       coordinate in the clustering. The value gives the conversion
%       factor from time (in seconds) to component amplitude (i.e.,
%       voltage).
%     extra_coord_func: if present, this should be a cell array of
%       function handles, one for each fitFile.  For each file, this
%       function is called with syntax func(spiketimes); the result is
%       appended as extra coordinates for clustering.  This can be useful
%       if, for example, you want to append stimulus information (see
%       stim2coord)
%     subcluster (default false): if true, then the "big" clustering
%       (using the full dimensionality of the data set) is followed up by
%       a "subclustering" of each cluster, in which PCA is used to do
%       dimensionality reduction on each and then a further clustering is
%       performed.  This is done recursively.
%       When subcluster is true, the following fields are relevant:
%         min_per_cluster (default 10): clusters below this size won't be
%           further sub-clustered 
%         ndims_pca (default 3): the number of dimensions to use in the
%           PCA analysis
%         maxpoints_dimred (default 4000): the maximum number of points
%           used when calculating projection directions for dimensionality
%           reduction
%     runlda (default false): if true, computes a t-statistic (using LDA)
%       to compute a measure of the separation between every pair of
%       clusters. The result is stored in a matrix called ldamtrx.
%     max_components (default Inf): the maximum number of fitting
%       components to retain
%     n_landmarks (default 5000): the maximum number of landmarks to use
%       when clustering
%     normalize (default false): if true, normalize the spikes to the max
%       absolute value  added by HAA
%
% See also: STIM2COORD, FIT_COMPONENTS, AUTOMERGE_CLUSTERS_LDA.
  
% Copyright 2008 by Timothy E. Holy
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'mode','minmax','max_components',Inf,'t2V',0,'max_snips_to_cluster',Inf,'n_landmarks',20000,'subcluster',false,'min_per_cluster',10,'ndims_pca',3,'maxpoints_dimred',4000,'runlda',false, 'normalize',false);
  fu = fit_utilities;
  
  if ischar(fitFiles)
    fitFiles = {fitFiles};
  end
  n_files = length(fitFiles);
  
  % Load the spike data
  fprintf('Loading the data...');
  nspikes = nan(1,n_files);
  [amplitudes_all,times_all,residual_keys, headers] = load_amplitudes_times(fitFiles,options.mode); %#ok<NASGU>
  fprintf('done\n');
  for fileIndex = 1:n_files
    nspikes(fileIndex) = length(times_all{fileIndex});
  end
  
  % Concatenate into one big pile
  switch options.mode
   case 'amplitudes'
    amplitudes_all = cat(2,amplitudes_all{:});
    % Truncate dimensionality if requested
    max_components = min(options.max_components,size(amplitudes_all,1));
    amplitudes_all = amplitudes_all(1:max_components,:);
   case 'minmax'
    amplitudes_all = cat(3,amplitudes_all{:});
    amplitudes_all = fu.snip2vec_by_c(amplitudes_all);
  end  

  % HAA:Normalize the amplitudes to maximum absolute value
  if options.normalize
      norm_factor = max(abs(amplitudes_all),[],1);
      amplitudes_all = amplitudes_all ./ repmat(norm_factor,size(amplitudes_all,1),1);
  end
  
  % Optionally append extra coordinates
  if isfield(options,'extra_coord_func')
    ec = cell(1,n_files);
    for fileIndex = 1:n_files
      ec{fileIndex} = options.extra_coord_func{fileIndex}(times_all{fileIndex});
    end
    amplitudes_all = [amplitudes_all; cat(2,ec{:})];
  end

  times_all = scantime2abstime(times_all,headers);
  times_all = cat(2,times_all{:});
  
  % Optionally include the spike time as a coordinate (allows for
  % changing spike amplitudes/shapes)
  if options.t2V
    a = [amplitudes_all; times_all*options.t2V];
  else
    a = amplitudes_all;
  end

  % Optionally do the clustering with only a subset of the spikes
  % (Can be used to increase performance)
  skip = ceil(size(a,2)/options.max_snips_to_cluster);
  if (skip < 1)
     skip = 1;
  end
  keepIndex = 1:skip:size(a,2);
  asub = double(a(:,keepIndex));
  
  
  % Do the clustering
  lminfo = choose_landmarks(asub,options.n_landmarks,struct('frac_random',min(95,options.n_landmarks)/100));
  fprintf('Clustering...');
  clust = msams(asub,lminfo,struct('flow_only_landmarks',true,'consolidate',false,'reflow',true));
  landmarkClust = clust(lminfo.landmark_xIndex);
  fprintf('done\n');
  
  % Apply the results of the clustering to all spikes
  if (skip > 1)
    fprintf('Applying results to all spikes...');
    % For each spike, find the closest landmark
    [dist,lmIndex] = mindist(a,lminfo.landmarks);
    % Update the lminfo structure appropriately
    lminfo.landmarkAssignment = lmIndex;
    lminfo.landmarkList = agglabel(lmIndex);
    lminfo.dist_to_closest_landmark = dist;
    lminfo.landmark_xIndex = keepIndex(lminfo.landmark_xIndex);
    % Assign each spike to the closest landmark's cluster
    spikeClust = landmarkClust(lmIndex);
    fprintf('done\n');
  else
    spikeClust = landmarkClust(lminfo.landmarkAssignment);
  end
  
  if options.subcluster
    fprintf('Subclustering %d clusters...\n',max(landmarkClust));
    options.subcluster_level = 0;
    spikeClust = subcluster(a,spikeClust,options);
  end
  
  if options.runlda
    clabel = agglabel(spikeClust);
    fprintf('Running the 1-d LDA analysis on %d clusters\n',length(clabel));
    nClusters = length(clabel);
    xlda = cell(1,nClusters);
    for i = 1:nClusters
      xlda{i} = cluster_fits_subset(a(:,clabel{i}),options.maxpoints_dimred);
    end
    ldamtrx = zeros(nClusters,nClusters);
    for i = 1:nClusters
      fprintf('%d...',i);
      x1lda = xlda{i};
      if (size(x1lda,2) == 1)
        ldamtrx(i,[1:i-1 i+1:nClusters]) = inf;
        ldamtrx([1:i-1 i+1:nClusters],i) = inf;
        continue;
      end
      for j = 1:i-1
        x2lda = xlda{j};
        if (size(x2lda,2) == 1)
          ldamtrx(i,j) = inf;
          ldamtrx(j,i) = inf;
          continue
        end
        evec = lda([x1lda x2lda],...
          [ones(1,size(x1lda,2)) ones(1,size(x2lda,2))+1]);
        p1 = evec'*a(:,clabel{i});
        p2 = evec'*a(:,clabel{j});
        t = abs(mean(p1)-mean(p2))/sqrt(var(p1) + var(p2));
        ldamtrx(i,j) = t;
        ldamtrx(j,i) = t;
      end
    end
    fprintf('\n');
  end
  
   
  % Save the results
  fprintf('Saving the results...');
  cumspikes = [0 cumsum(nspikes)]; %#ok<NASGU>
  mode = options.mode; %#ok<NASGU>
  varsToSave = {'landmarkClust','lminfo','spikeClust', ...
    'cumspikes','options','fitFiles','mode','residual_keys'};
  if options.runlda
    varsToSave{end+1} = 'ldamtrx';
  end
  save('-mat',sortfile,varsToSave{:});
  fprintf('done\n');
  return
  
  % Generate a unique identifier so these are identified as a group
%   rand('state', sum(100*clock));  % generate _different_ keys on each matlab startup
%   cluster_key = round(1e8*rand(1,1));
%   for fileIndex = 1:n_files
%     outfile = replace_extension(fitFiles{fileIndex},'.fitclusters');
%     rng = cumspikes(fileIndex)+1:cumspikes(fileIndex+1);
%     cluster_label = clust(rng);
%     %spiketimes = times_all(rng);
%     residual_key = residual_keys(fileIndex);
%     save('-mat',outfile,'cluster_label','residual_key','cluster_key','spiketimes');
%   end

% This function will perform PCA/clustering on each cluster until it no
% longer splits any more
function clusternew = subcluster(x,clusterold,options)
  clusternew = clusterold;
  [clabel,nlabel] = agglabel(clusterold);
  nClusters = length(nlabel);
  for i = 1:nClusters
    indent = repmat(' ',1,options.subcluster_level);
    fprintf('%s%d (%d)\n',indent,i,nlabel(i));
    if (nlabel(i) < options.min_per_cluster)
      continue;
    end
    xc = x(:,clabel{i});
    tclust = cluster_pca_msams(xc,options);
    if (max(tclust) > 1)
      % recursively refine any cluster that gets split
      newops = options;
      newops.subcluster_level = newops.subcluster_level + 1;
      tclust = subcluster(xc,tclust,newops);
    end
    % Label with a unique tag
    clusternew(clabel{i}) = tclust + i/(nClusters+1);
  end
  % Convert to integer labels
  [ulabel,tmp,clusternew] = unique(clusternew);

function clust = cluster_pca_msams(x,options)
  % Do PCA, possibly on a subset of the data
  d = size(x,1);
  xpca = cluster_fits_subset(x,options.maxpoints_dimred);
  pd = pca(xpca');
  pd = pd(:,1:min(d,options.ndims_pca));
  proj = pd'*x;
  % Do the clustering
  lminfo = choose_landmarks(proj,options.n_landmarks,struct('frac_random',min(95,options.n_landmarks)/100,'progress',false));
  landmarkClust = msams(proj,lminfo);
  clust = landmarkClust(lminfo.landmarkAssignment);

function xskip = cluster_fits_subset(x,maxpoints)
  [d,N] = size(x);
  skip = ceil(N/maxpoints);
  if (skip < 1)
    skip = 1;
  end
  xskip = x(:,1:skip:end);
  