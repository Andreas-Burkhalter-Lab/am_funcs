function cluster_fits(fitFiles,sortfile,options)
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'t2V',0,'max_snips_to_cluster',Inf,'n_landmarks',5000);
  
  if ischar(fitFiles)
    fitFiles = {fitFiles};
  end
  n_files = length(fitFiles);
  
  % Load the spike data
  nspikes = nan(1,n_files);
  [amplitudes_all,times_all,residual_keys] = load_amplitudes_times(fitFiles);
  for fileIndex = 1:n_files
    nspikes(fileIndex) = length(times_all{fileIndex});
  end
  % Concatenate into one big pile
  amplitudes_all = cat(2,amplitudes_all{:});
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
  asub = a(:,1:skip:end);
  
  % Do the clustering
  lminfo = choose_landmarks(asub,options.n_landmarks);
  landmarkClust = msams(asub,lminfo);
  
  % Apply the results of the clustering to all spikes
  if (skip > 1)
    % For each spike, find the closest landmark
    [dist,lmIndex] = mindist(a,lminfo.landmarks);
    % Assign each spike to the closest landmark's cluster
    spikeClust = landmarkClust(lmIndex);
  else
    spikeClust = landmarkClust(lminfo.landmarkAssignment);
  end
   
  % Save the results
  cumspikes = [0 cumsum(nspikes)];
  save('-mat',sortfile,'landmarkClust','lminfo','spikeClust','cumspikes','options','fitFiles','residual_keys');
  return
  
  % Generate a unique identifier so these are identified as a group
  rand('state', sum(100*clock));  % generate _different_ keys on each matlab startup
  cluster_key = round(1e8*rand(1,1));
  for fileIndex = 1:n_files
    outfile = replace_extension(fitFiles{fileIndex},'.fitclusters');
    rng = cumspikes(fileIndex)+1:cumspikes(fileIndex+1);
    cluster_label = clust(rng);
    %spiketimes = times_all(rng);
    residual_key = residual_keys(fileIndex);
    save('-mat',outfile,'cluster_label','residual_key','cluster_key','spiketimes');
  end
  