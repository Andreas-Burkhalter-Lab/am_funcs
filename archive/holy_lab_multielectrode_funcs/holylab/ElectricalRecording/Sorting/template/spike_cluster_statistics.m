function stats = spike_cluster_statistics(filename,options)
% SPIKE_CLUSTER_STATISTICS: compute metrics of cluster quality for single units
% This function computes a measure of autocorrelation cleanliness (Omega)
% and separability (dprime).
% Syntax:
%    spike_cluster_statistics(filename)
% When filename is a cluster_nav file, this opens the spike data, computes
% the metrics, and saves the results to the same file. This prepares the
% file for automated thresholding of spike quality
%    stats = spike_cluster_statistics(filename)
% instead returns the data to the command line (and does not alter the
% file). "stats" is a structure with the following fields:
%   Autocorrelation-related fields:
%   omega: a scalar indicating the "purity" of a given autocorrelation
%     function (ranges from 0 to 1, with 1 indicating no refractory
%     violations)
%   n_tac: a 1-by-n_cells cell array, each element containing the
%     autocorrelation histogram for the given cluster
%   Spike shape-related fields:
%   dprime: an n_cells-by-n_cells matrix indicating the 1-dimensional
%     discriminability between spike waveforms
%   nbrIDs: an n-cells-by-2 matrix containing the names of the two nearest
%     neighbor clusters as judged by dprime
%   proj: a cell array of length n_cells, each element containing the
%     2-by-N matrix of projections of waveform shape onto the directions
%     determined by LDA to best discriminate among the given cluster and
%     its 2 nearest neighbors (as indexed by nbrIDs)
%   clustbreaks: an n_cells-by-3 matrix indicating the index in proj of the
%     last spike of the cluster. So
%     proj{i}(:,1:clustbreaks(i,1)) contains the projections of the ith
%     cell,
%     proj{i}(:,clustbreaks(i,1)+1:clustbreaks(i,2)) contains the
%     projections of its closest neighbor, and 
%     proj{i}(:,clustbreaks(i,2)+1:clustbreaks(i,3)) contains the
%     projections of its next-nearest neighbor.
%   meanwaveforms is a n_channels-by-n_times-by-n_cells matrix containing
%     the mean waveform of each cluster.
%   options: a copy of the settings used in calling this function.
%
% A second input argument, options, is a structure which may have the
% following fields:
%   refractory period (default 0.02): in seconds
%   plateau (default 1): in seconds, a time longer than the burst-induced
%     autocorrelation peak.
%   n_bins (default 500): the number of time bins to use in returning n_tac
%   do_projections (default true): if false, the nbrIDs, proj, and
%     clustbreaks fields are not set in "stats"
%   do_meanwaveforms (default true): if false, the meanwaveforms field is
%     not calculated

% Copyright 2008 by Timothy E. Holy

  %% Argument parsing
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'refractory_period',0.02,'plateau',1,'n_bins',500,'do_projections',true,'do_meanwaveforms',true);
  saving_to_file = false;
  if (nargout == 0)
    saving_to_file = true;
  end
  
  %% Load the data
  clustdata = load('-mat',filename);
  n_cells = max(clustdata.spikeClust);
  % Put stuff here that checks to see if we need to back up to get the
  % spike times
  if isfield(clustdata,'preClusterNavFile')
    spikeClust = clustdata.spikeClust;
    clusterIDs = clustdata.clusterIDs(2:end);  % skip the 'noise' cluster
    clustdata = load('-mat',clustdata.preClusterNavFile);
    clustdata.spikeClust = spikeClust;
  else
    clusterIDs = cell(1,n_cells);
    for cellIndex = 1:n_cells
      clusterIDs{cellIndex} = ['c' num2str(cellIndex)];
    end
  end
  
  n_files = length(clustdata.fitFiles);
  
  % Load the raw data
  for fileIndex = 1:n_files
    ffdata(fileIndex) = load('-mat',clustdata.fitFiles{fileIndex}); %#ok<AGROW>
    h(fileIndex) = readheader(ffdata(fileIndex).fileToFit); %#ok<AGROW>
  end
  
  %% Determine which spikes belong to which cells in each file
  spikeClust = cell(1,n_files);
  groups = cell(n_files,n_cells);
  n_per_group = zeros(n_files,n_cells);
  for fileIndex = 1:n_files
    spikeClust{fileIndex} = clustdata.spikeClust(clustdata.cumspikes(fileIndex)+1:clustdata.cumspikes(fileIndex+1));
    [gtmp,npgtmp] = agglabel(spikeClust{fileIndex});
    ngtmp = length(npgtmp);
    groups(fileIndex,1:ngtmp) = gtmp;
    n_per_group(fileIndex,1:ngtmp) = npgtmp;
  end
  
  %% Do the temporal autocorrelation
  % Initialize storage for autocorrelation analysis of cells
  n_refract = zeros(n_files,n_cells);  % # refractory violations
  n_plateau = zeros(n_files,n_cells);  % # in plateau region
  n_tac = zeros(n_cells,options.n_bins);
  
  for cellIndex = 1:n_cells
    for fileIndex = 1:n_files
      t = ffdata(fileIndex).spiketimes(groups{fileIndex,cellIndex}) / h(fileIndex).scanrate;
      n_refract(fileIndex,cellIndex) = autocorrspike(t,options.refractory_period,1);
      n_plateau_tmp = autocorrspike(t,2*options.plateau,2);
      n_plateau(fileIndex,cellIndex) = n_plateau_tmp(2);
      n_tac(cellIndex,:) = n_tac(cellIndex,:) + autocorrspike(t,options.plateau,options.n_bins);
    end
  end
  % Calculate the "fraction clean" statistic on these quantities
  tratio = options.plateau/options.refractory_period;
  denom = sum(n_plateau,1);
  num = sum(n_refract,1);
  num(denom == 0) = 1; denom(denom == 0) = 1;
  omega = sqrt(1 - num ./ denom * tratio);
  omega(omega ~= real(omega)) = 0;
    
  %% Do the 1-d LDA
  % Calculate the stride for the 1-d LDA
  d = size(clustdata.lminfo.landmarks,1);     % # of dimensions used for clustering
  max_for_lda = 2000+d;                % don't use any more waveforms than this for 1-d LDA
  skip = floor(sum(n_per_group,1) / max_for_lda);
  skip(skip < 1) = 1;
  
  % Collect the waveforms
  wf = cell(1,n_cells);          % "waveform" for LDA
  wftmp = cell(1,n_files);
  nwf = zeros(1,n_cells);
  for cellIndex = 1:n_cells
    for fileIndex = 1:n_files
      wftmp{fileIndex} = fetch_wfdata(ffdata(fileIndex),groups{fileIndex,cellIndex},skip(cellIndex),clustdata.mode,clustdata.options);
    end
    wf{cellIndex} = cat(2,wftmp{:});
    nwf(cellIndex) = size(wf{cellIndex},2);
  end

  % Calculate the 1d LDAs
  dprime_mtrx = zeros(n_cells,n_cells);
  for cellIndex = 1:n_cells
    for cellIndex2 = cellIndex+1:n_cells
      group_label = [ones(1,nwf(cellIndex)), 1+ones(1,nwf(cellIndex2))];
      spike_data = [wf{cellIndex} wf{cellIndex2}];
      eigvec = lda(spike_data,group_label);
      proj_tmp = eigvec'*spike_data;
      % Quantify the overlap in terms of a d' value
      proj1 = proj_tmp(group_label == 1);
      proj2 = proj_tmp(group_label == 2);
      dprime = abs(mean(proj1) - mean(proj2)) / sqrt(var(proj1) + var(proj2));
      dprime_mtrx(cellIndex,cellIndex2) = dprime;
      dprime_mtrx(cellIndex2,cellIndex) = dprime;
    end
  end
  
  %% Compute the 2-d projections
  if options.do_projections
    proj = cell(1,n_cells);
    clustbreaks = zeros(n_cells,3);
    % Re-collect all of the waveform data
    for cellIndex = 1:n_cells
      for fileIndex = 1:n_files
        wftmp{fileIndex} = fetch_wfdata(ffdata(fileIndex),groups{fileIndex,cellIndex},1,clustdata.mode,clustdata.options);
      end
      wf{cellIndex} = cat(2,wftmp{:});
      nwf(cellIndex) = size(wf{cellIndex},2);
    end
    % For each cell, pick the two closest neighbors and do LDA
    nbrIDs = cell(n_cells,2);
    nbrIndex_all = zeros(n_cells,2);
    for cellIndex = 1:n_cells
      [dpsort,sortIndex] = sort(dprime_mtrx(cellIndex,:));
      nbrIndex = sortIndex(2:3);
      nbrIndex_all(cellIndex,:) = nbrIndex;
      nbrIDs(cellIndex,:) = clusterIDs(nbrIndex);
      group_label = [ones(1,nwf(cellIndex)), 1+ones(1,nwf(nbrIndex(1))), 2+ones(1,nwf(nbrIndex(2)))];
      spike_data = [wf{cellIndex} wf{nbrIndex(1)} wf{nbrIndex(2)}];
      eigvec = lda(spike_data,group_label);
      proj{cellIndex} = eigvec'*spike_data;
      clustbreaks(cellIndex,:) = cumsum(nwf([cellIndex nbrIndex]));
    end
  end

  %% Calculate the mean waveforms
  if options.do_meanwaveforms
    % Collect the mean amplitudes
    mean_amp = zeros(size(ffdata(1).amplitudes,1),n_cells);
    for cellIndex = 1:n_cells
      wftmp = cell(1,n_cells);
      for fileIndex = 1:n_files
        wftmp{fileIndex} = fetch_wfdata(ffdata(fileIndex),groups{fileIndex,cellIndex},1,'amplitudes',struct('normalize',false));
      end
      mean_amp(:,cellIndex) = mean(cat(2,wftmp{:}),2);
    end
    % Use the components to convert amplitudes to waveform shapes
    cF = load('-mat',ffdata(1).componentFile);
    wfc = cF.templates;
    sz = size(wfc);
    fu = fit_utilities;
    wfc = fu.snip2vec_by_c(wfc);
    meanwaveforms = wfc*mean_amp;
    meanwaveforms = fu.vec2snip_by_c(meanwaveforms,sz(1));
  end
  
  %% Prepare the output
  stats = struct('omega',omega,...
    'n_tac',n_tac,...
    'dprime',dprime_mtrx);
  if options.do_projections
    stats.nbrIndex = nbrIndex_all;
    stats.nbrIDs = nbrIDs;
    stats.proj = proj;
    stats.clustbreaks = clustbreaks;
  end
  if options.do_meanwaveforms
    stats.meanwaveforms = meanwaveforms;
  end
  stats.options = options;
  
  %% Save the data to disk
  if saving_to_file
    clustdata.stats = stats;
    save(filename,'-mat',clustdata);
  end
end


function wftmp = fetch_wfdata(ffdata,groups,skip,mode,options)
  if strcmp(mode,'minmax')
    wftmp = ffdata.minmax(:,:,groups(1:skip:end));
    sz = size(wftmp);
    wftmp = reshape(wftmp,[sz(1)*sz(2) sz(3)]);
  else
    wftmp = ffdata.amplitudes(:,groups(1:skip:end));
  end
  if isfield(options,'normalize') && options.normalize
    wftmp = wftmp ./ repmat(max(abs(wftmp),[],1),1,size(wftmp,2));
  end
end