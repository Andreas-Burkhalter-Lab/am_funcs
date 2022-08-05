function autosort_template(savedir,fitcompFiles,options)
% AUTOSORT_TEMPLATE: automatic spike sorting of multielectrode recordings
%
% To sort data and save to disk, use the following syntax:
%   autosort_template(savedir,fitcompFiles,options)
%   autosort_template(savedir)   (to resume after abort)
% where
%   savedir is the name of the directory to be created to contain the
%     results;
%   fitcompFiles is a cell array of .fitcomp files (see FIT_COMPONENTS)
%   options is a structure which may have the fields described in
%     AUTOSORT_TEMPLATE_OPTIONS.
%
% See also: AUTOSORT_TEMPLATE_OPTIONS, AUTOSORT_CALC_T2V.
  
% Copyright 2007 by Timothy E. Holy
% Based on AUTOSORT.
  
  if (nargin < 3)
    options = struct;
  end
  options = autosort_template_options(options);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Do the "big clustering" which decides how to parcel out all events
  % into "fake channels."  This part of the clustering is (by default)
  % done on _normalized events_, emphasizing the distribution of voltage
  % across electrodes (and of course, temporal shape info) rather than
  % the total amplitude. (Total amplitude will be used in subclustering
  % each fake channel---the main reason for this design is to make sure that
  % cases where spike amplitudes decrease during bursts at least end up
  % having all the reasonably-sized spikes appear in the same fake
  % channel, so that they are presented as a unit and users
  % can easily see what is happening.  It's much harder to see what's
  % happening if they end up in different fake channels.)
  %
  if exist(savedir,'dir') && (nargin == 1)
    load(fullfile(savedir,'overview'));
    n_files = length(sorthead);
    n_fake_channels = length(unique(landmarkClust));
    componentFileInfo = load('-mat',componentFile);
  else
    fprintf('Performing the "big sort" to create fake channels:\n');
    % Before we do any real work, make the saving directory. That way
    % permissions errors get noticed without delay.
    [status,msg] = mkdir(savedir);
    if ~status
      error(['Can''t make directory ' savedir ': ' msg]);
    end
    if ~isempty(msg)
      warning(msg);
    end
    
    if ischar(fitcompFiles)
      fitcompFiles = {fitcompFiles};
    end
    n_files = length(fitcompFiles);
    % Find the time at which recording started in each file
    for fileIndex = 1:n_files
      tmp = load('-mat',fitcompFiles{fileIndex},'fileToFit');
      header(fileIndex) = readheader(tmp.fileToFit);
      tstart(fileIndex) = datenum([header(fileIndex).date ' ' header(fileIndex).time])*(3600*24);
    end
    % Sort the files in temporal order
    [tstart,fileOrder] = sort(tstart);
    fitcompFiles = fitcompFiles(fileOrder);
    header = header(fileOrder);
    % Figure out how many snippets are in each fitcomp file; also load
    % information about related files
    n_spikes_per_file = zeros(1,n_files);
    componentFile = cell(1,n_files);
    residualFile = cell(1,n_files);
    for fileIndex = 1:n_files
      tmp = load('-mat',fitcompFiles{fileIndex},'spiketimes','componentFile','residualFile');
      n_spikes_per_file(fileIndex) = length(tmp.spiketimes);
      componentFile{fileIndex} = tmp.componentFile;
      if isfield(tmp,'residualFile')
        residualFile{fileIndex} = tmp.residualFile;
      end
    end
    componentFile = unique(componentFile);
    if (length(componentFile) > 1)
      error('All files must use the same component file');
    end
    componentFile = componentFile{1};
    componentFileInfo = load('-mat',componentFile);
    % Determine the fraction of spikes to load
    n_spikes = sum(n_spikes_per_file);
    subset_fraction = options.max_snips_to_cluster / n_spikes;
    if (subset_fraction > 1)
      subset_fraction = 1;
    end
    stride = round(1/subset_fraction);
    % Load the spike projections
    amplitudes = cell(1,n_files);
    for fileIndex = 1:n_files
      tmp = load('-mat',fitcompFiles{fileIndex},'amplitudes');
      amplitudes{fileIndex} = tmp.amplitudes(:,1:stride:end);
    end
    amplitudes = cat(2,amplitudes{:});
    tpl_proj = componentFileInfo.projections;
    if options.normalize
      ampl_norm = sqrt(sum(amplitudes.^2,1));
      amplitudes = amplitudes ./ repmat(ampl_norm,size(amplitudes,1),1);
      tpl_norm = sqrt(sum(tpl_proj.^2,1));
      tpl_proj = tpl_proj ./ repmat(tpl_norm,size(tpl_proj,1),1);
    end
    % Choose landmarks
    fprintf('  Choosing landmarks...\n');
    lminfo = choose_landmarks(amplitudes,options.n_landmarks, ...
			      struct('seed_landmarks',tpl_proj));
    % Do the "big" autoclustering
    fprintf('  Clustering (landmark progress, out of %d): ',options.n_landmarks);
    landmarkClust = msams(amplitudes,lminfo,struct('variable_metric',false,'show_progress',true));
    % Now we know how the landmarks group together, and the number of
    % distinct clusters will be the "channels" that the viewer sees in
    % CASS.
    n_fake_channels = length(unique(landmarkClust));
    % Create a fake sorthead structure array
    for fileIndex = 1:n_files
      [fcpathstr,fcfilebase,fcext] = fileparts(fitcompFiles{fileIndex});
      if isempty(fcpathstr) || fcpathstr(1) ~= filesep
        fcpathstr = [pwd filesep fcpathstr];
      end
      sorthead(fileIndex) = struct('type','fit_component',...
        'nscans',header(fileIndex).nscans,...
        'numch',n_fake_channels,...
        'channels',1:n_fake_channels,...
        'real_channels',componentFileInfo.channels,...
        'real_thresh',componentFileInfo.thresh,...
        'scanrate',header(fileIndex).scanrate,...
        'sniprange',componentFileInfo.sniprange,...
        'date',header(fileIndex).date,...
        'time',header(fileIndex).time,...
        'thresh',median(componentFileInfo.thresh)*ones(1,n_fake_channels),...
        'scalemult',1,...
        'scaleoff',0,...
        'polarity',0,...
        'voltageMin',header(fileIndex).voltageMin,...
        'voltageMax',header(fileIndex).voltageMax,...
        'dirname',savedir,...
        'fh_fitcomp',filehandle('filename',[fcfilebase fcext], ...
				'abspathstr',fcpathstr));
      % The main issue is that the thresh is handled in a cheesy way, but
      % it's not entirely clear what the _right_ way is...
      if ~isempty(residualFile{fileIndex})
        [respathstr,resfilebase,resext] = fileparts(residualFile{fileIndex});
        if isempty(respathstr) || respathstr(1) ~= filesep
          respathstr = [pwd filesep respathstr];
        end
        sorthead(fileIndex).fh_residual = filehandle('filename',[resfilebase resext],'abspathstr',respathstr);
      end
    end
    % Save the overview file
    save_overview_file(savedir);
  end  % if exist(savedir)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Now populate the individual fake channels with events. Another key
  % aspect to the design is that a single event can end up in more than
  % one fake channel---the idea being that when working on a given fake
  % channel, you can see "all" events that might possibly be related.
  % For each event, we'll calculate the closest landmarks, and then put a
  % copy of the event in each of the corresponding fake channels. But
  % we'll also retain a notion of how the assignment ranks among all of
  % its possible assignments.
  %
  fullpath = [savedir filesep 'chan' num2str(n_fake_channels)];
  lastfile = fullfile(fullpath,['time_amp' num2str(n_files) '.mat']);
  fileToFit = '';
  optionsForFit = struct;
  spiketimes = [];
  if ~exist(lastfile,'file')
    fprintf('\nAssigning events to fake channels\n');
    % Make subdirectories
    for channelIndex = 1:n_fake_channels
      mkdir([savedir filesep 'chan' ...
	     num2str(channelIndex)])
    end
    for fileIndex = 1:n_files
      fprintf('Parceling %s.',sorthead(fileIndex).fh_fitcomp.filename);
      load('-mat',fullfile(sorthead(fileIndex).fh_fitcomp));
      fprintf('.');
      % Compute distances between spike component amplitudes and the
      % landmarks
      amplitudes_sort = amplitudes;  % we'll need to save amps for reconstruction
      if options.normalize
        ampl_norm = sqrt(sum(amplitudes.^2,1));
        amplitudes_sort = amplitudes ./ repmat(ampl_norm,size(amplitudes,1),1);
      end
      sd = sqrdist(lminfo.landmarks,amplitudes_sort);
      fprintf('.');
      % Sort the distance
      [ssd,sortIndex] = sort(sd);
      fprintf('.');
      % Find the fake channels that the closest n landmarks go to
      spikeClust = landmarkClust(sortIndex(1:options.n_closest_landmarks,:));
      fprintf('.');
      % For each spike, dump it in all relevant "fake channel" bins
      fake_channel_spike_index = cell(1,n_fake_channels);
      fake_channel_spike_rank = cell(1,n_fake_channels);
      for candidateIndex = 1:options.n_closest_landmarks
        tmp = agglabel(spikeClust(candidateIndex,:));
        for chanIndex = 1:n_fake_channels
          [fake_channel_spike_index{chanIndex},indx] = unique([...
            fake_channel_spike_index{chanIndex} tmp{chanIndex}]);
          tmprank = [fake_channel_spike_rank{chanIndex},...
            repmat(candidateIndex,1,length(tmp{chanIndex}))];
          fake_channel_spike_rank{chanIndex} = tmprank(indx);
          % Highest rank (i.e., lowest rank #) will be retained this way
          % due to ordering provided by the call to "unique"
        end
      end
      fprintf('.');
      % Now store times & amplitudes in each fake channel (it's dangerous
      % to assume we can store everything from all files in RAM, yet
      % re-loading for each file takes time. So it would seem to be worth
      % splitting out into separate files for each fake channel
      for chanIndex = 1:n_fake_channels
        fn = fullfile([savedir filesep 'chan' num2str(chanIndex)],['time_amp' num2str(fileIndex)]);
        fc_spiketimes = spiketimes(fake_channel_spike_index{chanIndex});
        fc_amplitudes = amplitudes(:,fake_channel_spike_index{chanIndex});
        fc_ranks = fake_channel_spike_rank{chanIndex};
        save(fn,'fc_spiketimes','fc_amplitudes','fc_ranks','residual_key');
        sorthead(fileIndex).numofsnips(chanIndex) = length(fc_spiketimes);
      end
      fprintf('done\n');
    end
    % Save the overview file with the updated info
    save_overview_file(savedir);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Do the sorting in each fake channel. This simply results in a set of
  % projection directions, landmarks, and (preliminary) landmark grouping
  % decisions.
  % The spikes in each fake channel determine how the landmarks are
  % chosen and grouped, but the final assignment of each event to a
  % cluster will be made after the user has a chance to go over the
  % landmark grouping to see what's worth keeping, what needs to be
  % merged further, etc.  (This manual aspect of sorting will be done in
  % CASS.)
  %
  templates_mtrx = cat(2,componentFileInfo.templates{:});
  rand('state', sum(100*clock));  % generate _different_ keys on each matlab startup
  everRequiredUserInteraction = false;
  for chanIndex = 1:n_fake_channels
  %for chanIndex = 1:4
    fullpath = [savedir filesep 'chan' num2str(chanIndex)];
    if ~exist(fullfile(fullpath,'autosort_info.mat'),'file')
      fprintf('\nSub-clustering channel %d out of %d\n',chanIndex,n_fake_channels);
      % Determine the fraction of snippets to use
      n_spikes = zeros(1,n_files);
      for fileIndex = 1:n_files
        n_spikes(fileIndex) = sorthead(fileIndex).numofsnips(chanIndex);
      end
      n_spikes_tot = sum(n_spikes);
      stride = ceil(n_spikes_tot/options.max_snips_to_cluster);
      fprintf('  Loading snippets\n');
      % Load the snippets
      snips = cell(1,n_files);
      ranks = cell(1,n_files);
      for fileIndex = 1:n_files
        load(fullfile(fullpath,['time_amp' num2str(fileIndex)]));
        if isfield(sorthead,'fh_residual')
          [fh,requiredUserInteraction] = fixpath(sorthead(fileIndex).fh_residual,{'.'});
          if requiredUserInteraction && ~isequal(fh,-1);
            sorthead(fileIndex).fh_residual = fh;
            everRequiredUserInteraction = true;
          end
          mmres = merecmm(fullfile(fh),'tovolts',true);
          header = mmres.header;
          if ~isequal(str2double(key2value(header.wholeheader,'identifier key')),...
              residual_key)
            error(['The residual.merec file key does not match the' ...
              ' time_amp file''s residual_key, indicating that' ...
              ' these two files are probably not a valid pair']);
          end
          snips{fileIndex} = fetch_snippets_from_residual(...
            fc_amplitudes(:,1:stride:end),...
            templates_mtrx,...
            fc_spiketimes(1:stride:end),...
            mmres,...
            componentFileInfo.channels,...
            componentFileInfo.sniprange);
        else
          snips{fileIndex} = templates_mtrx * fc_amplitudes(:,1:stride:end);
        end
        ranks{fileIndex} = fc_ranks(1:stride:end);
      end
      snips = cat(2,snips{:});
      ranks = cat(2,ranks{:});
      n_snips = size(snips,2);
      % Find projection directions by LDA. We'll get the first half of the
      % dimensions using just the rank 1 objects, and the other half will
      % come from including all objects.
      fprintf('  Calculating projection directions\n');
      %pd = calc_pd_from_snips(snips,4,8,options);
      n_dims = ceil(options.n_dimensions_per_channel/2);
      n_preclusters = 2*n_dims;
      pd1 = calc_pd_from_snips(snips(:,ranks == 1),n_dims,n_preclusters,options);
      pd2 = calc_pd_from_snips(snips,n_dims,n_preclusters,options);
      pd3 = lda(snips,ranks);
      pd3 = pd3(:,1:n_dims);
      pd = [pd1 pd2 pd3];
      % Project the spikes
      proj = pd'*snips;
      % Do the "real" clustering: first find landmarks
      fprintf('  Choosing landmarks\n');
      n_landmarks = min(options.n_landmarks_channel,ceil(n_snips/10));
      lminfo = choose_landmarks(proj,n_landmarks);
      % Calculate the landmark centers in terms of spike waveform
      landmarkWaveform = zeros(size(snips,1),n_landmarks);
      for lmIndex = 1:n_landmarks
        landmarkWaveform(:,lmIndex) = ...
          mean(snips(:,lminfo.landmarkList{lmIndex}),2);
      end
      % Now do the clustering
      fprintf('  Clustering with %d landmarks: ',n_landmarks);
      clust = msams(proj,lminfo,struct('variable_metric',false, ...
        'show_progress',true));
      % Save the results
      proj_key = round(1e8*rand(1,1));
      sort_info = struct('channel',chanIndex,...
        'sniprange',componentFileInfo.sniprange,...
        'use_projection',true,...
        't2V',0,...
        'landmarkWaveform',landmarkWaveform,...
        'landmarkT',zeros(1,n_landmarks),...
        'landmarkClust',clust,...
        'projectDirections',pd,...
        'projectKey',proj_key);
      % Save the results
      save(fullfile(fullpath,'autosort_info.mat'),'sort_info');
    end
  end
  if everRequiredUserInteraction
    save_overview_file(savedir);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Pre-calculate the projections of the spike waveforms onto the
  % projection directions. This simply saves interactive time (the
  % projections will usually be low-dimensional, but the full waveforms
  % can be very high-dimensional, especially if there are many
  % electrodes, and the loading time is therefore significant).
  %
  everRequiredUserInteraction = false;
  for chanIndex = 1:n_fake_channels
    fullpath = [savedir filesep 'chan' num2str(chanIndex)];
    if ~exist(fullfile(fullpath,'proj.mat'),'file')
      fprintf('Calculating projections on channel %d out of %d...\n',chanIndex,n_fake_channels);
      tt = load(fullfile(fullpath,'autosort_info'));
      pd = tt.sort_info.projectDirections';
      projectKey = tt.sort_info.projectKey;
      snipProj = cell(1,n_files);
      snipmm = cell(1,n_files);
      for fileIndex = 1:n_files
        load(fullfile(fullpath,['time_amp' num2str(fileIndex)]));
        n_spikes = size(fc_amplitudes,2);
        using_residual = false;
        if isfield(sorthead,'fh_residual')
          [fh,requiredUserInteraction] = fixpath(sorthead(fileIndex).fh_residual,{'.'});
          if requiredUserInteraction && ~isequal(fh,-1);
            sorthead(fileIndex).fh_residual = fh;
            everRequiredUserInteraction = true;
          end
          mmres = merecmm(fullfile(fh),'tovolts',true);
          header = mmres.header;
          if ~isequal(str2double(key2value(header.wholeheader,'identifier key')),...
              residual_key)
            error(['The residual.merec file key does not match the' ...
              ' time_amp file''s residual_key, indicating that' ...
              ' these two files are probably not a valid pair']);
          end
          using_residual = true;
        end
        offset = 0;
        blocksize = 1000;
        while (offset < n_spikes)
          rng = offset+1:min(offset+blocksize,length(fc_spiketimes));
          offset = offset+blocksize;
          if using_residual
          snips = fetch_snippets_from_residual(...
            fc_amplitudes(:,rng),...
            templates_mtrx,...
            fc_spiketimes(rng),...
            mmres,...
            componentFileInfo.channels,...
            componentFileInfo.sniprange);
          else
            snips = templates_mtrx * fc_amplitudes(:,rng);
          end
          proj = pd * snips;
          snipProj{fileIndex} = [snipProj{fileIndex} proj];
          snipmm{fileIndex} = [snipmm{fileIndex},[min(snips); max(snips)]];
        end  % loop over blocks
      end % loop over files
      save(fullfile(fullpath,'proj'),'snipProj','snipmm','projectKey');
    end % if proj.mat exists
  end % loop over fake channels
  if everRequiredUserInteraction
    save_overview_file(savedir);
  end

  function save_overview_file(ovwfsavedir)
    overviewFileName = [ovwfsavedir filesep 'overview'];
    save(overviewFileName,'options','sorthead','lminfo','landmarkClust','componentFile');
  end
end

function pd = calc_pd_from_snips(snips,n_dims,n_preclusters,options)
  n_snips = size(snips,2);
  pdstride = ceil(n_snips/options.n_snips_projection);
  pdsnips = snips(:,1:pdstride:end);
  clustpd = kmeans_hard(pdsnips,n_preclusters);
  pd = lda(pdsnips,clustpd);
  pd = pd(:,1:n_dims);
end