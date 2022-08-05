function cluster_events_raw_manual(filename)
  load('-mat',filename);
  groups = agglabel(clust);
  n_groups = length(groups);
  
  final_groups = {};
  eventMean = zeros(size(eventWave,1),0);
  for groupIndex = 1:n_groups
    indx = {};
    evWtmp = eventWave(:,groups{groupIndex});
    pd = pca(evWtmp');
    %clusttmp = kmeans_hard(evWtmp,3);
    %pd = lda(evWtmp,clusttmp);
    proj = pd(:,1:2)'*evWtmp;
%     if (max(clust2{groupIndex}) > 1)
%       evWtmp = eventWave(:,groups{groupIndex});
%       pd = lda(evWtmp,clust2{groupIndex});
%       if (size(pd,2) > 1)
%         proj = pd(:,1:2)'*evWtmp;
%       else
%         proj = [pd'*evWtmp; randn(1,size(evWtmp,2))];
%       end
    subClustOps.fignum = figure;
    title([num2str(groupIndex) ' out of ' num2str(n_groups)]);
    [indx,aborted] = mdexplore2(proj,mat2cell(evWtmp,size(evWtmp,1),ones(1,size(evWtmp,2))),subClustOps);
    if aborted
      return
    end
    close(subClustOps.fignum);
%     end
    if isempty(indx)
      % We won't toss clusters at this stage; if user just quit, then select all. 
      % This also handles the case where there was only one group.
      indx = {1:length(groups{groupIndex})};
    end
    for newgroupIndex = 1:length(indx)
      final_groups{end+1} = groups{groupIndex}(indx{newgroupIndex});
      eventMean(:,end+1) = mean(eventWave(:,final_groups{end}),2);
    end
  end
  
  % Sort them so that the "best" clusters come out first
  mxEventMean = max(abs(eventMean));
  [smxEM,sortIndex] = sort(mxEventMean,2,'descend');
  final_groups = final_groups(sortIndex);
  eventMean = eventMean(:,sortIndex);

  % Now we convert the indices into eventWave into real snippet times. At
  % this point, snippet times are the only thing that matters.
  sniprange = options.sniprange; % [-10 30];
  for idxCluster=1:length(final_groups)
    tclust = sort(final_groups{idxCluster});
    % NOTE: tclust is column indices to eventWave and eventTime
    spikeTimes=eventTime(tclust); % NOTE: spikeTimes is col indices to allWave

    % remove points near the boundary:
    spikeTimes = spikeTimes(spikeTimes > -sniprange(1)*2+1 & spikeTimes < header.nscans-sniprange(2)*2);
    % tclust = tclust(tclust > -sniprange(1)+1 & tclust < size(eventWave,2)-sniprange(2));

    % segment spikeTimes (NOT tclust)
    spikeTimes=segment_indices(spikeTimes); % continuous events are treated as one event
    spikeTimes=spikeTimes(:,1); % use the ealiest time to represent the event

    % NOTE: to save memory, don't extract/save snippets
    % snippets{idxCluster}=snip; % snippets{idxCluster} is all snippets for cluster_idxCluster
    sniptimes{idxCluster}=spikeTimes; % sniptimes{idxCluster} is the times of all snippets in cluster_idxCluster
  end % for, each cluster

  if(isfield(options, 'fileToSave'))
    fileToSave=options.fileToSave;
  else
    fileToSave=replace_extension(merecFilename, '.raw_cluster');
  end

  if(~isempty(options.saveToDir))
    fileToSave=replace_parent_dir(fileToSave, options.saveToDir);
  end
  % TODO: may need save the median and/or mean for padding purpose
  save(fileToSave, ... % 'snippets',
    'channels', 'merecFilename', 'options', ...
    'labels', 'berryLabels', ...
    'sniptimes', 'eventMean', ...
    'thresh', 'medv', ...
    '-mat');
