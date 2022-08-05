function cluster_events_raw(merecFilename, options)
% cluster_events_raw(merecFilename, options)
% Detect events on groups of channels on the array, and split the events up
% into "raw clusters" that contain similar events. Each "event" is
% considered to occur at a single time-slice, roughly corresponding to the
% peak of an action potential; it does not contain any information about
% the temporal shape of action potentials. The main idea is to use the
% cross-electrode distribution of amplitude to make a first cut at
% separating events into groups, the "raw clusters."
%
% Later steps will use these raw clusters to make further refinements based
% on the full waveforms. Having this initial split has the benefit of
% allowing the PCA analysis to "concentrate" on small differences between
% similar waveforms.
%
% NOTE: here the "cluster" means clustering typical data to create templates.
%       Don't confuse it w/ using templates to cluster spikes to cells
%       (fit_template() does that job).
% PRE:
%    merecFilename: a *.merec file to generate raw clusters
%    options: (optional) a struct that has fields:
%       channels: a struct that has 2 fields, 'what' and 'value'. If you
%         want to specify particular channel numbers then set
%              'what'='specific' and 'value'=[channel list]
%         If you want to use the conventional half-fields, then set
%              'what'='half' and 'value'=1 or 2 (depending on which half
%                    you want)
%       blackList: a vector of channels you want them excluded from sorting.
%       timeRangeToCluster: a nRanges-by-2 matrix, where each row contains
%         the [start stop] time (in seconds) for a contiguous range of
%         scans used to look for events. If you leave this out, or specify
%         [0 -1], the entire file will be used.
%       timeRangeToCalcThreshold=[0 3]: the data range in sec used to
%         calculate the threshold on each channel;
%       sniprange=[-16 47]: the snip range
%       fileToSave: if this field exists, save the result to this file,
%                   otherwise save it as the_input_main_name.raw_cluster.
% NOTE: options has more options than listed here. I only list what the end
%       user may care to know.
% HISTORY: 
%       ?           (ZG)    wrote it

   if(nargin==1)
      options=struct;
   end
   options=fillOptions(options);

   merecFilename=make_abs_path(merecFilename);
   
   memm = merecmm(merecFilename,'tovolts',true,'contiguous',true);
   allChannels=memm.channels;
   % channels = setdiff(allChannels(1:30), options.channelsToSkip);

   switch(options.channels.what)
      case 'half'
         if(options.channels.value==1)
            channels=get_hda_chan(1);
         else
            channels=get_hda_chan(2);
         end
      case 'specific'
         channels=options.channels.value;
   end % switch,
   channels=intersect(allChannels, channels); % this will get rid of channels not recorded.
   channels=setdiff(channels, options.blackList);
   % channels get sorted also.
   %channels=sort(channels); 

   % calc corresponding labels
   header=memm.header;
   allLabels=cellstr(split_label(key2value(header.wholeheader, 'label list')));
   [tt, channelIndices]=intersect(allChannels, channels); % NOTE: channels were sorted by setdiff().
   % NOTE: channelIndices is indices to allChannels/allLabels
   labels=allLabels(channelIndices);
   % because channelIndices is indices to allChannels, it is also Berry Lab's
   % notation (NOTE: allChannels is orderd by Berry Lab's notation)
   berryLabels=channelIndices; % TODO: this works only when berrylab layout is used and all channels are recorded.
   % NOTE: the reliable way is to refer channels by labels instead of by
   %       comedi channel numbers or berryLabels.

   ttScanRange=round(header.scanrate*options.timeRangeToCalcThreshold);
   if(ttScanRange(1)<=0) 
      ttScanRange(1)=1; 
   end
   if(ttScanRange(2)>header.nscans)
      ttScanRange(2)=header.nscans;
   end
   wave = memm(channels, ttScanRange); % NOTE: each row is for each channel; each col is for each scan

   % Compute the median and threshold for each channel
   [d,N] = size(wave); % d: #channels; N: #scans
   medv = median(wave,2); % got a column vector. That is, median for each channel
   allWave = wave - repmat(medv,1,N); % allWave: normalized wave
   adv = mean(abs(allWave),2); % got a column vector
   thresh = 6*adv;  % a column vector. That is, threashold for each channel
   
   blockSize=options.blockDuration*header.scanrate;
   
   eventTimes={};
   eventWaves={};
   allScanRangesToCluster=round(options.timeRangeToCluster*header.scanrate); 
   for idxScanRangeToCluster=1:size(allScanRangesToCluster,1)
      scanRangeToCluster=allScanRangesToCluster(idxScanRangeToCluster, :);
   
      scanRangeToCluster(1)=scanRangeToCluster(1)+1; % +1: matlab is 1-base numbering.
      if(scanRangeToCluster(1)>header.nscans)
         error('tried to cluster beyond file end');
      end
      if(scanRangeToCluster(2)<0 || scanRangeToCluster(2)>header.nscans)
         scanRangeToCluster(2)=header.nscans;
      end

      nBlocks=ceil((diff(scanRangeToCluster)+1)/blockSize);
      for idxBlock=1:nBlocks
         scanNumFrom=(idxBlock-1)*blockSize+scanRangeToCluster(1);
         scanNumTo=min(scanNumFrom+blockSize-1, scanRangeToCluster(2));
         wave= memm(channels, [scanNumFrom scanNumTo]);
         [d,N] = size(wave); % d: #channels; N: #scans
         thresh_mtrx = repmat(thresh,1,N);
      
         allWave = wave - repmat(medv,1,N); % allWave: normalized wave for all times
      
         % Find the points in time that exceeded the threshold
         absNormalizedWave=abs(allWave);
         isbig = absNormalizedWave > thresh_mtrx; % just big enough
      
         % 3-pts maxima
         if(options.isUsePeakCriteria)
            isPeak= absNormalizedWave(:,[1 1:end-1]) < absNormalizedWave  &  absNormalizedWave >= absNormalizedWave(:,[2:end end]);
            isbig=isbig & isPeak; % big peak
         end

         % don't use sum() which will convert to double from logical
         if(d>1)
            keep_timeslice = max(isbig); % those times that something happened on at least one electrode
         else
            keep_timeslice = isbig;
         end
         
         % NOTE: an event is a time when >=1 one channel is beyond threshold
      
         eventTime = find(keep_timeslice); % record the times. NOTE: eventTime is the col indices to allWave

         % isMergeAdjacent=0;
         if(options.isMergeAdjacent)
            ibreak = find(diff(eventTime) >= 20); % TODO: make it a param. As a rule of thumb, pick a value > -sniprange(1)+maxshift
            eventTime=eventTime([1 ibreak+1]);
         end

         eventWave=double(allWave(:,eventTime));
      
         eventTimes{end+1}=eventTime+scanNumFrom-1;
         eventWaves{end+1}=eventWave;
      end % for each block
   end % for, each scan range to cluster
   
   % concatenate eventTimes and eventWave
   eventTime=cat(2, eventTimes{:});
   eventWave=cat(2, eventWaves{:});
   
   % to make debug easier:
   seed = 1; % 931316785; 
   rand('seed',seed); randn('seed',seed);

   
   % NOTE: done w/ prepdata.m
   
   [landmarkClusters, landmarkFinalPos, landmarks, landmarkPos]=clustdata(eventWave, struct('nLandmarks', 5000));
   % [cclust,         cfinal,           clabel,    c          ]=clustdata(eventWave);
   
   % check_temporal(landmarkClusters, landmarkFinalPos, landmarks, landmarkPos)
   % check_temporal(cclust,           cfinal,           clabel,    c);

   % clust2wf(cclust,           labels, berryLabels, clabel,    ts,         c,           vs)
   % clust2wf(landmarkClusters, labels, berryLabels, landmarks, eventTime,  landmarkPos, allWave)

   sniprange = options.sniprange; % [-10 30];
   snipMedian=repmat(medv, 1, diff(sniprange)+1);
   for idxCluster=1:length(landmarkClusters)
      tt = landmarks(landmarkClusters{idxCluster}); % events assigned to cluster_idxCluster
      tclust = sort(cat(2,tt{:})); % merge the indices of events assigned to the cluster
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
      fileToSave=replace_extension(merecFilename, '.raw_cluster_auto');
   end
   [pth,basename,ext] = fileparts(fileToSave);
   if isempty(ext)
     ext = '.raw_cluster_auto';
   end
   fileToSave = fullfile(pth,basename,ext);
   
   if(~isempty(options.saveToDir))
      fileToSave=replace_parent_dir(fileToSave, options.saveToDir);
   end
   % TODO: may need save the median and/or mean for padding purpose
   save(fileToSave, ... % 'snippets', 
        'channels', 'merecFilename', 'options', ...
	'landmarkClusters', 'labels', 'berryLabels', 'landmarks', ...
	'sniptimes', 'landmarkPos', 'landmarkFinalPos', ...
	'thresh', 'medv', ...
	'-mat');

     
function options=fillOptions(options)     
   if(~isfield(options, 'channels'))
      options.channels.what='half'; % other valid string: specific
      options.channels.value=1; % first half
   end

   default_options('blackList', []);
   
   if(~isfield(options, 'timeRangeToCalcThreshold'))
      options.timeRangeToCalcThreshold=[0 3]; % TODO: default too short?
   end
   if(~isfield(options, 'blockDuration'))
      options.blockDuration=10; % 10sec
   end
   if(~isfield(options, 'isUsePeakCriteria'))
      options.isUsePeakCriteria=1;
   end
   if(~isfield(options, 'isMergeAdjacent'))
      options.isMergeAdjacent=1; 
   end
   if(~isfield(options, 'sniprange'))
      options.sniprange=[-16 47];
   end
   if(~isfield(options, 'timeRangeToCluster'))
      options.timeRangeToCluster=[0 -1]; % this option specifies which part of file
                                         % you want to use to generate templates
                                         % (instead of the whole file).
                                         % It's nRangesx2 matrix. Each row is a range in second.
                                         % Default to the whole file.
   end
   if(~isfield(options, 'saveToDir'))
      options.saveToDir='';
   end
 