function cluster_events_raw_auto(merecFilename, options)
% cluster_events_raw_auto(merecFilename, options)
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
% This function expects to you then run "cluster_events_raw_manual" so that
% you can help determine whether further refinement of raw clusters is
% possible.  Then, you can run "cluster_events_fine" to make use of the
% temporal information.
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
%       saveToDir:  the name of the directory to which the results should
%                   be saved
%       fileToSave: if this field exists, save the result to this file,
%                   otherwise save it as the_input_main_name.raw_cluster_auto.
% NOTE: options has more options than listed here. I only list what the end
%       user may care to know.
%
% See also: CLUSTER_EVENTS_RAW_MANUAL, CLUSTER_EVENTS_RAW.

% HISTORY: 
%       ?           (ZG)    wrote it
%     2007-07-19    (TEH)   Adapted from cluster_events_raw: switch to
%                           msams, save eventWave for use in
%                           cluster_events_raw_manual

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
   fprintf('Finding events...');
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
         % Subtract the bias
         allWave = wave - repmat(medv,1,size(wave,2));
         % Find the timing of events
         [eventTime,peakVal] = findpeaks_multichan(allWave,thresh,options);
         if options.peakOverRange
           % extract the peak values on all electrodes over a time span
	   % Note: this is probably a bad idea
           tRange = -options.adjacentNScans:options.adjacentNScans;
           tEdge = [1 size(wave,2)] + [1 -1]*options.adjacentNScans;
           keepFlag = (eventTime >= tEdge(1)) & (eventTime <= tEdge(2));
           eventTime = eventTime(keepFlag);
           peakVal = peakVal(keepFlag);
           nEvents = length(eventTime);
           eventWave = zeros(size(allWave,1),nEvents);
           for eventIndex = 1:nEvents
             thisRange = tRange+eventTime(eventIndex);
             if peakVal(eventIndex) < 0
               eventWave(:,eventIndex) = min(allWave(:,thisRange),[],2);
             else
               eventWave(:,eventIndex) = max(allWave(:,thisRange),[],2);
             end
           end
         else
           eventWave = allWave(:,eventTime);
         end
         % Store the events
         eventTimes{end+1}=eventTime+scanNumFrom-1;
         eventWaves{end+1}=eventWave;
      end % for each block
   end % for, each scan range to cluster
   
   % concatenate eventTimes and eventWave
   eventTime=cat(2, eventTimes{:});
   eventWave=double(cat(2, eventWaves{:}));
   fprintf('found %d events.\n',length(eventTime));
   
   % to make debug easier:
   seed = 1; % 931316785; 
   rand('seed',seed); randn('seed',seed);

   
   % NOTE: done w/ prepdata.m
   
   % Do the "big" clustering
   n_landmarks = 5000;
   fprintf('Choosing %d landmarks\n',n_landmarks);
   lminfo = choose_landmarks(eventWave,n_landmarks);
   fprintf('Performing automatic clustering\n');
   % Note when using MEX file it doesn't show progress
   clustlm = msams(eventWave,lminfo,struct('show_progress',true,'variable_metric',false));
   clust = clustlm(lminfo.landmarkAssignment); % This is the raw-cluster # assigned to each event
   
   % Now go through and sub-cluster each cluster---this will be useful for
   % doing LDA on the individual clusters, so that they can be more easily
   % refined
   groups = agglabel(clust);
   n_groups = length(groups);
   for groupIndex = 1:n_groups
     fprintf('Sub-clustering group %d (%d events) out of %d...\n',groupIndex,length(groups{groupIndex}),n_groups);
     clust2{groupIndex} = msams(eventWave(:,groups{groupIndex}),struct('variable_metric',false)); % Don't use landmarks this time
   end

   %[landmarkClusters, landmarkFinalPos, landmarks, landmarkPos]=clustdata(eventWave, struct('nLandmarks', 5000));
   % [cclust,         cfinal,           clabel,    c          ]=clustdata(eventWave);
   
   % check_temporal(landmarkClusters, landmarkFinalPos, landmarks, landmarkPos)
   % check_temporal(cclust,           cfinal,           clabel,    c);

   % clust2wf(cclust,           labels, berryLabels, clabel,    ts,         c,           vs)
   % clust2wf(landmarkClusters, labels, berryLabels, landmarks, eventTime,  landmarkPos, allWave)

   
   if(isfield(options, 'fileToSave'))
      fileToSave=options.fileToSave;
   else
      fileToSave=replace_extension(merecFilename, '.raw_auto');
   end
   
   if(~isempty(options.saveToDir))
      fileToSave=replace_parent_dir(fileToSave, options.saveToDir);
   end
   % TODO: may need save the median and/or mean for padding purpose
   save(fileToSave, ... % 'snippets', 
     'channels', 'merecFilename', 'options', ...
     'header',...
     'labels', 'berryLabels', ...
     'thresh', 'medv', ...
     'eventWave','eventTime',...
     'clust','clust2',...    
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
   options = default(options,'adjacentNScans',5);
   options = default(options,'peakOverRange',false);
   options = default(options,'polarity',-1);
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
 