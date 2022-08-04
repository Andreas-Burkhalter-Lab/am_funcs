function sort_info_o = autosort_AM(sorthead,savedir,options)
% AUTOSORT: automatic spike sorting
%%% AM edited 8/10/15 on vivid
%%% -added options to suppress outputs

% To sort data and save to disk, use the following syntax:
%   autosort(sorthead,savedir,options)
% To sort data and return the result (works with a single channel only),
% use the following syntax:
%   sort_info = autosort(sorthead,options)
% where
%   sorthead is a sorting header (see snipfile2sortheader);
%   savedir is the name of the directory to be created to contain the
%     results (if saving);
%   options is a structure which may have the fields described in
%     AUTOSORT_OPTIONS.  In addition, it may have the fields:
%      channels: the list of channels to sort (default: all)
%
% See also: AUTOSORT_OPTIONS, AUTOSORT_CALC_T2V.
  
% Copyright 2005 by Timothy E. Holy
% 2005-10-30 RCH - complete hack job to make null sets within tsnip
%                  come out with dimensions that don't confuse cat
% 2006-03-23 RCH - added capability for reading in timemarkers added 
%                  by merec and automatically creating them as timemarkers
%                  seen by cass (appear automatically, but can be delted if
%                  desired)
% 2006-06-27 RCH - fixed timemarker readin system so that doesn't mess up
%                  when first marker is user defined
  
    %% AM 8/201/15 added options
warn_overwrite = 0; % display warning when trying to mkdir a directory that already exists
display_channel = 1; % print current channel being processed
display_mean_n = 0; % print 'Mean n' values from adaptive_meanshift.m
progressbar_on = 0; % display progress processing current channel
warn_dot_product = 0; %%%%%%%%%%%%%%%%%%%%% not implemented yet

%%

  saving = 1;
  if (nargout > 0)
    saving = 0;
    if (nargin < 2)
      options = struct;
    else
      options = savedir;
    end
  else
    if (nargin < 3)
      options = struct;
    end
  end
  
  if ~isfield(options,'channels')
    options.channels = unique([sorthead.channels]);
  end
  options = autosort_options(options);
  nFiles = length(sorthead);
  
  if (~saving && length(options.channels) > 1)
    error(['If returning results to memory, must specify a single ' ...
           'channel'])
  end
  
  % Sort the files in temporal order of recording
  fileStartTime = sortheader_absolute_starttime(sorthead);
  [fileStartTime,reorderIndex] = sort(fileStartTime);
  sorthead = sorthead(reorderIndex);
  
  % Prepare the output files & directories
  if saving
    % Make the saving directory
    
        %% AM added 8/201/15
    if ~warn_overwrite
        warning('off','MATLAB:MKDIR:DirectoryExists');
    end
    %%
    
    [status,msg] = mkdir(savedir);
    if ~status
      error(['Can''t make directory ' savedir ': ' msg]);
    end
    if ~isempty(msg) && warn_overwrite %% AM edited 8/201/15
      warning(msg);
    end
    % Save the overview file
    overviewFileName = [savedir filesep 'overview'];
    save(overviewFileName,'sorthead','options');
    % Make subdirectories
    for channelIndex = 1:length(options.channels)
      mkdir([savedir filesep 'chan' ...
             num2str(options.channels(channelIndex))]);
    end
    
    warning('on','MATLAB:MKDIR:DirectoryExists'); %% AM added 8/201/15
    
  end

  % Preparatory stuff for progress bar
  n_Rfactor = length(options.Rfactor);
  progressMax = options.n_replicates*n_Rfactor + 3;  

  %
  % Get started: loop over channels
  %
  channelIndexCompleteFlag = zeros(size(options.channels));
  for channelIndex = 1:length(options.channels)
      
     %% AM 8/10/15 added conditional
    if display_channel
        fprintf('Current channel: %d\n',options.channels(channelIndex));
    end
    %%    
    
    % Let the user know that something is being done (progress bar)
        
    %% AM 8/10/15 added conditional
    if progressbar_on
        progstr = sprintf('Clustering progress on channel %d (%d/%d)',...
          options.channels(channelIndex),...
          channelIndex,...
          length(options.channels));
        progress_bar(struct('progress',0,...
          'max',progressMax,...
          'what',progstr));
    end
      %%
    
    % Focus on one channel, and load the spike times in units of scan numbers
    shc = sortheader_importchan(sorthead,options.channels(channelIndex));
    nsnips_per_file = [shc.numofsnips];
    nsnips_total = sum(nsnips_per_file);
    
    %
    % Check to make sure we have enough snippets to proceed
    %
    if (nsnips_total < options.min_snips_to_cluster)
      continue
    end
    channelIndexCompleteFlag(channelIndex) = 1;
    
    %
    % Determine the tradeoff between time and voltage (a measure of
    % recording stability)
    %
    if ~isfield(options,'t2V')
      %[t2V,options] = autosort_calc_t2V(shc,options);
      [t2V,options] = options.t2V_func(shc,options);
    else
      if isscalar(options.t2V)
        t2V = options.t2V;
      else
        t2V = options.t2V(channelIndex);
      end
    end
    
    %% AM 8/10/15 added conditional
    if progressbar_on
        progress_bar(struct('progress',1,...
          'max',progressMax,...
          'what',progstr));
    end
    %%%
    
    %
    % We'll do the clustering on a subset of waveforms. Get them.
    %
    sniprange = reshape([shc.sniprange],2,length(shc));
    sniprange = unique(sniprange','rows');
    if (size(sniprange,1) > 1)
      error(['Snippets do not have the same number of samples in each ' ...
             'file']);
    end
    snipLength = diff(sniprange)+1;
    max_snips_in_mem = min(round(options.max_snip_memsize/snipLength),...
                           options.max_snips_to_cluster);
    subset_fraction = min(1,max_snips_in_mem/nsnips_total);
    stride = round(1/subset_fraction);
    snipIndex = cell(1,nFiles);
    for fileIndex = 1:nFiles
      % Select evenly-spaced spikes
      snipIndex{fileIndex} = 1:stride:nsnips_per_file(fileIndex);
    end
    snip = sortheader_readsnips(shc,snipIndex);
    tsnip = sortheader_absolute_sniptime(shc,snipIndex);
    snipall = cat(2,snip{:});
    snipOrig = snipall;  % save for when we calculate waveforms (projection)
    n = length(tsnip);
    for nth = 1:n
      if size(tsnip{nth},1) == 0
        tsnip{nth} = (tsnip{nth})';
      end
    end
    tsnipall = cat(2,tsnip{:});
    tsnipall = tsnipall-fileStartTime(1);  % Convert to relative to start of
                                           % experiment (just for convenience)
    nsnips = length(tsnipall);
        
    %% AM 8/10/15 added conditional
    if progressbar_on
        progress_bar(struct('progress',2,...
          'max',progressMax,...
          'what',progstr));
    end
    %%


    tOptions = options;   % This is meant to hold options specific to
                          % this channel. Generic options get
                          % stored in options itself.

    %
    % Consider projection/dimensionality reduction
    %
    projectDirections = [];
    if (options.use_projection || options.ploteach || options.plotclust)
      % Either we want to reduce dimensionality (to reduce noise), or
      % we're visualizing things and need to have a way of meaningfully
      % getting down to 2-d.
      projfrac = min(1,options.n_snips_projection/nsnips);
      projstride = round(1/projfrac);
            
    %% AM 8/10/15 added
    options.projection_func = @pca_sort_bootstrap_AM; % allows for suppressing progressbar
    %%
      
      [projectDirections,options] = ... %% AM 8/10/15 added third input
          options.projection_func(snipall(:,1:projstride:end),options,progressbar_on);
      if options.use_projection
        snipall = projectDirections'*snipall;  % project the waveforms
      end
      if (options.ploteach || options.plotclust)
        pdim = min(2,size(projectDirections,2)); % visualize in at most 2-d
        if options.use_projection
          % We've already projected
          tOptions.plotproject = ...
              [eye(pdim) zeros(pdim,size(projectDirections,2)-pdim+1)];
          % The extra [0;0] is for time, i.e., time will not be displayed
        else
          tOptions.plotproject = [projectDirections(:,1:pdim)',...
                              repmat(0,pdim,1)];
        end
      end
    end
    % Append spike time as a sorting coordinate
    snip0 = snipall;
    snipall = [snipall; tsnipall*t2V];
    
            
    %% AM 8/10/15 added conditional
    if progressbar_on
        progress_bar(struct('progress',3,...
                            'max',progressMax,...
                            'what',progstr));
    end
    %%
    
    if (options.ploteach || options.plotclust)
      shg   % Bring a figure to the fore
    end

    % Give the user the chance to adjust t2V
      if options.interactive_landmarking
        landmark_repeat = 1;
        while landmark_repeat
          [clandmarkIndex,clandmarkPosition] = ...
                    options.landmark_func(snipall,100);   %%%%% hard-coded # landmarks
          % This only makes sense if we're projecting, so don't worry
          % about it. Also, it's hard-wired to show t on the x-axis
          [landmark_repeat,t2V] = landmarkgui(clandmarkPosition([end 1],:),...
                                            snipall([end 1],:),...
                                            clandmarkIndex,...
                                            t2V);
          if landmark_repeat
            snipall = [snip0; tsnipall*t2V];
          end
        end
      end


    %
    % Do the clustering
    %
    % User might want to try several times, using different sets of
    % landmarks, so loop over the number of replicates
    landmarkClust = cell(1,options.n_replicates);
    for repeatIndex = 1:options.n_replicates
      % Establish landmarks, and determine how our subset of waveforms
      % should be assigned to landmarks (using K-means)
      [clandmarkIndex,clandmarkPosition] = ...
          options.landmark_func(snipall,options.n_landmarks);
      % Calculate the average spike waveform associated with each
      % landmark
      nLandmarks = size(clandmarkPosition,2);
      clabel = agglabel(clandmarkIndex);
      clandmarkWaveform = zeros(snipLength,nLandmarks);
      for i = 1:nLandmarks
        clandmarkWaveform(:,i) = mean(snipOrig(:,clabel{i}),2);
      end
      % Sort the landmarks in order of peak voltage
      cmaxWaveform = max(abs(clandmarkWaveform));
      [smax,p] = sort(cmaxWaveform,2,'descend');
      clandmarkPosition = clandmarkPosition(:,p);
      clandmarkWaveform = clandmarkWaveform(:,p);
      clandmarkIndex = p(clandmarkIndex);

      landmarkWaveform{repeatIndex} = clandmarkWaveform; % save for posterity

      % Convert the time associated with each landmark
      if t2V
        landmarkT{repeatIndex} = clandmarkPosition(end,:)/t2V;
      else
        landmarkT{repeatIndex} = zeros(1,nLandmarks);
      end
      
      % Flow the landmarks uphill, to figure out how they should group
      % together
      for RIndex = 1:n_Rfactor   % allow a range of smoothing scales
        tOptions.plottitle = sprintf('Channel %d, repeat %d, Rfactor %g',...
                                     options.channels(channelIndex),...
                                     repeatIndex,...
                                     options.Rfactor(RIndex));
        tOptions.R = options.Rfactor(RIndex);    
        
    %% AM 8/10/15 added 
        options.cluster_func = @adaptive_meanshift_AM;
    %%
        
        clandmarkClust = options.cluster_func(snipall,...
                                   clandmarkPosition,...
                                   tOptions,display_mean_n); % AM 8/10/15 added last input
        landmarkClust{repeatIndex}(RIndex,:) = clandmarkClust;
        if options.plotclust
          autosort_plotclust(clandmarkClust,clandmarkIndex,snipall, ...
                             tOptions)
        end
                
    %% AM 8/10/15 added conditional
    if progressbar_on
        progress_bar(struct('progress',3+(repeatIndex-1)*n_Rfactor+RIndex,...
                            'max',progressMax,...
                            'what',progstr));
    end
    %%
      end  % for RIndex = 1:n_Rfactor
    end  % for repeatIndex = 1:options.n_replicates
    
    % In case there are timemarkers already in the .merec files, 
    % take them out of the sorthead, convert them into a form that can
    % be recognized by cass as it opens, and include them in sort_info
    if isfield(sorthead,'timemark')
        count = 1;
        timeMarker = [];
        scanrate = sorthead(1).scanrate;
        for nthFile = 1:nFiles
            marks = sorthead(nthFile).timemark;
            if isfield(marks,'time')
                nMarks = length(marks.time);
                for nthMark = 1:nMarks
                    mark.fileIndex = nthFile;
                    mark.fileTime = scanrate*marks.time(nthMark);
                    mark.clustNum = NaN;
                    mark.action = 'mark_all';
                    mark.comment = {marks.comment{nthMark}};
                    if count == 1;
                        timeMarker = mark;
                    else
                        timeMarker(count) = mark;
                    end
                    count = count+1;
                end
            end
        end
    else
        timeMarker = [];
    end
    
    % We have results! Output to file, if desired
    sort_info = struct('channel',options.channels(channelIndex),...
                       'sniprange',sniprange,...
                       'use_projection',options.use_projection,...
                       't2V',t2V,...
                       'landmarkWaveform',landmarkWaveform,...
                       'landmarkT',landmarkT,...
                       'Rfactor',options.Rfactor,...
                       'landmarkClust',landmarkClust,...
                       'timeMarker',timeMarker);
    if options.use_projection
      for i = 1:length(sort_info)
        sort_info(i).projectDirections = projectDirections;
      end
    end
    if saving
      save_sort_info(sort_info,[savedir filesep ...
        'chan' num2str(sort_info(1).channel)...
        filesep 'autosort_info']);
    end
      
    %
    % Finally, need to assign each waveform in the _whole_ snippet file
    % (not just the sorted subset) to a landmark. But there's a separate
    % function for that: autosort_apply. Do this after the user has had a
    % chance to inspect or choose among the clustering results.
    %
  end  % for channelIndex = 1:length(options.channels)

  if saving
    % Update with any changes in options that were set by called
    % sub-functions (e.g., t2V, projection, etc.) and update the channel
    % list
    options.channels = options.channels(logical(channelIndexCompleteFlag));
    save(overviewFileName,'sorthead','options')
  end
  if (nargout > 0)
    sort_info_o = sort_info;
  end
      

function autosort_plotclust(clust,yclose,x,options)  
  % Map cluster numbers from landmarks to points  
  clust = clust(yclose);
  % Do the plotting
  clf
  co = get(gca,'ColorOrder');
  hold on
  clabel = agglabel(clust);
  nclust = length(clabel);
  for i = 1:nclust
    colindx = mod(i-1,size(co,1))+1;
    clindx = clabel{i};
    xp = options.plotproject*x;
    if (size(xp,1) > 1)
      plot(xp(1,clindx),xp(2,clindx),'.','Color',co(colindx,:));
    else
      plot(xp(clindx),zeros(size(clindx)),'.','Color',co(colindx,:));
    end
  end
  hold off
  axis equal
  if isfield(options,'plottitle')
    title(options.plottitle);
  end
  if (options.plotpause || options.plotclustpause)
    pause
  else
    drawnow
  end
  

