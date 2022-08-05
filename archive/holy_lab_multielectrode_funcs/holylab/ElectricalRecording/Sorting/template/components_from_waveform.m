function [componentsm,options] = components_from_waveform(v,options)
% COMPONENTS_FROM_WAVEFORM: pick a set of spike waveforms from raw data
% Syntax:
%   components = components_from_waveform(v,options)
% where
%   v is an n_channels-by-nscans matrix (usually of type 'single')
%   options must have the following 2 fields:
%     sniprange (e.g., [-7 50]): the displacement from the peak of a
%       spike to its [beginning end];
%     polarity: -1, 1, or 0, for downward, upward, or either.
%   options can also have some other fields that tweak the performance
%   of this algorithm; see documentation in FINDPEAKS_MULTICHAN and this
%   code.
%
% See also: COMPONENTS_FROM_MEREC.
  
% Copyright 2008 by Timothy E. Holy
  
  %% Input parsing & initialization
  [nchan,N] = size(v);
  if (nargin < 2)
    options = struct;
  end
  required_fieldnames = {'sniprange','polarity'};
  missing_fieldnames = setdiff(required_fieldnames,fieldnames(options));
  if ~isempty(missing_fieldnames)
    fprintf('Missing fields:\n');
    for i = 1:length(missing_fieldnames)
      fprintf('%s\n',missing_fieldnames{i});
    end
    error('Must have all required fieldnames');
  end
  %options = default(options,'max_snips',Inf,'isMergeAdjacent',true,'adjacentNScans',5,'thresh_factor',5,'skip_baselineshift',false,'min_per_group',5,'nIters',3,'maxPCA',4000,'nGroupPCA',2,'maxChanPCA',4,'badFitFactor',4);
  options = default(options,'max_snips',Inf,'isMergeAdjacent',true,'adjacentNScans',5,'thresh_factor',5,'skip_baselineshift',false,'min_per_group',5,'nIters',3,'badFitFactor',4);
  
  % Get handles to a variety of utilities for spike sorting
  fu = fit_utilities;

  %% Compute some statistics on the input waveform
  if ~options.skip_baselineshift
    vmed = median(v,2);
    v = v - repmat(vmed,1,N);
  end
  vnoise = mean(abs(v),2); % robust estimation of noise
  
  if ~isfield(options,'thresh')
    % Compute thresholds
    options.thresh = options.thresh_factor*vnoise;
  end
  
  %% Find spike times
  fprintf('Finding spike times...');
  t = findpeaks_multichan(v,options.thresh,options);
  % Toss any times that are too close to edges
  keepFlag = fu.flag_times_in_range(t,[1 N],options.sniprange);
  t = t(keepFlag);
  % Choose a subset of times, if there are too many
  if (length(t) > options.max_snips)
    skip = ceil(options.max_snips/length(t));
    t = t(1:skip:end);
  end
  n_snips = length(t);
  rng = options.sniprange(1):options.sniprange(2);
  fprintf('done\n');
  
  %% Check snippets to make sure none have NaNs
  % NaNs could be used to encode undefined regions/boundaries
  visnan = isnan(v);
  if any(visnan(:))
    fprintf('Validating spikes...');
    badFlag = false(1,n_snips);
    for i = 1:n_snips
      nansnip = visnan(:,t(i) + rng);
      if any(nansnip(:))
        badFlag(i) = true;
      end
    end
    t(badFlag) = [];
    fprintf('done');
  end
  
  n_snips = length(t);
  tshift = t + options.sniprange(1);
  fprintf('%d spike snippets were cut.\n',n_snips);
  
  %% Do clustering on peak amplitudes
  eventPeak = double(v(:,t));
  n_landmarks = 20000;
  fprintf('Performing automatic clustering on peak amplitudes:\n');
  clust = fu.run_clustering(eventPeak,n_landmarks);
  groups = pickgroups(clust,options);
  pA = zeros(nchan,length(groups));
  for i = 1:length(groups)
    pA(:,i) = mean(eventPeak(:,groups{i}),2);
  end
  keepFlag = validate_cluster_amplitude(pA,options.thresh);
  groups = groups(keepFlag);
  if any(~keepFlag)
    fprintf('%d clusters discarded based on amplitude.\n',sum(~keepFlag));
  end

  %% Snip out the waveforms associated with each cluster
  [templates,templates_dev] = find_templates_from_v(v,t,options.sniprange,groups);
  
  %% Fit the waveform and extract min/max pairs for each event
  fprintf('Fitting the waveform with first-pass templates:\n');
  a = fu.fit_waveform(v,tshift,templates,true);
  
  %% Calculate the residual
  fprintf('Calculating the residual...');
  res = fu.calculate_residual(v,tshift,templates,a);
  fprintf('done\n');
    
  %% Get min/max data on each snippet
  fprintf('Getting min/max pairs...');
  pp = fu.minmax_using_residual(res,tshift,templates,a);
  fprintf('done\n');
  
  %% Cluster using min/max values
  % This should be more accurate than the first peak-based clustering
  fprintf('Clustering min/max pairs...');
  clust = fu.run_clustering(double(fu.snip2vec_by_c(pp)),n_landmarks);
  groups = pickgroups(clust,options);
  fprintf('done\n');
  
  %% Snip out the waveforms associated with each cluster
  templates = find_templates_from_resid(res,t,options.sniprange,templates,a,groups);
  pA = squeeze(templates(:,-options.sniprange(1)+1,:));
  keepFlag = validate_cluster_amplitude(pA,options.thresh);
  templates = templates(:,:,keepFlag);
  if any(~keepFlag)
    fprintf('%d clusters discarded based on amplitude.\n',sum(~keepFlag));
  end

  %% Iteratively improve the templates, using cases where the residual has
  %% large error to find new templates
  % Note that we deliberately do not allow the introduction of new event
  % times (which is done in fit_components); the idea is that we want to be
  % certain that our templates represent "real" events, and the best way to
  % do that is to build them from well-isolated spikes.
  for iter = 1:options.nIters
    % Fit the waveform with these templates
    fprintf('Fitting the waveform with improved templates:\n');
%    try
      % Put this in a try-catch block in case we run into memory troubles
      a = fu.fit_waveform(v,tshift,templates,true);
%    catch
%      break
%    end

    % Calculate the residual
    fprintf('Calculating the residual...');
    res = fu.calculate_residual(v,tshift,templates,a);
    fprintf('done\n');
    
    % Calculate the error
    fprintf('Calculating the error...');
    err = fu.residual_errors(res,t,options.sniprange);
    fprintf('done\n');
    
    % Find the regions with significant fitting error
    mederr = median(err);
    derr = err - mederr;
    aderr = mean(abs(derr));
    errthresh = aderr*options.badFitFactor;
    badFlag = derr > errthresh;
    
    % Get the min/max voltages per channel for these regions
    fprintf('Getting data on %d regions with high error...',sum(badFlag));
    pp = fu.minmax_using_residual(res,tshift(badFlag),templates,a(:,badFlag));
    fprintf('done\n');
    
    % Re-cluster. Do this in a way that tends to not accept clusters with
    % few points.
    fprintf('Reclustering these regions...');
    newclust = fu.run_clustering(double(fu.snip2vec_by_c(pp)),n_landmarks,struct('min_to_check',10));
    fprintf('done (%d new clusters found)\n',max(newclust));
    newops = options;
    newops.min_per_group = 0;
    groups = pickgroups(newclust,newops);
    % Check the self-consistency of the snippets contributing to each
    % template by collecting the stddev
    [newtemplates,newerr] = find_templates_from_resid(res,t(badFlag),options.sniprange,templates,a(:,badFlag),groups);
    keepFlag = newerr - mederr < 3*errthresh; % keep only ones that are fairly self-consistent
    % Also check that none of these are too small in their amplitude
    pA = squeeze(newtemplates(:,-options.sniprange(1)+1,:));
    keepFlag = keepFlag & validate_cluster_amplitude(pA,options.thresh);
    n_keep = sum(keepFlag);
    if n_keep > 0
      templates(:,:,end+1:end+sum(keepFlag)) = newtemplates(:,:,keepFlag);
      fprintf('Accepted %d new templates on the basis of their self-consistency.\n',sum(keepFlag));
    else
      fprintf('No additional good templates found, quitting.\n');
      break
    end
  end
  componentsm = templates;
end


function [templates,templates_dev] = find_templates_from_v(v,t,sniprange,groups)
  nchan = size(v,1);
  n_groups = length(groups);
  sniplen = diff(sniprange)+1;
  fprintf('Collecting templates for %d groups...',n_groups);
  templates = zeros(nchan,sniplen,n_groups);
  templates_dev = templates;
  fu = fit_utilities;
  for groupIndex = 1:n_groups
    thisGroup = groups{groupIndex};
    snip = fu.collect_snippets(v,t(thisGroup),sniprange);
    % Use the median so we are not influenced as much by overlapping spikes
    templates(:,:,groupIndex) = median(snip,3);
    dsnip = snip - repmat(templates(:,:,groupIndex),[1 1 size(snip,3)]);
    templates_dev(:,:,groupIndex) = mean(abs(dsnip),3);
  end
  fprintf('done\n');
end

function [templates,stddev] = find_templates_from_resid(res,t,sniprange,components,a,groups)
  nchan = size(res,1);
  n_groups = length(groups);
  sniplen = diff(sniprange)+1;
  fprintf('Collecting templates for %d groups...',n_groups);
  templates = zeros(nchan,sniplen,n_groups);
  stddev = zeros(1,n_groups);
  fu = fit_utilities;
  for groupIndex = 1:n_groups
    thisGroup = groups{groupIndex};
    snip = fu.reconstruct_snippets_using_residual(res,t(thisGroup)+sniprange(1),components,a(:,thisGroup));
    meansnip = mean(snip,3);
    templates(:,:,groupIndex) = meansnip;
    if (nargout > 1)
      n_snips = length(thisGroup);
      dtmp = snip - repmat(meansnip,[1 1 n_snips]);
      stddev(groupIndex) = sqrt(sum(dtmp(:).^2)/n_snips);
    end
  end
  fprintf('done\n');
end

function groups = pickgroups(clust,options)
  [groups,n_per_group] = agglabel(clust);
  groups = groups(n_per_group > options.min_per_group);
  n_groups = length(groups);
  n_discarded = length(n_per_group)-n_groups;
  if (n_discarded > 0)
    fprintf('%d clusters found (%d others discarded for having too few members).\n',n_groups,n_discarded);
  else
    fprintf('%d clusters found.\n',n_groups);
  end    
end

function keepFlag = validate_cluster_amplitude(pA,thresh)
  keepFlag = any(abs(pA) > repmat(thresh,1,size(pA,2)),1);
end

%   %% Find the channels with the most information about each group.
%   % These are loosely defined as the channels with the largest mean peak
%   % amplitudes or largest variance in peak amplitudes. We'll use this to
%   % reduce the PCA to only a few channels.
%   pA = zeros(nchan,n_groups);  % will hold mean peak amplitude
%   rA = zeros(nchan,n_groups);  % will hold the range of peak amplitudes
%   sA = zeros(nchan,n_groups);  % will hold the std. dev. of peak amps.
%   for i = 1:n_groups
%     pA(:,i) = mean(eventPeak(:,groups{i}),2);
%     rA(:,i) = range(eventPeak(:,groups{i}),2);
%     sA(:,i) = std(eventPeak(:,groups{i}),[],2);
%   end
%   threshrep = repmat(abs(options.thresh(:)),1,n_groups);
%   suprathresh = abs(pA) > threshrep | rA > 2*threshrep;
%   no_suprathresh = sum(suprathresh,1) == 0;
%   if any(no_suprathresh)
%     % There were no channels above threshold, choose the biggest
%     nstIndex = find(no_suprathresh);
%     for i = nstIndex
%       [mpA,mpAIndex] = max(pA(:,i));
%       suprathresh(mpAIndex,i) = true;
%     end
%   end
%   toomany_suprathresh = (sum(suprathresh) > options.maxChanPCA);
%   if any(toomany_suprathresh)
%     % Choose the top few channels
%     stIndex = find(toomany_suprathresh);
%     for i = stIndex
%       [ssA,sortIndex] = sort(sA(:,i),'descend');
%       suprathresh(:,i) = false;
%       suprathresh(sortIndex(1:options.maxChanPCA),i) = true;
%     end
%   end
%   
%       % For each group, re-snippet and do PCA on the resulting waveforms
%     % (this should approximately yield spike waveforms with no overlaps)
%     subGroups = {};
%     for groupIndex = 1:n_groups
%       fprintf('Resnippeting group %d...',groupIndex);
%       thisGroup = groups{groupIndex};
%       n_snips = length(thisGroup);
%       thisChan = suprathresh(:,groupIndex);
%       thisNChan = sum(thisChan);
%       snip = zeros(thisNChan,sniplen,n_snips);
%       for snipIndex = 1:n_snips
%         thisT = t(thisGroup(snipIndex));
%         snip(:,:,snipIndex) = res(thisChan,thisT + rng);
%       end
%       snip = snip2vec(snip) + snip2vec(templates(thisChan,:,:)) * a(:,thisGroup);
%       fprintf('done\n');
%       % Do the PCA and keep the first n directions
%       skip = floor(length(thisGroup)/options.maxPCA);
%       if (skip < 1)
%         skip = 1;
%       end
%       snipPCA = snip(:,1:skip:end);
%       pd = pca(snipPCA');
%       pd = pd(:,1:options.nGroupPCA);
%       proj = pd'*snip;
%       % Subcluster
%       fprintf('Subclustering %d snippets:\n',n_snips);
%       lminfo = choose_landmarks(proj,n_landmarks);
%       % Note when using MEX file it doesn't show progress
%       clustlm = msams(proj,lminfo,struct('show_progress',true));
%       clust = clustlm(lminfo.landmarkAssignment);
%       [subg,n_per_subg] = agglabel(clust);
%       subg = subg(n_per_subg > options.min_per_group);
%       n_subg = length(subg);
%       for i = 1:n_subg
%         subg{i} = thisGroup(subg{i});
%       end
%       fprintf('%d sub-clusters found.\n',n_subg);
%       subGroups(end+1:end+n_subg) = subg;
%     end
% 
%     templates = findtemplates(v,t,options.sniprange,subGroups);