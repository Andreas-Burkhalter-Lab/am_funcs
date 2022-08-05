function stim_spike(ephysin,channel,tags,options)
% STIM_SPIKE: analyze stimulus responses as fcn of spike shape w/o sorting
% Syntax:
%   stim_spike(ephysin,channel,tags)
%   stim_spike(ephysin,channel,tags,options)
% where
%   ephysin is a tagged ephys structure array, broken down into trials,
%     and provided in temporal order;
%   channel is the channel number you're interested in analyzing;
%   tags is a cell array of tag strings that label the responses you want
%     to examine, in the sequence that you want to see the response;
%   options is a structure which may have the following fields:
%     n_landmarks (default 20): the number of landmarks to use in
%       splitting up the space of spike waveforms.
% The output is a set of figures, one per cycle, in which a 2-d
% scatterplot of spike waveform projections is overlain with a set of
% graphs, one per landmark (and centered at the landmark position). Each
% graph plots the number of spikes associated with that landmark as a
% function of the stimulus # (the tag index).
% 
% See also: AUTOSORT, CASS.
  
% Copyright 2005 by Timothy E. Holy
  
  %channel_index = find(ephysin.channel == channel);
  %if isempty(channel_index)
  %  error('Channel %d not found\n',channel);
  %end
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'n_landmarks')
    options.n_landmarks = 20;
  end
  
  ephyschan = ephyssubchan(ephysin,channel);
  for i = 1:length(tags)
    trials{i} = strmatch(tags{i},{ephyschan.tag},'exact');
    n_trials(i) = length(trials{i});
  end
  % Organize the data into cycles
  all_trials = sort(cat(1,trials{:}));
  cycle_index = organize_by_cycle(ephyschan(all_trials));
  total_trials = length(all_trials);
  % Load the snippets
  for i = 1:total_trials
    etmp = ephysfetch(ephyschan(all_trials(i)),'snippets');
    snp{i} = etmp.snippets{1};
    %t_snp{i} = (ephyschan(i).sniptimes{1} - ephyschan(i).scanrange(1)) / ...
    %    ephyschan(i).scanrate;
    n_snp(i) = size(snp{i},2);
  end
  snp_all = cat(2,snp{:});
  %t_snp_all = cat(2,t_snp{:});
  [pd,sv] = pca(snp_all');  % Find the projection directions by PCA
  pd = pd(:,1:2);
  proj = pd'*snp_all;
  [idx,cntr,rmsd] = kmeans_hard(proj,options.n_landmarks); % Divide the space using landmarks
  n_landmarks = size(cntr,2);  % In case kmeans tosses landmarks
  % For each trial, compute the number of spikes assigned to each
  % landmark
  indx_end = cumsum(n_snp);
  indx_start = [0 indx_end(1:end-1)]+1;
  for i = 1:length(snp)
    [clabel,nlabel] = agglabel(idx(indx_start(i):indx_end(i)));
    if (length(nlabel) < n_landmarks)
      nlabel(n_landmarks) = 0;
    end
    n_stim_landmark(i,:) = nlabel;
  end
  % Organize the spike # data into cycles
  for i = 1:length(cycle_index)
    n_stim_landmark_cycle{i} = nan(length(tags),n_landmarks);
    for j = 1:length(cycle_index{i})
      % Go in tag-order, rather than temporal order
      cindex = cycle_index{i}(j);
      tag_index = strmatch(ephyschan(all_trials(cindex)).tag,tags);
      n_stim_landmark_cycle{i}(tag_index,:) = n_stim_landmark(cindex,:);
    end
    % Now plot
    figure('name',['Cycle ' num2str(i)])
    %subplot(1,2,1)
    hax_base = gca;
    % Spike projection scatterplot
    rng = indx_start(cycle_index{i}(1)):indx_end(cycle_index{i}(end));
    hscatter = plot(proj(1,rng),proj(2,rng),'b.','MarkerSize',1);
    % Plot the landmarks
    %line(cntr(1,:),cntr(2,:),'Marker','o','MarkerFaceColor','b', ...
    %     'LineStyle','none');
    axis tight
    xlim = get(gca,'XLim');
    ylim = get(gca,'YLim');
    % Plot the nspikes vs. stimulus in the same physical arrangement as the
    % landmarks
    %hax_base = subplot(1,2,2);
    set(hax_base,'XLim',xlim,'YLim',ylim);
    panelsize = mean(rmsd);
    for landmarkIndex = 1:n_landmarks
      [alimx,alimy] = data2norm(cntr(1,landmarkIndex)+[-1 1]*panelsize,...
        cntr(2,landmarkIndex)+[-1 1]*panelsize,...
        hax_base);
      hax = axes('position',[alimx(1) alimy(1) diff(alimx) diff(alimy)]);
      hline = plot(n_stim_landmark_cycle{i}(:,landmarkIndex),'r.',...
        'LineStyle','-','LineWidth',2);
      axis tight
      ylim = get(gca,'YLim');
      set(gca,'YLim',[0 1.01*ylim(2)])
      set(hax,'Visible','off')
      set(hline,'Visible','on')
    end
  end
    