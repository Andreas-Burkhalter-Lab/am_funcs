function cluster_correlations(shc,clust,options)
% shc: channel-specific sortheader 
% clust: a cell array of cluster numbers, one vector/file
% options has the following fields:
%   corr_tmax: compute correlations up to this time shift;
%   corr_tmin: the size of the smallest bin; if you set this equal to 0,
%     it will default to a single scan;
%   n_log_bins: the number of bins, spaced logarithmically in time, to use in
%     displaying correlations
%   selectedIndex: the index numbers of selected clusters
%   timeMarker: if supplied, will be used to restrict analysis to "valid"
%     intervals for each cell (i.e., it will respect "start_cell" and
%     "stop_cell" commands)
%
% History:
%   ?           (TH)    wrote it
%   2007-04-13  (RCH)   added ability to control where figures pop up
%   2007-07-10  (TEH)   added support for timeMarkers
  
  % Split indices into separate clusters
  nfiles = length(shc);
  clustIndex = cell(nfiles,1);
  clustn = zeros(nfiles,1);
  for i = 1:nfiles
    [tIndex,nperlabel] = agglabel(clust{i}+2);  % +2 for delete = -1
    clustIndex(i,1:length(tIndex)) = tIndex;
    clustn(i,1:length(tIndex)) = nperlabel;
  end
  max_clusters = size(clustIndex,2);
  % Implement timeMarker support: discard spike times that fall outside of
  % of the valid intervals specified by the timeMarkers
  if isfield(options,'timeMarker')
    intervals = timemarkers2intervals(options.timeMarker,[shc.nscans],max_clusters-2);
    for fileIndex = 1:nfiles
      for cellIndex = 1:size(intervals,1)
        this_interval = intervals{cellIndex,fileIndex};
        this_clust_index = clustIndex{fileIndex,cellIndex+2}; %+2 for delete=-1
        this_t = shc(fileIndex).sniptimes(this_clust_index)';
        ikeep = false(size(this_clust_index));
        for intervalIndex = 1:size(this_interval,1)
          ikeep = ikeep | (this_t >= this_interval(intervalIndex,1) & ...
            this_t < this_interval(intervalIndex,2));
        end
        clustIndex{fileIndex,cellIndex+2} = this_clust_index(ikeep);
      end
    end
  end
  % Keep only selected, non-empty clusters
  nperclust = sum(clustn,1);
  selectedIndex = intersect(options.selectedIndex+1,find(nperclust)); %+1 for delete
  clustIndex = clustIndex(:,selectedIndex);
  clustn = clustn(:,selectedIndex);
  nclusters = size(clustIndex,2);
  

  % Set up the bin boundaries
  if (options.corr_tmin == 0)
    options.corr_tmin = 1/max([shc.scanrate]);  % set to a single scan
  end
  nbins = options.corr_tmax/options.corr_tmin;
  tbreak = (1:nbins)*options.corr_tmin;
  tcenter = tbreak - options.corr_tmin/2;
  % Map these linearly-spaced bins to log-spaced bins
  tbreak_log = [0 logspace(log10(tbreak(1)),log10(tbreak(end)),options.n_log_bins)];
  logbinIndex = cell(1,options.n_log_bins);
  logbinBreak = nan(1,options.n_log_bins);
  for i = 1:options.n_log_bins
    logbinIndex{i} = find(tcenter >= tbreak_log(i) & ...
                          tcenter < tbreak_log(i+1));
    if ~isempty(logbinIndex{i})
      logbinBreak(i) = tbreak(logbinIndex{i}(end));
    end
  end
  logbinIndex(isnan(logbinBreak)) = [];
  logbinBreak(isnan(logbinBreak)) = [];
  logbinBreak = [0 logbinBreak];
  n_log_bins = length(logbinIndex);

  %
  % Autocorrelations
  %
  hfig_ac = figure;
  if isfield(options,'fig_positions')
      if isstruct(options.fig_positions)
          set(hfig_ac,'Position',options.fig_positions.autocorr);
      end
  end
  setappdata(hfig_ac,'shc',shc);
  setappdata(hfig_ac,'clustIndex',clustIndex);
  %setappdata(hfig_ac,'selectedIndex',selectedIndex);
  n_per_row = 4;
  n_rows = ceil(nclusters/n_per_row);
  n_cols = min(n_per_row,nclusters);
  for i = 1:nclusters
    n_per_bin = zeros(1,nbins);
    for k = 1:nfiles
      tclust = double(shc(k).sniptimes(clustIndex{k,i}))/shc(k).scanrate;
      tclust=sort(tclust);
      n_per_bin = n_per_bin + ...
          autocorrspike(tclust,options.corr_tmax,nbins); 
    end
    rate = zeros(1,n_log_bins);
    rate1 = zeros(1,n_log_bins);
    for j = 1:n_log_bins
      rate1(j) = 1/diff(logbinBreak([j j+1]));
      rate(j) = sum(n_per_bin(logbinIndex{j})).*rate1(j);
    end
    subplot(n_rows,n_cols,i)
    col = unique_color(selectedIndex(i)-1,max_clusters);
    setappdata(gca,'clusterColor',col);
    [x,y] = stairs(logbinBreak,[0 rate]);
    hpatch = patch(x([1:end end]),[y; 0],col,...
      'EdgeColor',col);
    set(gca,'XScale','log','TickDir','out','XLim',x([1 end]))
    % Sometimes it makes bad decisions about tick marks along the X axis.
    % Check for this
    xtick = get(gca,'XTick');
    if (length(xtick) < 2)
      ulx10 = unique(round(log10(x(find(x)))));
      set(gca,'XTick',10.^ulx10);
    end
    ylim = get(gca,'YLim');
    % Draw the 1/bin scale marks
    line(sqrt(logbinBreak(1:end-1).*logbinBreak(2:end)),...
         rate1,...
         'Color','k',...
         'LineStyle','none',...
         'Marker','.');
    setappdata(gca,'cluster',i);
    setappdata(gca,'dataSource',hfig_ac);
    % Draw the selector line
    hline = line(x([4 4]),ylim,'Color','k','LineStyle','--');
    drag_line(hline,struct('type','v','onDragDone',@cc_showwaveforms));
    if options.set_correlation_yaxes_without_fastest
        halfway_pt = ceil(length(y)/2);
        long_half = y(halfway_pt:end);
        max_val = max(long_half);
        ylim(2) = 1.1*max_val;
    end
    set(gca,'YLim',ylim);
  end
  %hfig_cc = figure;


  %
  % Cross-correlations
  %
  if (nclusters > 4 || nclusters < 2)
    return  % Too much to fit on cross-correlation plots
  end
  hfig_cc = figure;
  if isfield(options,'fig_positions')
     if isstruct(options.fig_positions)
       set(hfig_cc,'Position',options.fig_positions.crosscorr);
     end
  end
  setappdata(hfig_cc,'shc',shc);
  setappdata(hfig_cc,'clustIndex',clustIndex);
  %setappdata(hfig_ac,'selectedIndex',selectedIndex);
  n_rows = nclusters;
  n_cols = nclusters;
  for i = 1:nclusters
    for j = i+1:nclusters
      n_per_bin = zeros(1,2*nbins);
      for k = 1:nfiles
        tclust1 = double(shc(k).sniptimes(clustIndex{k,i}))/ ...
                  shc(k).scanrate;
        tclust2 = double(shc(k).sniptimes(clustIndex{k,j}))/shc(k).scanrate;
        tclust1=sort(tclust1);
        tclust2=sort(tclust2);
        n_per_bin = n_per_bin + ...
          crosscorrspike(tclust1,tclust2,options.corr_tmax,2*nbins);
      end
      rateafter = zeros(1,n_log_bins);
      ratebefore = zeros(1,n_log_bins);
      rate1 = zeros(1,n_log_bins);
      for k = 1:n_log_bins
        rate1(k) = 1/diff(logbinBreak([k k+1]));
        rateafter(k) = sum(n_per_bin(logbinIndex{k}+nbins)).*rate1(k);
        ratebefore(k) = sum(n_per_bin(nbins-logbinIndex{k}+1)).*rate1(k);
      end
      ratebefore = ratebefore(end:-1:1); % reverse, since the bin indices got reversed
      hax1 = subplot(n_rows,n_cols,(i-1)*nclusters + j);
      hax2 = subplot(n_rows,n_cols,(j-1)*nclusters + i);
      col1 = unique_color(selectedIndex(i)-1,max_clusters);
      col2 = unique_color(selectedIndex(j)-1,max_clusters);
      setappdata(hax1,'clusterColor',col1);
      setappdata(hax2,'clusterColor',col2);
      [x1,y1] = stairs(logbinBreak,[ 0 rateafter]);
      [x2,y2] = stairs(logbinBreak,[ 0 ratebefore(end:-1:1)]);
      %%%y2 = [0; y2(length(y2):-1:2)];
      hpatch = patch(x1([1:end end]),[y1; 0],col1,...
        'EdgeColor',col1,'Parent',hax1);
      hpatch = patch(x2([1:end end]),[y2; 0],col2,...
        'EdgeColor',col2,'Parent',hax2);
      set(hax1,'XScale','log','TickDir','out','XLim',x1([1 end]))
      set(hax2,'XScale','log','TickDir','out','XLim',x2([1 end]))
      %%%set(hax2,'XDir','reverse')
      % Sometimes it makes bad decisions about tick marks along the X axis.
      % Check for this
      xtick = get(gca,'XTick');
      if (length(xtick) < 2)
        ulx10 = unique(round(log10(x(find(x)))));
        set([hax1 hax2],'XTick',10.^ulx10);
      end
      ylim1 = get(hax1,'YLim');
      ylim2 = get(hax2,'YLim');
      % Draw the 1/bin scale marks
%       hline = line(sqrt(logbinBreak(1:end-1).*logbinBreak(2:end)),...
%         rate1,...
%         'Color','k',...
%         'LineStyle','none',...
%         'Marker','.',...
%         'Parent',hax1);
%       copyobj(hline,hax2);
      setappdata(hax1,'cluster',i);
      setappdata(hax2,'cluster',j);
      setappdata(hax1,'dataSource',hfig_cc);
      setappdata(hax2,'dataSource',hfig_cc);
      % Draw the selector line
      %hline = line(x([4 4]),ylim,'Color','k','LineStyle','--');
      %drag_line(hline,struct('type','v','onDragDone',@cc_showwaveforms));
      if options.set_correlation_yaxes_without_fastest
        halfway_pt = ceil(length(y)/2);
        long_half1 = y1(halfway_pt:end);
        long_half2 = y2(1:halfway_pt);
        ylim1(2) = 1.1*max(long_half1);
        ylim2(2) = 1.1*max(long_half2);
        set(hax1,'YLim',ylim1)
        set(hax2,'YLim',ylim2)
      end
    end
  end
  %hfig_cc = figure;

  
function cc_showwaveforms(sender,event_args)
  % Show the reconstructed waveforms for snippets with sniptimes below
  % the user selection
  % First, get relevant data
  hline = sender;
  hax = get(sender,'Parent');
  hsource = getappdata(hax,'dataSource');
  clustIndex = getappdata(hsource,'clustIndex');
  shc = getappdata(hsource,'shc');
  col = getappdata(hax,'clusterColor');
  [nfiles,nclusters] = size(clustIndex);
  cluster = getappdata(hax,'cluster');
  xdata = get(hline,'XData');
  thresh = xdata(1);
  isac = (length(cluster) == 1);
  figure
  tot_pairs = 0;
  for k = 1:nfiles
    if isac
      % Compute autocorrelation to get the relevant spikes
      tIndex = clustIndex{k,cluster};
      t = double(shc(k).sniptimes(tIndex))/shc(k).scanrate;
      [tac,pairIndex] = autocorrspike(t, thresh);
      pairIndex = tIndex(pairIndex);
      % To read in the snippets, we need to put their indices in sorted
      % order.  However, we want to keep track of the permutations
      % required, so that we can get them back correctly
      [uIndex,tmp,uMapping] = unique(pairIndex(:));
      tsnips = sortheader_readsnips(shc(k),uIndex');
    else
      % Basically do the same thing for cross correlations
    end
    % Do reconstruction of pairs
    snipx = (shc(k).sniprange(1):shc(k).sniprange(2))/shc(k).scanrate;
    %col1 = unique_color(cluster(1)+coloffset,nclusters+coloffset);
    col1 = col(1,:);
    if isac
      col2 = col1;
    else
      col2 = col(2,:);
    end
    for i = 1:length(tac)
      line(snipx,tsnips(:,uMapping(2*i-1)),'Color',col1);
      line(snipx+tac(i),tsnips(:,uMapping(2*i)),'Color',col2);
    end
    tot_pairs = tot_pairs+length(tac);
  end
  title([num2str(tot_pairs) ' pairs'])
      