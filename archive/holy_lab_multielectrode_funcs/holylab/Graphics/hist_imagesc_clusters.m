function options = hist_imagesc_clusters(clustc,Xc,options)
% hist_imagesc_clusters: display abundance of clusters across datasets
% This function plots the fraction of observations that fall into each
% cluster for multiple datasets. It also plots the cluster mean as an
% image.
%
%
% Syntax:
%   ret = hist_imagesc_clusters(clustc,Xc,options)
% where
%   clustc is a 1-by-n_groups cell array, each element being a vector
%     containing the cluster assignment of the observations in the
%     corresponding data set
%   Xc is a 1-by-n_groups cell array, each element being a matrix that
%     contains all the observations in a given data set (each should have
%     the same number of rows, but will have different columns).
%   options may have the following fields:
%     scaling (default 'relative'): 'absolute' plots the total number of
%       observations by cluster/dataset, 'relative' plots the fraction of
%       total observations in each dataset across clusters. 'relative' is
%       particularly useful if your datasets have different numbers of
%       points in them but this fact is not in itself meaningful.
%     histax and imageax: use these two fields to supply particular axes
%       for the histogram and image plots, respectively (default is to
%       create a new figure).
%     repeats (default empty cell): supply a cell specifying which datasets are
%       repeats and will be merged together, error bar will be added in
%       this case. e.g.{[1 2 3],[4 5 6],[7 8 9]} will merge dataset Xc{1}, Xc{2},
%       Xc{3} together as 3 repeats.The same as dataset 4,5,6 and 7,8,9.
%
% On output, ret is a structure holding the options settings and handles
% for various components of the plot.
%
% See also: imagesc_clust, imagesc_clust_multiaxis.

% Copyright 2011 by Timothy E. Holy, modified by Pei S. Xu (Jan 2012)
  
  %% Parse arguments
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'scaling','relative');
  if ~isfield(options,'histax') && ~isfield(options,'imageax')
    figure
    options.histax = subplot(2,1,1);
    options.imageax = subplot(2,1,2);
  end
  
  %% Compute the average of each cluster & numbers of examples
  n_groups = length(clustc);
  n_clust = max(cellfun(@max,clustc));
  d = size(Xc{1},1);
  n = zeros(n_groups,n_clust);
  mu = zeros(d,n_clust);
  for groupIndex = 1:n_groups
    [clabel,ntmp] = agglabel(clustc{groupIndex});
    n(groupIndex,1:length(ntmp)) = ntmp;
    for clustIndex = 1:length(clabel)
      mu(:,clustIndex) = mu(:,clustIndex) + sum(Xc{groupIndex}(:,clabel{clustIndex}),2);
    end
  end
  ntot = sum(n,1);
  mu = bsxfun(@rdivide,mu,ntot);
  
  switch options.scaling
    case 'absolute'
      % Do nothing
    case 'relative'
      n = bsxfun(@rdivide,n,sum(n,2));
    otherwise
      error('options.scaling is not recognized');
  end
  
  %% Merge repeats if there are
  options = default(options,'repeats',{});
  if ~isempty(options.repeats)
    n_groups = length(options.repeats);
    n_rep_mean = zeros(n_groups, n_clust);
    n_rep_std = zeros(n_groups, n_clust);
    for groupIdx = 1:n_groups
      n_rep_mean(groupIdx,:) = mean(n(options.repeats{groupIdx},:),1);
      n_rep_std(groupIdx,:) = std(n(options.repeats{groupIdx},:),1);
    end
  end
  %% Plot the histogram
  if isempty(options.repeats)
    options.hbar = bar(options.histax,1:n_clust,n');
  else
    axes(options.histax);
    options.hbar = barweb(n_rep_mean',n_rep_std', [], [], [], [], [], [], [], [], 1, 'axis');
  end
  %% Plot the means
  options.himg = image(mu,'Parent',options.imageax,'CDataMapping','scaled');
  set(options.histax,'XLim',get(options.imageax,'XLim'))
  set([options.histax options.imageax],'TickDir','out','Box','off');
end  
  
    