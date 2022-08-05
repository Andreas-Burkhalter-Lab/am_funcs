function options = imagesc_clusters_multiaxis(clustc,Xc,options)
% imagesc_clusters_multiaxis: show clusters across datasets
% This function plots the observations that fall into each cluster for
% multiple datasets. It places each cluster in a separate axis to ensure
% alignment across datasets, and uses a uniform width so that number of
% observations can be quantified.
%
% Syntax:
%   ret = imagesc_clusters_multiaxis(clustc,Xc,options)
% where
%   clustc is a 1-by-n_groups cell array, each element being a vector
%     containing the cluster assignment of the observations in the
%     corresponding data set
%   Xc is a 1-by-n_groups cell array, each element being a matrix that
%     contains all the observations in a given data set (each should have
%     the same number of rows, but will have different columns).
%   options may have the following fields:
%     gap_horiz and gap_vert: specify the size of the gap between the axes,
%       in normalized units (i.e., 1.0 means full figure width).
%     centered (default true): if true, each axis will be centered; if
%       false, left-justified
%     clim (default min/max): the colorscaling limits for each axis
%
% On output, ret is a structure holding the options settings and handles
% for various components of the plot.
%
% See also: imagesc_clust, imagesc_clust_multiaxis.

% Copyright 2011 by Timothy E. Holy
    
  %% Parse arguments
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'hax')
    options.hax = gca;
  end
  
  n_groups = length(clustc);
  n_clust = max(cellfun(@max,clustc));
  d = size(Xc{1},1);
  options = default(options,'gap_horiz',0.2/n_clust,'gap_vert',0.5/n_groups,'centered',true);
  if ~isfield(options,'clim')
    Xmin = inf;
    Xmax = -inf;
    for i = 1:n_groups
      thisXmin = min(Xc{i}(:));
      thisXmax = max(Xc{i}(:));
      Xmin = min(Xmin,thisXmin);
      Xmax = max(Xmax,thisXmax);
    end
    options.clim = [Xmin Xmax];
  end

  %% Count # in each cluster by group
  n = zeros(n_groups,n_clust);
  clabel = cell(1,n_groups);
  for i = 1:n_groups
    [clabel{i},ntmp] = agglabel(clustc{i});
    n(i,1:length(ntmp)) = ntmp;
    if (length(ntmp) < n_clust)
      clabel{i}{n_clust} = [];
    end
  end
  
  %% Calculate the horizontal splits
  nmax = max(n,[],1);
  tot = sum(nmax);
  gap = round(options.gap_horiz*tot);
  npergroup = [nmax; repmat(gap,1,n_clust)];
  ncum = cumsum(npergroup(:));
  ncum = ncum(1:end-1);  % no gap needed after last cluster
  hsplit = ncum/ncum(end);
  hsplit = hsplit(1:end-1);  % no split needed at end of last axis
  
  %% Create the axes
  vsplit = SplitAxesEvenly(n_groups,options.gap_vert);
  if ~isempty(vsplit)
    options.hax = SplitGrid(hsplit,vsplit,[],[],options.hax);
  else
    options.hax = SplitHoriz(hsplit,repmat([1 0],1,n_clust),options.hax);
  end
  
  %% Do the plots
  for i = 1:n_groups
    for j = 1:n_clust
      xc = 1:length(clabel{i}{j});
      if options.centered
        xc = xc + round((nmax(j)-length(xc))/2);
      end
      options.himg(i,j) = image(xc,1:d,Xc{i}(:,clabel{i}{j}),'Parent',options.hax(i,j),'CDataMapping','scaled');
    end
  end
  for j = 1:n_clust
    set(options.hax(:,j),'XLim',[0.5 nmax(j)+0.5]);
  end
  set(options.hax,'CLim',options.clim,'Box','off','TickDir','out')
end
  
  
  