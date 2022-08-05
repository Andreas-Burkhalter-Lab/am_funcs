function [himgo,hlineso,htexto,sortIndex] = imagesc_clusters(X,clust,options)
% IMAGESC_CLUSTERS: plot data as image, sorted by cluster identity
% Syntax:
%   imagesc_clusters(X,clust,options)
%   [himg,hlines,htext] = imagesc_clusters(...) 
%   [himg,hlines,htext,sortIndex] = imagesc_clusters(...) 
% where
%   X is the data you want to plot, a d-by-N matrix, where data points are
%     the columns of X;
%   clust is a 1-by-N matrix giving the cluster # associated with each
%     column of X;
%   options is a structure which may have the following fields:
%     showtext (default true): if true, the cluster number will be
%       printed above each cluster of points
%     showlines (default true): if true, dashed lines will be used to
%       separate the columns in a cluster
%     clim: if supplied, overrides the default clims defined by the data;
%     permute (default false): if true, randomizes the order within each
%       cluster
%     symmetrical (default false): if true, picks a color scheme that is
%       symmetric around 0;
%     + any other options for colormap_pm (if symmetrical is true)
% and himg, hlines, and htext are handles to the image, lines, and text,
% respectively. sortIndex gives the ordering of columns in the output with
% respect to the input; this is particularly useful if permute=true.
%
% See also: colormap_pm.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'showtext',true,'showlines',true,'permute',false,'symmetrical',false);
  clim = [min(X(:)) max(X(:))];
  if options.symmetrical
    clim = [-1 1]*max(abs(clim));
  end
  if isfield(options,'clim')
    clim = options.clim;
  end
  d = size(X,1);
  
  hlines = [];
  htext = [];
  
  %% Plot the matrix
  hax = gca;
  [clabel,nlabel] = agglabel(clust);
  if options.permute
    for i = 1:length(clabel)
      rp = randperm(nlabel(i));
      clabel{i} = clabel{i}(rp);
    end
  end
  sortIndex = cat(2,clabel{:});
  himg = imagesc(X(:,sortIndex),clim);
  if (options.symmetrical)
    cm = colormap_pm(X,options);
    colormap(cm)
  end
  set(hax,'TickDir','out','Box','off')
  
  %% Show the cluster separation lines
  nlcum = cumsum(nlabel);
  n_clust = length(nlabel);
  if options.showlines
    col = 'w';
    if isfield(options,'background')
      if (ischar(options.background) && isequal(options.background,'w')) || ...
          isequal(options.background,[1 1 1])
        col = 'k';
      end
    end
    hlines = zeros(1,n_clust-1);
    for i = 1:n_clust-1
      hlines(i) = line([1 1]*nlcum(i)+0.5,[0.5 d+0.5],'Color',col,'LineStyle','--');
    end
  end
  ncenter = ([0 nlcum(1:end-1)]+1 + nlcum)/2;
  if options.showtext
    htext = zeros(1,n_clust);
    for i = 1:n_clust
      htext(i) = text(ncenter(i),0,num2str(i),'Parent',hax,'HorizontalAlignment','center');
    end
  end
  setappdata(hax,'nlcum',nlcum); % store cluster boundary data for any GUI applications

  %% Output only if requested
  % This prevents command-line output if user doesn't provide semicolon
  if (nargout > 0)
    himgo = himg;
    hlineso = hlines;
    htexto = htext;
  end