function hline_out = clusterplot(hax,x,clust,options)
% CLUSTERPLOT: show clusters as colored dots
% Syntax:
%   clusterplot(x,clust)
%   hline = clusterplot(hax,x,clust)
%   hline = clusterplot(...,options)
% where
%   x is a d-by-N matrix containing the data points. The behavior depends
%     on the dimensionality d:
%       For d = 1, a histogram is plotted
%       For d = 2 or d > 3, the first two dimensions are plotted as a
%         scatter plot. Points in different clusters are shown as dots; any
%         clusters containing only a single point are shown with a circle.
%         You can click on the individual points, and the cluster # is
%         displayed on the command line.
%       For d = 3, a three-dimensional scatter plot is shown.
%   clust is 1-by-N, containing the cluster # assigned to each point
%   hax is an axis handle in which you want the plot to appear
%   options is a structure which may have the following fields:
%     color: an n-by-3 matrix specifying the RGB color of each cluster on
%       as rows. If the colors are exhausted, it cycles through them again.
%       By default, unique_color is used to 
%     multidraw (default false): if true, the screen is updated after
%       drawing each cluster; 
%       this can help show the boundaries in cases where there are many
%       clusters.
%

% Copyright 2007-2009 Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end
  if ~isscalar(hax) || ~ishandle(hax(1))
    if (nargin > 2)
      options = clust;
    end
    clust = x;
    x = hax;
    figure
    hax = gca;
  end
  %col = get(gca,'ColorOrder');
  %n_colors = size(col,1);
  %if (max(clust) > n_colors)
  %  warning('There are more clusters than colors available');
  %end
  
  lclust = agglabel(clust);
  n_clust = length(lclust);
  if ~isfield(options,'color')
    options.color = zeros(n_clust,3);
    for i = 1:n_clust
      col = unique_color(i,n_clust);
      if isequal(col,get(hax,'Color'))
        col = 1 - col;  % for black backgrounds
      end
      options.color(i,:) = col;
    end
  end
  options = default(options,'multidraw',false);
  
  d = size(x,1);
  if (d == 1)
    xc = linspace(min(x),max(x),100);
  else
    xr = [min(x(1,:)) max(x(1,:))];
    yr = [min(x(2,:)) max(x(2,:))];
    set(hax,'XLim',xr,'Ylim',yr)
  end
  if (nargout > 0)
    hline_out = nan(1,n_clust);
  end
  for i = n_clust:-1:1
    col = options.color(mod(i-1,size(options.color,1))+1,:);
    if (d > 1)
      if (d == 2 || d > 3)
        hline = line(x(1,lclust{i}),x(2,lclust{i}),'LineStyle','none',...
          'Color',col,'Marker','.','Parent',hax,'ButtonDownFcn',sprintf('disp(''cluster %d'')',i));
      else
        hline = line(x(1,lclust{i}),x(2,lclust{i}),x(3,lclust{i}),...
          'LineStyle','none','Color',col,'Marker','.','Parent',hax);
      end
      if (length(lclust{i}) == 1)
        set(hline,'Marker','o')
      end
      if (nargout > 0 && ~isempty(hline))
        hline_out(i) = hline;
      end
      if options.multidraw
        drawnow % plot them one at a time so one can tell the boundary extents
      end
    else
      n = hist(x(lclust{i}),xc);
      hb = bar(hax,xc,n,1);
      hold(hax,'on')
      set(hb,'FaceColor',col,'EdgeColor','none');
    end
  end
 


 
