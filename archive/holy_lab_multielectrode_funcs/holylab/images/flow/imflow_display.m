function hlout = imflow_display(varargin)
% imflow_display: show a vector flow field on a regular pixel grid (2D uc
% only). Use imflow_display3D for displaying 3D uc.
%
% Syntax:
%   imflow_display(dX)
%   imflow_display(hax,dX)
%   hline = imflow_display(...)
% where
%   dX is a nrows-by-ncols-by-2 matrix, dX(i,j,:) is the flow vector at
%     pixel i,j;
%   hax (default gca) is the axis to use for display
% and
%   hline is a 2-vector, the first containing the handle for the "line of
%     lines" illustrating the flow, the second for the (default blue) dots
%     marking the end point of each flow vector.

% Copyright 2009-2010 by Timothy E. Holy

  if (isscalar(varargin{1}) && ishandle(varargin{1}))
    hax = varargin{1};
    dx = varargin{2};
  else
    hax = gca;
    dx = varargin{1};
  end
  sz = size(dx);
  sz_spatial = sz(1:end-1);
  n_dims = length(sz_spatial);
  x = cell(1,n_dims);
  for i = 1:n_dims
    x{i} = 1:sz_spatial(i);
  end
  X = cell(1,n_dims);
  [X{:}] = ndgrid(x{:});
  start = cat(ndims(dx),X{:});
  finish = start+dx;
  start = reshape(start,[prod(sz_spatial) n_dims]);
  finish = reshape(finish,[prod(sz_spatial) n_dims]);
  
  linevertices = [start finish nan(size(start))];
  xx = cell(1,n_dims);
  for i = 1:n_dims
    xx{i} = linevertices(:,i:n_dims:end)';
  end
  if (n_dims == 2)
    hl = line(xx{2}(:),xx{1}(:),'Parent',hax,'Color','k');
    hl(2) = line(finish(:,2),finish(:,1),'Parent',hax,'Color','b','LineStyle','none','Marker','.');
    set(hax,'YDir','reverse','DataAspectRatio',[1 1 1])
  elseif (n_dims == 3)
    hl = line(xx{2}(:),xx{1}(:),xx{3}(:),'Parent',hax,'Color','k');
    hl(2) = line(finish(:,2),finish(:,1),finish(:,3),'Parent',hax,'Color','b','LineStyle','none','Marker','.');
  end
  if (nargout > 0)
    hlout = hl;
  end