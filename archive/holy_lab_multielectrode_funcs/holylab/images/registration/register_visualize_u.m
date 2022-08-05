function register_visualize_u(varargin)
% REGISTER_VISUALIZE_U: view a deformation in 2d or 3d
% Syntax:
%   register_visualize_u(u)
%   register_visualize_u(u,rmg_params)
%   register_visualize_u(u,rmg_params,skip)
%   register_visualize_u(hax,...)
%
% where
%   u is the deformation
%   rmg_params is the output of register_multigrid_options; you need to
%     supply this if you want to scale the axes in units of physical space
%     (it uses the "pixel_spacing" fields of the image_grid)
%   skip (default 1) sets the interval between plotted points (an integer,
%     make larger than 1 to plot fewer points)
%   hax is an axis handle, if you want to direct the plot to a particular
%     axis.
%
% In the generated plot, the blue dot marks the destination to which the
% other end of each line is moved.

% Copyright 2010 by Timothy E. Holy

  % Parse inputs
  curarg = 1;
  if isscalar(varargin{1}) && ishandle(varargin{1})
    hax = varargin{1};
    curarg = 2;
  else
    hax = gca;
  end
  u = varargin{curarg};
  if (length(varargin) == curarg)
    % There was no rmg_params supplied, just do a grid plot with default
    % scaling
    if isempty(u)
      error('Can''t plot empty u without rmg_params');
    end
    pixel_spacing = ones(1,ndims(u)-1);
  else
    rmg_params = varargin{curarg+1};
    if ~isempty(u)
      szu = size(u);
      clevel = find(cellfun(@(x)isequal(x,szu(1:end-1)),{rmg_params.image_grid.sz}));
      pixel_spacing = rmg_params.image_grid(clevel).pixel_spacing;
    else
      pixel_spacing = rmg_params.image_grid(1+rmg_params.g_level_gap).pixel_spacing;
    end
  end
  skip = 1;
  if (length(varargin) > curarg+1)
    skip = varargin{curarg+2};
  end
  
  % Extract size information from u
  if isempty(u)
    % User supplied empty u, which corresponds to no deformation. We have
    % to get size information from rmg_params.
    sz = rmg_params.image_grid(1+rmg_params.g_level_gap).sz;
    n_dims = length(sz);
    u = zeros([sz n_dims]);
  else
    sz = size(u);
    n_dims = sz(end);
    sz = sz(1:n_dims);
  end
  
  % Generate two grids: undeformed and deformed
  g0 = register_g0(sz);
  g0 = cat(n_dims+1,g0{:});
  g = g0+u;
  if (skip > 1)
    x = cell(1,n_dims+1);
    for dimIndex = 1:n_dims;
      x{dimIndex} = 1:skip:sz(dimIndex);
    end
    x{n_dims+1} = 1:n_dims;
    g0 = g0(x{:});
    g = g(x{:});
    sz = size(g);
    sz = sz(1:n_dims);
  end
  
  % Pack in a form that plot can accept
  g0 = reshape(g0,[prod(sz) n_dims]);
  g0 = g0 .* repmat(pixel_spacing,[size(g0,1) 1]);
  g = reshape(g,[prod(sz) n_dims]);
  g = g .* repmat(pixel_spacing,[size(g,1) 1]);
  lines = cell(1,n_dims);
  for dimIndex = 1:n_dims
    lines{dimIndex} = [g(:,dimIndex) g0(:,dimIndex)]';
  end
  
  % Do the plot
  if (n_dims == 2)
    plot(hax,lines{:},'k');
    line(g0(:,1),g0(:,2),'parent',hax,'color','b','linestyle','none','marker','.')
  elseif (n_dims == 3)
    plot3(hax,lines{:},'k');
    line(g0(:,1),g0(:,2),g0(:,3),'parent',hax,'color','b','linestyle','none','marker','.')
    zlabel('coord 3')
  end
  xlabel('coord 1')
  ylabel('coord 2')
  axis equal