function g = register_expandg(gin,sz,method)
% REGISTER_EXPANDG: interpolate g to a larger image
% Syntax:
%   g = register_expandg(gin,sz)
%   g = register_expandg(gin,sz,method)
% where
%   gin is the deformation of space at one resolution scale (see
%     REGISTER_G0 for a detailed explanation);
%   sz is the image size at the next resolution scale;
%   method (default: 'spline') specifies the interpolation method;
%     use 'iminterp' if you want to use fast linear interpolation that is
%     light on memory requirements.
% and
%   g is the interpolated deformation for the larger scale.
%
% See also: REGISTER_G0.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 3)
    method = 'spline';
  end
  n_dims_g = length(gin);
  n_dims_each = ndims(gin{1});
  sz_pad = sz;                 % fill out dimensions ...
  sz_pad(end+1:n_dims_g) = 1;  % ...with 1s if needed
  if (all(size(gin{1}) == 1))
    % g is constant over space, so just copy it over all spatial
    % locations (but scaling it up to the new coordinates)
    g = register_g0(sz_pad);
    for dimIndex = 1:n_dims_g
      m = sz_pad(dimIndex);
      g{dimIndex} = g{dimIndex} + m*(gin{dimIndex}-1);
    end
    return
  end

  % g varies over space, so we have to do interpolation. But some
  % dimensions may be singletons, which makes interpolation go bad. We
  % deal with these by converting these coordinates to their final size.
  sz_g = size(gin{1});
  sz_g_pad = sz_g;
  sz_g_pad(end+1:n_dims_g) = 1;
  is_singleton = (sz_g_pad == 1);
  sz_g_filled = sz_g_pad;
  sz_g_filled(is_singleton) = sz_pad(is_singleton);
  rep_factor = sz_g_filled ./ sz_g_pad;  % replication factor
  if any(sz_g_pad == 1)
    for dimIndex = 1:n_dims_g
      if (sz_g_pad(dimIndex) > 1)
        gin{dimIndex} = repmat(gin{dimIndex},rep_factor);
      else
        rep_dim = sz_g_filled;
        rep_dim(dimIndex) = 1;
        m = rep_factor(dimIndex);
        gin{dimIndex} = repmat(m*(gin{dimIndex}-1),rep_factor) + ...
          repmat(reshape([1:rep_factor(dimIndex)],sz_g_filled ./ rep_dim),rep_dim);
      end
    end
  end
  x = cell(1,n_dims_g);
  y = cell(1,n_dims_g);
  g = cell(1,n_dims_g);
  % Set up the coordinates
  decfactors = round(sz_pad./sz_g_filled);
  for dimIndex = 1:n_dims_g
    %x{dimIndex} = 1+round((decfactors(dimIndex)-1)/2):decfactors(dimIndex):sz_pad(dimIndex);
    %y{dimIndex} = 1:sz_pad(dimIndex);
    if strcmp(method,'iminterp')
        y{dimIndex} = single(linspace(1,sz_g_filled(dimIndex),sz_pad(dimIndex)));
    else
        m = (sz_pad(dimIndex)-1)/(sz_g_filled(dimIndex)-1);
        x{dimIndex} = (0:sz_g_filled(dimIndex)-1)*m;
        y{dimIndex} = 0:sz_pad(dimIndex)-1;
    end
  end
  [y{:}] = ndgrid(y{:});
  % Do the interpolation; scale each coordinate appropriately
  for dimIndex = 1:n_dims_g
%     m = decfactors(dimIndex);
%     b = -round((decfactors(dimIndex)-1)/2);
%     g_tmp = m*gin{dimIndex}+b;   % Scaled g
%     g{dimIndex} = interpn(x{:},g_tmp,y{:},method) +
%     floor((x{dimIndex}(1)-1 - sz_pad(dimIndex)+x{dimIndex}(end))/2);
    m = (sz_pad(dimIndex)-1)/(sz_g_filled(dimIndex)-1);
    b = 0;
    g_tmp = m*(gin{dimIndex}-1)+b;
    if strcmp(method,'iminterp')
        g{dimIndex} = iminterp(g_tmp,y{:})+1;
    else
        g{dimIndex} = interpn(x{:},g_tmp,y{:},method)+1;
    end
  end
  