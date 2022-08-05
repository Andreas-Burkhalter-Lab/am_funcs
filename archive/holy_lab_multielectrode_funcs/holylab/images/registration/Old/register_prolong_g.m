function gout = register_prolong_g(gin,sz_out)
  sz = size(gin{1});
  dimflag = sz > 1;
  sz = sz(dimflag);
  sz_out_keep = sz_out(dimflag);
  n_dims = length(sz);
  % Perhaps this should just call register_prolong; the downside of this
  % is that one would need to allocate storage for the entire size of g,
  % not just for one component of g. So for now we are keeping this function.
  x = cell(1,n_dims);
  for dimIndex = 1:n_dims
    x{dimIndex} = single(0.5:0.5:sz_out_keep(dimIndex)/2);
    % Now fix up right boundary to prevent NaNs
    x{dimIndex}(end) = x{dimIndex}(end) - 2*eps(x{dimIndex}(end));
  end
  X = cell(1,n_dims);
  if (n_dims > 1)
    [X{:}] = ndgrid(x{:});
  else
    X = x;
  end
  %uout = 2*iminterp(u,X{:},'extrap');
  g0in = register_g0(sz);
  g0out = register_g0(sz_out_keep);
  for dimIndex = 1:n_dims
    % This version insures we do nearest-neighbor extrapolation, which
    % seems likely to cause less trouble than linear extrapolation. It
    % works in this case, or should, because we are less than 1 pixel
    % beyond the edge.
    [u,w] = iminterp(gin{dimIndex}-g0in{dimIndex},X{:});
    gout{dimIndex} = 2*u + g0out{dimIndex};
  end
  return
  
  % Fix boundaries
  n_dims = ndims(uout);
  sz_out = size(uout);
  colons = repmat({':'},1,n_dims);
  for dimIndex = 1:n_dims
    if (sz_out(dimIndex) > 1)
      for i = 1:3
        colons_tmp{i} = colons;
        colons_tmp{i}{dimIndex} = i;
      end
      uout(colons_tmp{1}{:}) = 2*uout(colons_tmp{2}{:}) - uout(colons_tmp{3}{:});
      if (sz_out(dimIndex) > 2*size(u,dimIndex))
        for i = 1:3
          colons_tmp{i}{dimIndex} = sz_out(dimIndex) - i;
        end
        uout(colons_tmp{1}{:}) = 2*uout(colons_tmp{2}{:}) - uout(colons_tmp{3}{:});
      end
      for i = 1:3
        colons_tmp{i}{dimIndex} = sz_out(dimIndex) - i+1;
      end
      uout(colons_tmp{1}{:}) = 2*uout(colons_tmp{2}{:}) - uout(colons_tmp{3}{:});
    end
  end
  % Now scale to next size
  uout = 2*uout;
  