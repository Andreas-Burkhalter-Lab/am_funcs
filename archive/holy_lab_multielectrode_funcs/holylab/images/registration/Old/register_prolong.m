function uout = register_prolong(uin,sz_out)
  if ~iscell(uin)
    uin = {uin};
  end
  sz = size(uin{1});
  dimflag = sz > 1;
  sz = sz(dimflag);
  sz_out_keep = sz_out(dimflag);
  n_dims = length(sz);
  x = cell(1,n_dims);
  for dimIndex = 1:n_dims
    x{dimIndex} = single(0.5:0.5:sz_out_keep(dimIndex)/2);
    % Fix up right boundary to prevent NaNs
    x{dimIndex}(end) = x{dimIndex}(end) - 2*eps(x{dimIndex}(end));
  end
  X = cell(1,n_dims);
  if (n_dims > 1)
    [X{:}] = ndgrid(x{:});
  else
    X = x;
  end
  %uout = 2*iminterp(u,X{:},'extrap');
  for uIndex = 1:length(uin)
    % This version insures we do nearest-neighbor extrapolation, which
    % seems likely to cause less trouble than linear extrapolation. The
    % extrapolation only along the edge works in this case, or should,
    % because with prolongation we are less than 1 pixel beyond the edge.
    [uout{uIndex},w] = iminterp(uin{uIndex},X{:});
  end
