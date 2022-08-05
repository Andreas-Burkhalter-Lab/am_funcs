function gOut = register_reduceg(gIn,decfactors)
  n_dims = length(gIn);
  sz_big = size(gIn{1});
  sz_big(end+1:n_dims) = 1;
  g0_big = register_g0(sz_big);
  sz_small = size(imreduce(gIn{1},decfactors));
  sz_small(end+1:n_dims) = 1;
  g0_small = register_g0(sz_small);
  for dimIndex = 1:n_dims
      dg = gIn{dimIndex} - g0_big{dimIndex};
      dg =dg/decfactors(dimIndex);
      dg = imreduce(dg,decfactors);
      gOut{dimIndex} = dg + g0_small{dimIndex};
  end
