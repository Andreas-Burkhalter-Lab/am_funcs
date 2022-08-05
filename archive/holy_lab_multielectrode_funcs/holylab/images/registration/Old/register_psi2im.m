function im = register_psi2im(psi,options)
  if options.covariant
    im = psi.^2;
  else
    im = psi;
  end
  