function psi = register_im2psi(im,options)
  if options.covariant
    psi = sqrt(im);
  else
    psi = im;
  end
  