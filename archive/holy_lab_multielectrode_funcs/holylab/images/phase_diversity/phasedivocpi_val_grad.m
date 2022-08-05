function sout = phasedivocpi_val_grad(inputparams,outputparams)
% PHASEDIVOCPI_VAL_GRAD: core computations for homogenous OCPI phase diversity
% Note: you probably want to call ocpi_pdwrapper instead of this function.
%
% Syntax:
%  sout = phasedivocpi_val_grad(inputparams,outputparams)
% where the fields of inputparams are (see also OCPI_PDWRAPPER for better
% documentation):
% phik (the aberration functions) is framesize-by-K, and is wrapped so
%   that the center of the PSF  is in the upper left corner (see
%   filter_wrap)
% v (the z-profile of the illumination) is stacklen, and is wrapped so
%   that z=0 corresponds to the first element
% H0 (the pupil) is framesize, and is wrapped
% rho (the pupil radius) is framesize, and is wrapped
% c0 is a scalar, scaled to accept integer-valued zs
%
% f (the base object) is framesize-by-stacklen
%
% Iobsk is framesize-by-stacklen-by-K
% mask is framesize-by-stacklen
    
% Note: normalize_s is unfinished!

  phik = inputparams.phik;
  v = inputparams.v;
  H0 = inputparams.H0;
  rho = inputparams.rho;
  c0 = inputparams.c0;
  
  compute_Iobs = isfield(inputparams,'f');
  compute_err = compute_Iobs && isfield(inputparams,'Iobsk');
  inputparams = default(inputparams,'normalize_s',false);
  outputparams = default(outputparams,'err',compute_err,'PSFs',~compute_Iobs,'Iobs_calc',false,'fgrad',false,'vgrad',false);

  imclass = class(phik);
  sz = size(phik);
  if (length(sz) < 3)
    sz(3) = 1;
  end
  K = sz(3);  % the # of diversity images
  framesz = sz([1 2]);
  stacklen = length(v);   % the # of images in the stack
  imsz = [framesz stacklen];

  z = (0:stacklen-1);
  zwrap = z > stacklen/2;
  z(zwrap) = z(zwrap) - stacklen;

  cz = c0*z;
  
  % Compute the pupil functions and PSFs for the different diversity images
  % This is different from pd_forward_model_2d because of the z-term
  Hk = zeros([imsz K],imclass);  % Pupil function
  hk = zeros([imsz K],imclass);  % PSF for coherent illumination
  for Kindx = 1:K
    for stackindx = 1:stacklen
      tmp = H0 .* exp(i * (phik(:,:,Kindx) + cz(stackindx)*rho.^2));
      Hk(:,:,stackindx,Kindx) = tmp;
      hk(:,:,stackindx,Kindx) = ifft2(tmp);
    end
  end
  sk = hk .* conj(hk);  % PSF for incoherent illumination
  for stackindx = 1:stacklen
    % Include the illumination
    sk(:,:,stackindx,:) = sk(:,:,stackindx,:) * v(stackindx);
  end
  snorm = ones(1,K);
  if inputparams.normalize_s
    for Kindx = 1:K
      sK = sk(:,:,:,Kindx);
      snorm(Kindx) = sum(sK(:));
      sk(:,:,:,Kindx) = sK / snorm(Kindx);
    end
  end
  
  % Return the result?
  if outputparams.PSFs
    sout.PSFs = sk;
  end
  
  % Check for an early return (PSF only)
  if ~compute_Iobs
    return
  end
  
  % Convolve the PSFs with the base object.
  Sk = zeros([imsz K],imclass);
  for Kindx = 1:K
    Sk(:,:,:,Kindx) = fftn(sk(:,:,:,Kindx));
  end
  Iobs_calc = zeros([imsz K],imclass);
  I0fft = fftn(inputparams.f);
  for Kindx = 1:K
    Iobs_calc(:,:,:,Kindx) = ifftn(Sk(:,:,:,Kindx) .* I0fft);
  end
  
  % Return the result?
  if outputparams.Iobs_calc
    sout.Iobs_calc = Iobs_calc;
  end

  % Check for an early return (observed images only)
  if ~compute_err
    return
  end
  
  if isfield(inputparams,'mask')
    mask = inputparams.mask;
  else
    mask = ones(imsz,imclass);
  end
  
  % Compute the residual and the value of the functional
  dIk = zeros([imsz K],imclass);
  dIkfft = zeros([imsz K],imclass);
  val = 0;
  for Kindx = 1:K
    dIk_tmp = Iobs_calc(:,:,:,Kindx) - inputparams.Iobsk(:,:,:,Kindx);
    dIk_tmp = dIk_tmp .* mask;
    dIk(:,:,:,Kindx) = dIk_tmp;
    dIkfft_tmp = fftn(dIk_tmp);
    dIkfft(:,:,:,Kindx) = dIkfft_tmp;
    tmp = dIkfft_tmp .* conj(dIkfft_tmp);
    val = val + sum(tmp(:));
  end
  sout.val = val/numel(inputparams.f);
  
  % Compute the gradient with respect to the object
  if outputparams.fgrad
    fgrad = zeros(imsz,imclass);
    for Kindx = 1:K
      fgrad = fgrad + ifftn(conj(Sk(:,:,:,Kindx)) .* dIkfft(:,:,:,Kindx));
    end
    sout.fgrad = 2*fgrad;
  end
  
  if outputparams.vgrad || outputparams.phikgrad
    % Compute the flipped-correlation of the residual with the base object
    I0wrap = filter_wrap(inputparams.f,imsz);
    I0fftc = conj(fftn(I0wrap));
    IcI = zeros([imsz K],imclass);
    for Kindx = 1:K
      IcItmp = ifftn(dIkfft(:,:,:,Kindx) .* I0fftc);
      IcI(:,:,:,Kindx) = filter_wrap(IcItmp,imsz);
    end
  end

  % Compute the gradient with respect to the illumination PSF
  if outputparams.vgrad
    vgrad = zeros(imsz,imclass);
    for Kindx = 1:K
      vgrad = vgrad + hk(:,:,:,Kindx) .* conj(hk(:,:,:,Kindx)) .* ...
	      IcI(:,:,:,Kindx);
    end
    vgrad = sum(vgrad,1);
    vgrad = sum(vgrad,2);
    sout.vgrad = reshape(2*vgrad(:),size(v));
  end
    
  % Compute the gradient with respect to phik (here as a
  % general function of pupil coordinates; other functions might project
  % this onto a Zernike basis, for example)
  if outputparams.phikgrad
    sout.phikgrad = zeros([framesz K],imclass);
    for Kindx = 1:K
      tmp = fft2(hk(:,:,:,Kindx) .* IcI(:,:,:,Kindx));
      tmp = conj(tmp);
      for stackindx = 1:stacklen
        tmp(:,:,stackindx) = tmp(:,:,stackindx) * v(stackindx);
      end
      tmp = tmp .* Hk(:,:,:,Kindx);
      sout.phikgrad(:,:,Kindx) = -4 * imag(sum(tmp,3)) / numel(H0);
    end
  end
end
