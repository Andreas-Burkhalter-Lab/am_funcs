function [Hk,sk,im] = pd_forward_model_2d(phi,H,object)
% PD_FORWARD_MODEL_2D: compute PSFs given phase aberrations
%
% Given a pupil function and a phase aberration, compute the PSF for
% incoherent illumination.  To facilitate use in phase-diversity
% calculations, one can supply a set of K phase aberrations and
% (optionally) K different pupil functions.
%
% If you also supply an underlying object, this function can also
% calculate the resulting images.
%
% Here's the optical model (notation from Paxman et al, 1992 JOSA): each
% image is produced by convolution of the underlying object with a
% point-spread function s_k, which can be
% written as
%      s_k = |h_k|^2,
% where h_k is the optical transfer function for coherent illumination.
% The fourier transform H_k = fft2(h_k) satisfies
%    H_k = |H| exp(i * phi_k).
%
% Syntax:
%   [Hk,sk] = pd_forward_model_2d(phi,H)
%   [Hk,sk,im] = pd_forward_model_2d(phi,H,object)
% where
%   phi is a m-by-n-by-K array containing the phase aberrations; 
%   H is an m-by-n array containing the magnitude of the pupil
%     function (often this will be 1 inside the pupil and 0
%     outside). Optionally, this can be m-by-n-by-K if you need different
%     pupils for different images.
%   object (optional) is an m-by-n matrix containing the
%     pixel-intensities of the underlying object (_not_ the
%     diffraction-limited image)
% and
%   Hk is m-by-n-by-K, containing the Fourier transforms of the
%     coherent-illumination PSFs;
%   sk is m-by-n-by-K, containing the real-space incoherent-illumination
%     PSFs;
%   im is m-by-n-by-K, containing the output images with incoherent
%     illumination (available only if object is supplied).
%
% See also: PDPENALTY.

% Copyright 2008-2009 by Timothy E. Holy

  sz = size(phi);
  if (ndims(phi) > 2)
    imsz = sz(1:end-1);
    K = sz(end);  % the # of diversity images
  else
    K = 1;
    imsz = sz;
  end
  Hsz = size(H);
  multiple_pupils = false;
  if (length(Hsz) == length(imsz) && Hsz(end) == K)
    multiple_pupils = true;
  end
  
  % Compute the pupil functions for the different diversity images
  Hk = zeros([imsz K],class(phi));  % Pupil function
  hk = zeros([imsz K],class(phi));  % PSF for coherent illumination
  for indx = 1:K
    if multiple_pupils
      Htmp = H(:,:,indx);
    else
      Htmp = H;
    end
    Hk(:,:,indx) = Htmp .* exp(i * (phi(:,:,indx)));
    hk(:,:,indx) = ifft2(Hk(:,:,indx));
  end
  % Compute the transfer functions for incoherent illumination
  sk = hk .* conj(hk);  % PSF for incoherent illumination
  
  if (nargout > 2)
    % Calculate the images
    im = zeros([imsz K],class(object));
    objfft = fft2(object);
    for indx = 1:K
      im(:,:,indx) = ifft2(fft2(sk(:,:,indx)).*objfft);
    end
  end
end
