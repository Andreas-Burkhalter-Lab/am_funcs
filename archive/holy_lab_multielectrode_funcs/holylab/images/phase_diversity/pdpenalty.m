function [val,grad,object] = pdpenalty(phi,im,H)
% PDPENALTY: value and gradient of least square error for phase diversity
%
% Adapted fron Paxman et al, JOSA A 1992.
% Given K "views" (images), this function computes the error in a model
% in which each view arises as an image of an underlying two-dimensional
% object, assuming incoherent illumination. Each image is allowed to have
% a different phase aberration in the pupil plane.
% To facilitate optimization, this function also computes the gradient of
% the error with respect to the parameters (the phase aberrations).
%
% The model: each image is produced by convolution of the underlying
% object with a point-spread function s_k, which can be written as
%      s_k = |h_k|^2,
% where h_k is the optical transfer function for coherent illumination.
% The fourier transform H_k = fft2(h_k) satisfies
%    H_k = |H| exp(i * phi_k).
%
% Syntax:
%   [val,grad] = pdpenalty(phi,im,H)
% where
%   phi is a m-by-n-by-K array containing the current guess for the phase
%     aberrations; 
%   im is an m-by-n-by-K array containing the K images;
%   H is an m-by-n array containing the magnitude of the pupil
%     function (often this will be 1 inside the pupil and 0
%     outside). Optionally, this can be m-by-n-by-K if you need different
%     pupils for different images.
% and
%   val is the square error
%   grad is the gradient of the square error with respect to phi, and
%     hence is m-by-n-by-k
%
% Using the chain rule, you can place various restrictions or parametric
% forms on phi_k, so this function is a general foundation for such
% optimization problems.
%
% Optionally, the syntax
%   [val,grad,object] = pdpenalty(phi,im,H)
% also returns the "object estimate," the underlying object that produces
% the least square error.
% 
% See also: PD_DEMO, PDOPT_ZERNIKE, PD_FORWARD_MODEL_2D.

% Copyright Diwakar Turaga & Timothy E. Holy, 2008-2009

  sz = size(im);
  sz(end+1:3) = 1;  % make sure it has at least 3 entries
  imsz = sz(1:end-1);
  K = sz(end);  % the # of diversity images
  
  [Hk,sk] = pd_forward_model_2d(phi,H);

  % Compute the Fourier transforms of the PSFs
  Sk = zeros([imsz K],class(phi));
  for indx = 1:K
    Sk(:,:,indx) = fft2(sk(:,:,indx));
  end
  S2tot = sum(Sk .* conj(Sk),3);
  ukeep = S2tot > 10*eps;  % only retain pixels in pupil plane that have power
  ukeepv = ukeep(:);

  % Compute the Fourier transforms of the acquired images
  Dk = zeros(sz,class(im));
  for indx = 1:K
    Dk(:,:,indx) = fft2(im(:,:,indx));
  end
  D2tot = sum(Dk .* conj(Dk),3);
  DdotS = sum(Dk .* conj(Sk),3);
  
  % Compute the value of the penalty
  num = DdotS .* conj(DdotS);
  val = -sum(num(ukeepv) ./ S2tot(ukeepv)) + sum(D2tot(:));
  val = val / prod(imsz);
  
  % Compute the gradient of the penalty
  if (nargout > 1)
    Zk = zeros([K imsz],class(phi));
    coef1 = S2tot .* DdotS;
    coef2 = DdotS .* conj(DdotS);
    for indx = 1:K
      num = coef1 .* conj(Dk(:,:,indx)) - coef2 .* conj(Sk(:,:,indx));
      Zk(indx,ukeepv) = num(ukeepv) ./ S2tot(ukeepv).^2;
    end
    Zk = permute(Zk,[2 3 1]);
    ZconvH = zeros([imsz K],class(phi));
    for indx = 1:K
      ZconvH(:,:,indx) = fft2(ifft2(Zk(:,:,indx)) .* ...
			      ifft2(conj(Hk(:,:,indx))));
    end
    grad = 4*imag(Hk .* ZconvH) / prod(imsz);
  end
  
  % Return the value of the image object
  if (nargout > 2)
    F = zeros(size(DdotS),class(im));
    F(ukeepv) = DdotS(ukeepv)./S2tot(ukeepv);
    object = ifft2(F);
  end
end
