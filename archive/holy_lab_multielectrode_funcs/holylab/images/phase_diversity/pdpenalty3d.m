function [val,grad,object] = pdpenalty3d(phi,im,H,pdk, Hz)

%calculates phase diversity penalty
% Adapted fron Paxman et al, JOSA A 1992.
%
% Adapted from pdpenalty.m
% This version takes in two aberrated images and calculates the
% value and gradient
%
% INPUT
% phi = current estimate of unknown aberration
% im = contains K phase diverse stacks
% H = pupil function
% pdk = K known phase aberrations
% Hz = 3d pupil function ... obtained by applying defocus to pupil function
%                            see Hanser et al, J. Microscopy, vol 216, pg
%                            32, year 2004.
%
% OUTPUT
% val = val of the function we are trying to minimize
% grad = gradient of the function
% object = idealized stack
% 
% Diwakar Turaga,
%
% See also: PDPENALTY, PDOPT_ZERNIKE3D, PD_DEMO3D

  sz = size(im);
  imsz = sz(1:end-1);
  K = sz(end);  % the # of diversity images
  
  % Compute the pupil functions and PSFs for the different diversity images
  Hk = zeros([imsz K],class(H));  % Pupil function
  hk = zeros([imsz K],class(H));  % PSF for coherent illumination
  for indx1 = 1:K
      
      for indx2 = 1:sz(3)
          Hk(:,:,indx2, indx1) = Hz(:,:,indx2) .* exp(i * (phi + pdk(:,:,indx1)));
          hk(:,:,indx2, indx1) = ifft2(Hk(:,:,indx2, indx1));
      end
  end
  
  sk = hk .* conj(hk);  % PSF for incoherent illumination
  Sk = zeros([imsz K],class(sk));
  for indx = 1:K
      Sk(:,:,:, indx) = fftn(sk(:,:,:, indx));
  end
  
  S2tot = sum(Sk .* conj(Sk),4);
  ukeep = S2tot > eps;
  ukeepv = ukeep(:);

  % Compute the transforms of the acquired images
  Dk = zeros(sz,class(sk));
  for indx = 1:K
        Dk(:,:,:, indx) = fftn(im(:,:,:, indx));
  end
  
  D2tot = sum(Dk .* conj(Dk),4);
  DdotS = sum(Dk .* conj(Sk),4);
  
  % Compute the value of the penalty
  num = DdotS .* conj(DdotS);
  val = -sum(num(ukeepv) ./ S2tot(ukeepv)) + sum(D2tot(:));
  
  % Compute the gradient of the penalty
  if (nargout > 1)
    Zk = zeros([K imsz],class(im));
    coef1 = S2tot .* DdotS;
    coef2 = DdotS .* conj(DdotS);
    for indx = 1:K
      num = coef1 .* conj(Dk(:,:,:,indx)) - coef2 .* conj(Sk(:,:,:,indx));
      Zk(indx,ukeepv) = num(ukeepv) ./ S2tot(ukeepv).^2;
    end
    Zk = permute(Zk,[2 3 4 1]);
    ZconvH = zeros([imsz K],class(im));
    
    for indx = 1:K
            
        ZconvH(:,:,:, indx) = fftn(ifftn(Zk(:,:,:, indx)) .* ...
                 ifftn(conj(Hk(:,:,:, indx))));

    end
    grad = 4*imag(sum(Hk .* ZconvH,4));
  end
  
  % Return the value of the image object
  if (nargout > 2)
    F = zeros(size(DdotS), class(im));
    F(ukeepv) = DdotS(ukeepv)./S2tot(ukeepv);

    object = ifftn(F);
    

  end
end
