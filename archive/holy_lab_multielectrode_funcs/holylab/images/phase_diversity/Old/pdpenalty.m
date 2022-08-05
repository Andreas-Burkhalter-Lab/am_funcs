function [val,grad,object] = pdpenalty(phi,im,H,pdk)
  sz = size(im);
  imsz = sz(1:end-1);
  K = sz(end);  % the # of diversity images
  
  % Compute the pupil functions and PSFs for the different diversity images
  Hk = zeros([imsz K],'single');  % Pupil function
  hk = zeros([imsz K],'single');  % PSF for coherent illumination
  for indx = 1:K
    Hk(:,:,indx) = H .* exp(i * (phi + pdk(:,:,indx)));
    hk(:,:,indx) = ifft2(Hk(:,:,indx));
  end
  sk = hk .* conj(hk);  % PSF for incoherent illumination
  Sk = zeros([imsz K],'single');
  for indx = 1:K
    Sk(:,:,indx) = fft2(sk(:,:,indx));
  end
  S2tot = sum(Sk .* conj(Sk),3);
  ukeep = S2tot > eps;
  ukeepv = ukeep(:);

  % Compute the transforms of the acquired images
  Dk = zeros(sz,'single');
  for indx = 1:K
    Dk(:,:,indx) = fft2(im(:,:,indx));
  end
  D2tot = sum(Dk .* conj(Dk),3);
  DdotS = sum(Dk .* conj(Sk),3);
  
  % Compute the value of the penalty
  num = DdotS .* conj(DdotS);
  val = -sum(num(ukeepv) ./ S2tot(ukeepv)) + sum(D2tot(:));
  
  % Compute the gradient of the penalty
  if (nargout > 1)
    Zk = zeros([K imsz],'single');
    coef1 = S2tot .* DdotS;
    coef2 = DdotS .* conj(DdotS);
    for indx = 1:K
      num = coef1 .* conj(Dk(:,:,indx)) - coef2 .* conj(Sk(:,:,indx));
      Zk(indx,ukeepv) = num(ukeepv) ./ S2tot(ukeepv).^2;
    end
    Zk = permute(Zk,[2 3 1]);
    ZconvH = zeros([imsz K],'single');
    for indx = 1:K
      ZconvH(:,:,indx) = fft2(ifft2(Zk(:,:,indx)) .* ...
			      ifft2(conj(Hk(:,:,indx))));
    end
    grad = 4*imag(sum(Hk .* ZconvH,3));
  end
  
  % Return the value of the image object
  if (nargout > 2)
    F = zeros(size(DdotS),'single');
    F(ukeepv) = DdotS(ukeepv)./S2tot(ukeepv);
    object = ifft2(F);
  end
end
