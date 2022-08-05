function [Zcoefs,phi] = pdopt_zernike(im,H,pdk,rho,theta,Zcoefs)
  Zval = zernike_values(rho,theta,H,length(Zcoefs));
  optfun = @(Zc) pdpenalty_zernike(Zc,im,H,pdk,Zval);
  %minops = optimset('GradObj','on','DerivativeCheck','on','DiffMinChange',1e-4,'Display','iter');
  minops = optimset('GradObj','on','DiffMinChange',1e-4,'Display','iter');
  Zcoefs = fminunc(optfun,Zcoefs,minops);
  if (nargout > 1)
    phi = Zcoefs2phi(Zcoefs,Zval);
  end
end

function Zval = zernike_values(rho,theta,H,n_coefs)
% Compute the Zernike polynomials
% The first 3 Zernike's are omitted, since there is nothing in the data
% that allows us to determine them (and they do not affect image
% quality, only image position)
  rhosz = size(rho);
  rho(rho > 1) = 0;
  Zval = zeros([rhosz n_coefs],'single');
  for indx = 1:n_coefs
    Zval(:,:,indx) = reshape(zernfun2(indx+3,rho(:),theta(:)),size(rho)) .* H;
  end
end

function phi = Zcoefs2phi(Zcoefs,Zval)
  n_coefs = length(Zcoefs);
  sz = size(Zval);
  phi = zeros(sz(1:end-1),'single');
  for indx = 1:n_coefs
    phi = phi + Zcoefs(indx)*Zval(:,:,indx);
  end
end

function [val,grad] = pdpenalty_zernike(Zcoefs,im,H,pdk,Zval)
  phi = Zcoefs2phi(Zcoefs,Zval);
  if (nargout < 2)
    val = pdpenalty(phi,im,H,pdk);
  else
    [val,gradpix] = pdpenalty(phi,im,H,pdk);
    % Compute the projection of the gradient onto the Zernike
    % coefficients
    n_coefs = length(Zcoefs);
    grad = zeros(1,n_coefs);
    for indx = 1:n_coefs
      tmp = gradpix .* Zval(:,:,indx);
      grad(indx) = sum(tmp(:));
    end
  end
end
