function [Zcoefs,phi] = pdopt_zernike3d(im,H,pdk,rho,theta,Zcoefs, Hz)

% this function takes in two stacks of unknown aberrations and
% calculates the unknown zerneke aberrations
%
% INPUT
% im = 2-stacks of aberrated images
% H = pupil function
% rho, theta = positions on the image
% Zcoefs = starting guess for the various coefs
% Typically we are trying to calculate the first 8 zernekes starting from
% defocus (z = 4-11)
% Since the two images have two unknown defocuses and the same unknown
% higher order aberrations, the size of Zcoefs should be 8 + #images-1;
% typically 8 +1 = 9;
% Hz = 3d pupil function
%
% OUTPUT
% Zcoefs = optimized Zcoefs
% phi = optimized phi (including defocus)
%
% this script is currently hard-coded to work only for 2 unknown images, i.e. K = 2.
% Since that is the most likely use of the script, it should be fine
% Diwakar Turaga, 
%
% See also: PDPENALTY3D, PD_DEMO3D, PDOPT_ZERNIKE2 

  Zval = zernike_values(rho,theta,H,length(Zcoefs));
  optfun = @(Zc) pdpenalty_zernike(Zc,im,H,pdk,Zval, Hz);
  %minops = optimset('GradObj','on','DerivativeCheck','on','DiffMinChange',1e-4,'Display','iter');
  minops = optimset('GradObj','on','DiffMinChange',1e-4,'Display','iter', 'MaxIter', 30);
  %minops = optimset('GradObj','on','DiffMinChange',1e-4,'Display','iter');
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
  Zval = zeros([rhosz n_coefs],class(H));
  for indx = 1:n_coefs
    Zval(:,:,indx) = reshape(zernfun2(indx+3,rho(:),theta(:)),size(rho)) .* H;
  end
end

function phi = Zcoefs2phi(Zcoefs,Zval)
Zcoefs
  n_coefs = length(Zcoefs);
  sz = size(Zval);
  phi = zeros(sz(1:end-1),class(Zval));
  for indx = 1:n_coefs
    phi = phi + Zcoefs(indx)*Zval(:,:,indx);
  end
end

function [val,grad] = pdpenalty_zernike(Zcoefs,im,H,pdk,Zval, Hz)
  phi = Zcoefs2phi(Zcoefs,Zval);
  if (nargout < 2)
    val = pdpenalty3d(phi,im,H,pdk, Hz);
  else
    [val,gradpix] = pdpenalty3d(phi,im,H,pdk, Hz);
    % Compute the projection of the gradient onto the Zernike
    % coefficients
    n_coefs = length(Zcoefs);
    grad = zeros(1,n_coefs);
    for indx = 1:n_coefs
      Zval3d = repmat(Zval(:,:,indx), [1 1 size(gradpix,3)]);
      tmp = gradpix .* Zval3d;
      grad(indx) = sum(tmp(:));
    end
  end
end