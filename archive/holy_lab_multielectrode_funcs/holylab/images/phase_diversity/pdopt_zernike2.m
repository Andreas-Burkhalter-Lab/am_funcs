function [Zcoefs,psi] = pdopt_zernike2(im,H,rho,theta,Zcoefs)

% this function takes in two images of unknown defocus difference and
% calculates the unknown zerneke aberrations
%
% INPUT
% im = 2-images of unknown defocuses images
% H = pupil function
% rho, theta = positions on the image
% Zcoefs = starting guess for the various coefs
% Typically we are trying to calculate the first 8 zernekes starting from
% defocus (z = 4-11)
% Since the two images have two unknown defocuses and the same unknown
% higher order aberrations, the size of Zcoefs should be 8 + #images-1;
% typically 8 +1 = 9;
%
% OUTPUT
% Zcoefs = optimized Zcoefs
% phi = optimized phi (including defocus)
%
% this script is currently hard-coded to work only for 2 unknown images, i.e. K = 2.
% Since that is the most likely use of the script, it should be fine
% Diwakar Turaga, 8-20-2008
%
% See also: PDPENALTY2, PD_DEMO2, PDOPT_ZERNIKE 



  Zval = zernike_values(rho,theta,H,length(Zcoefs));
  optfun = @(Zc) pdpenalty_zernike(Zc,im,H,Zval);
  %minops = optimset('GradObj','on','DerivativeCheck','on','DiffMinChange',1e-4,'Display','iter');
  minops = optimset('GradObj','on','DiffMinChange',1e-4,'Display','iter', 'MaxIter', 50);
  Zcoefs = fminunc(optfun,Zcoefs,minops);
  if (nargout > 1)
    psi = Zcoefs2psi(Zcoefs,Zval);
  end
end

function Zval = zernike_values(rho,theta,H,n_coefs)
% Compute the Zernike polynomials
% The first 3 Zernike's are omitted, since there is nothing in the data
% that allows us to determine them (and they do not affect image
% quality, only image position)
  rhosz = size(rho);
  rho(rho > 1) = 0;
  %Zval = zeros([rhosz n_coefs],'single');
  Zval = zeros([rhosz n_coefs]);

  % unknown aberration of the system
  for indx = 1:n_coefs-1
    Zval(:,:,indx) = reshape(zernfun2(indx+3,rho(:),theta(:)),size(rho)) .* H;
  end
end



function [val,grad] = pdpenalty_zernike(Zcoefs,im,H,Zval)

  psi = Zcoefs2psi(Zcoefs,Zval);
  if (nargout < 2)
    val = pdpenalty2(psi,im,H);
  else
    [val,gradpix] = pdpenalty(psi,im,H);
    n_coefs = length(Zcoefs);
    
    grad = zeros(1,n_coefs);
    
    tmp = gradpix(:,:,1).*Zval(:,:,1);
    grad(1) = sum(tmp(:));
    
    tmp = gradpix(:,:,2).*Zval(:,:,1);
    grad(2) = sum(tmp(:));
    
    for indx = 3:n_coefs
        tmp1 = gradpix(:,:,1).*Zval(:,:,indx-1);
        tmp2 = gradpix(:,:,2).*Zval(:,:,indx-1);
        
        grad(indx) = sum(tmp1(:)) + sum(tmp2(:));
    end
    
  end
end


function psi = Zcoefs2psi(Zcoefs,Zval)

no_coeffs = length(Zcoefs);

psi(:,:,1) = Zcoefs(1).*Zval(:,:,1);
psi(:,:,2) = Zcoefs(2).*Zval(:,:,1);


for indx = 3:no_coeffs
    psi(:,:,1) = psi(:,:,1) + Zcoefs(indx).*Zval(:,:,indx-1);
    psi(:,:,2) = psi(:,:,2) + Zcoefs(indx).*Zval(:,:,indx-1);
end

end
