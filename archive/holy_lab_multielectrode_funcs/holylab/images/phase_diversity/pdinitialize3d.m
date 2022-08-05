function [H,pdk,rho,theta] = pdinitialize3d(NA,lambda,n, imsz,pixel_spacing,zk)
  % NA is numerical aperture of the objective, n is
  % the refractive index of the immersion medium.
  % lambda is wavelength of light
  % imsz is the number of pixels along each coordinate
  % pixel_spacing should be supplied in physical (image-space) units,
  % e.g., 4 microns for a camera with 4 micron pixels.
  % zk is the actual amount of focal displacement, measured in image-space
  % units
  %
  % NOTE: I have used an arbitrary aberration and not defocus to make pdk.
  % This is due to the fact that 3d psf itself has defocus in it, and I
  % haven't figured out how applied defocus aberration affects 3d psf
  % calculations
  %
  % Set up the relationships between physical space, fourier space, and
  % the aperture function
  %
  % See also: PDPENALTY3D, PD_DEMO3D, PDOPT_ZERNIKE3D 
 
  
  n_dims = 2;
  coords = cell(1,n_dims);
  for indx = 1:n_dims
    coords{indx} = 0:imsz(indx)-1;
  end
  X = cell(1,n_dims);  % coordinates in the fourier plane
  [X{:}] = ndgrid(coords{:});
  for indx = 1:n_dims
    % Set up circular boundary conditions with origin at corner
    X{indx} = modwrap(X{indx},imsz(indx))/(imsz(indx)*pixel_spacing(indx));
  end
  [theta,rho] = cart2pol(X{2},X{1});
  rhopupil = NA/lambda;
  
  theta = double(theta);
  rho = double(rho);
  
%   theta = single(theta);
%   rho = single(rho);
  
  H = rho.^2 <= rhopupil^2;
  %H = single(H);
  H = double(H);
  
  rholim = rho;
  rholim(rho > 1) = 0;
  abr = reshape(zernfun2(5,rholim(:),theta(:)),size(rho)); % a Z(3,-1) aberration
  pdk = zeros([size(rho) length(zk)],class(H));
colons = repmat({':'},1,n_dims);

for indx = 1:length(zk)
    pdk(colons{:},indx) = H .* abr .* zk(indx); % this is an arbitrary aberration...NOT defocus. 
    % this is due to fact that 3d psf assumes defocus aberration, so I
    % haven't figured out how to combine both 3d psf and defocus aberration
end


end

function xo = modwrap(xi,modval)
% xo = xi % modval, except return in range [-modval/2,modval/2]
  xo = mod(xi,modval);
  toobig = xo > modval/2;
  xo(toobig) = xo(toobig) - modval;
end