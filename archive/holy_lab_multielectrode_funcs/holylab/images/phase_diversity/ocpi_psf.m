function [psf_fm, psf_cyl, psf_ocpi] = ocpi_psf(imsz, zFrame, pixel_spacing, Zcoefs, n, lambda_det, NA_det, lambda_ill, NA_grin, laser_e2, slit_x, grin_thickness)

% calculating 3d psf and otf
% adapted from hanser et al, J. Microscopy, vol 216, pg 32, year 2004.
% see ocpi_psf_script for details


zFrame = abs(zFrame);
nFrame = length(zFrame);
%----------------------------------------------
% First dealing with the illumination to get the cylindrical component of
% the ocpi psf

slit_red = slit_x/(laser_e2/2);
NA_ill = (slit_x/grin_thickness)*NA_grin;

n_dims = 2;
coords = cell(1,n_dims);
for indx = 1:n_dims
    coords{indx} = 0:imsz(indx)-1;
end
X = cell(1,n_dims);  % coordinates in the fourier plane
[X{:}] = ndgrid(coords{:});

for indx = 1:n_dims
    % Set up circular boundary conditions with origin at corner

    xo = mod(X{indx},imsz(indx));
    toobig = xo > imsz(indx)/2;
    xo(toobig) = xo(toobig) - imsz(indx);

    X{indx} = xo/(imsz(indx)*pixel_spacing(indx));
    %X{indx} = xo*pixel_spacing(indx);
end


[theta,rho] = cart2pol(X{2},X{1});
rhopupil = NA_ill/lambda_ill;

H_ill = rho.^2 <= rhopupil^2;
H_ill = unwrap(H_ill);

mid_x = round(imsz(1)/2);
laser_sigma = round(sum(H_ill(mid_x,:))/2);

H = zeros(imsz);
H(mid_x,mid_x) = 1;
Hin = imfilter_gaussian(single(H), [laser_sigma laser_sigma]);

laser_sigma_slit = round(laser_sigma * slit_red);
Hcyl = zeros(size(Hin));
Hcyl(mid_x - laser_sigma_slit : mid_x +laser_sigma_slit, :) = 1;

Hincyl = Hin.*Hcyl;

hcyl = ifft2(Hincyl);
psf_cyl = hcyl.*conj(hcyl);
psf_cyl = unwrap(psf_cyl);

psf_cyl_1d = squeeze(psf_cyl(:, round(imsz(1)/2)));

psf_cyl_1d = imresize(psf_cyl_1d, [nFrame 1]);

%-------------------------------------------------------------------------
% Now to deal with the detector

% Set up the relationships between physical space, fourier space, and
% the aperture function
n_dims = 2;
coords = cell(1,n_dims);
for indx = 1:n_dims
    coords{indx} = 0:imsz(indx)-1;
end
X = cell(1,n_dims);  % coordinates in the fourier plane
[X{:}] = ndgrid(coords{:});

for indx = 1:n_dims
    % Set up circular boundary conditions with origin at corner

    xo = mod(X{indx},imsz(indx));
    toobig = xo > imsz(indx)/2;
    xo(toobig) = xo(toobig) - imsz(indx);

    X{indx} = xo/(imsz(indx)*pixel_spacing(indx));
    %X{indx} = xo*pixel_spacing(indx);
end


[theta,rho] = cart2pol(X{2},X{1});
rhopupil = NA_det/lambda_det;

H = rho.^2 <= rhopupil^2;

Zval = zernike_values(rho, theta, H, length(Zcoefs));
phi = Zcoefs2phi(Zcoefs,Zval);
Hk = H .* exp(i * phi);


%---------------------------------------------
% calculating 3d psf sedat style

coef = n/lambda_det;
kz = (coef^2 - rho.^2).^0.5;

%Hz = zeros([size(H), nFrame], 'single');
psf_fm = zeros([size(Hk), nFrame], 'single');
psf_ocpi = psf_fm;

for indx = 1:nFrame
    Hz = Hk .* exp(i * 2 * pi * zFrame(indx) .* kz);
    hz = ifft2(Hz);
    sz = hz .* conj(hz);
    psf_fm(:,:,indx) = unwrap(sz);
    psf_ocpi(:,:,indx) = psf_fm(:,:,indx)*psf_cyl_1d(indx);
end


end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function xo = unwrap(xi)

imsize = size(xi);

xo(1:imsize(1)/2, 1:imsize(2)/2) = imrotate(xi(1:imsize(1)/2, 1:imsize(2)/2), 180);
xo(1:imsize(1)/2, imsize(2)/2+1:imsize(2)) = imrotate(xi(1:imsize(1)/2, imsize(2)/2+1:imsize(2)), 180);
xo(imsize(1)/2+1:imsize(1), 1:imsize(2)/2) = imrotate(xi(imsize(1)/2+1:imsize(1), 1:imsize(2)/2), 180);
xo(imsize(1)/2+1:imsize(1), imsize(2)/2+1:imsize(2)) = imrotate(xi(imsize(1)/2+1:imsize(1), imsize(2)/2+1:imsize(2)), 180);

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
    Zval(:,:,indx) = reshape(zernfun2(indx+3,rho(:),theta(:), 'norm'),size(rho)) .* H;
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
    