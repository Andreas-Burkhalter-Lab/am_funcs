% Example of two-dimensional phase diversity.
%
% This is an example script which produces two aberrated images differing
% by a defocus. From these, it estimates the underlying object and the
% phase aberration.
% 
% See also: PDPENALTY, PDOPT_ZERNIKE2.

% Copyright 2009 by Timothy E. Holy & Diwakar Turaga

%% Prepare the underlying object
% Use this next line for a "real" image
im = imread('cameraman.tif');

% Or, uncomment the next 2 lines for a single "bead" (point)
%im = zeros([256 256],'single');
%im(128,128) = 1;

% Convert to a floating-point type
%im0 = single(im);
im0 = double(im);

% Set how much noise you want in the acquired images. 0.01 would correspond
% to noise that is 1% of the mean intensity.
% noisefactor = 0.01;
noisefactor = 0;

%% Specify the parameters of the imaging system
NAred = 0.5/1.34;  % NAred = NA / n, NA = Numerical Aperture, n = refractive index of immersion fluid
M = 20;   % magnification
lambda = 0.5;   % wavelength (the units here are microns, whatever you choose needs to be the same everywhere else)
defoc = 20*lambda;   % Amount that one image is defocused by, in terms of the object defocus (the camera would be shifted by M^2*defoc)
pixel_spacing = [4 4];   % Spacing between pixels in the detector, in the same units as wavelength

%% Initialize the objects needed for phase diversity
% H = the aperture function
% pdk = the 2 known phase functions (the first phase will zero, the second
%       the amount of defocus in the 2nd image)
% rho = the normalized radius within the aperture (imagesc(H.*rho) will
%       show you the radius; it's "wrapped" due to Fourier transform
%       issues)
% theta = the angle for positions within the aperture
[H,pdk,rho,theta] = pdinitialize(NAred,M,lambda,size(im),pixel_spacing,[0 defoc]*M^2);

%% Calculate and display the "diffraction-limited" image
h = ifft2(H);
H2 = fft2(h .* conj(h));  % FT of PSF for incoherent illumination
imfft = fft2(im0);
imdl = ifft2(imfft .* H2);

figure; imshowsc(im); title('Object')
figure; imshowsc(imdl); title('Diffraction-limited image')

%% Define and display a phase aberration
rholim = rho;
rholim(rho > 1) = 0;
phi = reshape(zernfun2(7,rholim(:),theta(:)),size(rho)) .* H; % a Z(3,-1) aberration
phi = 2*phi;

% phi2 = reshape(zernfun2(9,rholim(:),theta(:)),size(rho)) .* H;
% phi2 = -5*phi2;

figure; imagesc(phi); title('Phase aberration'); colorbar

%% Calculate and display the observed images
% Unlike the diffraction-limited image, this incorporates the
% aberration phi and any "aberrations" corresponding to the diversity
% images
imk = zeros(size(pdk),class(im0));
K = size(pdk,3);
for indx = 1:K
  Hk = H .* exp(1i * (phi + pdk(:,:,indx)));  % The coherent PSF in Fourier coords
  hk = ifft2(Hk);  % The coherent PSF in image coordinates
  Hk2 = fft2(hk .* conj(hk));   % The incoherent PSF in Fourier coords
  imk(:,:,indx) = ifft2(imfft .* Hk2);  % The noiseless image
  if (indx == 1)
    % Define the noise in terms of the observed image
    tmp = imk(:,:,1);
    noiseamplitude = noisefactor * mean(tmp(:));
  end
  imk(:,:,indx) = imk(:,:,indx) + noiseamplitude*randn(size(im));
end

figure; for indx = 1:K; subplot(1,K,indx); imshowsc(imk(:,:,indx)); end; suptitle('Diversity images')


%% Set up the initial guess for fitting the aberration
% Zcoefs = zeros(1,9); % Two image-specific defocus and unknown aberrations from Z(2,0) up to Z(4,0)
Zcoefs = rand(1,9)+0.1; % random starting guess seems to work the best. Something about starting from zeros makes it get stuck in local minima

%% Optimize the fit
% This uses pdopt_zernike2; a more "modern" and flexible approach is to use
% phi_parametrizations and calcphi2d.
[Zcoefs,psi] = pdopt_zernike2(imk,H,rho,theta,Zcoefs);
Zcoefs

%% Display the results
[val,grad,object] = pdpenalty(psi,imk,H);  % needed to get the object
% Note: if you have noise, the object will look horrible. But try smoothing
% it (e.g., with imfilter_gaussian), you will see that you've gotten
% something useful out.
figure; imshowsc(object); title('Estimate of original object')
figure; imagesc(phi); title('Estimate of aberration')

