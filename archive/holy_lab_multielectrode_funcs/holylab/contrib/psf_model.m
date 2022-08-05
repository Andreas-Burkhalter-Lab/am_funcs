function psf = psf_model(lambda, numAper, rindexObj, dxy, dz, ...
    xysize, nslices, varargin)
% Generate the unaberrated or aberrated PSF using the P.A. Stokseth model
%
% All units must be consistent, i.e., either um or nm.
%
% Inputs:
% lambda: Excitation/Emission Wavelength
% numAper: Numerical aperture of the microscope
% rindexObj: objective's intended immersion refractive index
% dxy: Pixel sizes referenced to the object plane
% dz: Axial Optical axis Z sampling or slicing width
% xysize: Size of the desired image (# of pixels)
% nslices: Number of slices desired
% Optional inputs:
% rindexSp: actual immersion refractive index (default: rindexObj)
% depth: Depth under the coverslip (default: 0)
%
% The output is a 3-array containing psf(x,y,z), with the origin at the
% center of the array.

% P. A. Stokseth (1969). Properties of a defocused optical system. 
% J. Opt. Soc. Am. A 59:1314-1321. 

% Extracted and modified from "wfmpsf" by Praveen Pankajakshan
% Timothy E. Holy, 2012

  % Parsing
  rindexSp = rindexObj;
  depth = 0;
  if ~isempty(varargin)
    rindexSp = varargin{1};
  end
  if length(varargin) > 1
    depth = varargin{2};
  end
  
  % initializing
  pupil = zeros(xysize, xysize);
  psf = pupil;
  
  kxcord = (1:xysize) - (xysize+1)/2;
  defocus = ((1:nslices) - (nslices+1)/2)*dz;
  
  % Pupil space pixel dimensions dkx, dky
  dkxy = (2*pi)/(xysize*dxy);
    
  %Calculated the wavelength of light inside the objective lens and specimen
  lambdaObj = lambda/rindexObj;
  lambdaSp = lambda/rindexSp;
  
  % Calculate the wave vectors in vaccuum, objective and specimens
  k0 = 2*pi/lambda;
  kObj = 2*pi/lambdaObj;
  kSp = 2*pi/lambdaSp;
  
  % Radius of the pupil function disk
  %kMax = 4*xysize*((dxy*numAper)/lambda)^2;
  kMax = (2*pi*numAper)/(lambda*dkxy);
  % Generate the pupil function amplitude
  kycord = kxcord;
  [kx, ky] = meshgrid(kxcord, kycord);
  k = sqrt(kx.^2+ky.^2);
  pupil = (k< kMax);
  
  % Calculate the sine of the semi-aperture angle in the objective lens
  sinthetaObj = (k.*(dkxy))/kObj;
  sinthetaObj(sinthetaObj>1) = 1;
  % Calculate the cosine of the semi-aperture angle in the objective lens
  costhetaObj = eps+sqrt(1-(sinthetaObj.^2));
  
  % Calculate the sine of the semi-aperture angle in the specimen
  sinthetaSp = (k.*(dkxy))/kSp;
  sinthetaSp(sinthetaSp>1) = 1;
  % Calculate the cosine of the semi-aperture angle in the specimen
  costhetaSp = eps+sqrt(1-(sinthetaSp.^2));
  
  % Defocus Phase calculation
  phid = (sqrt(-1)*kObj).*costhetaObj;
  % Spherical aberration phase calculation
  phisa = (sqrt(-1)*k0*depth).*((rindexSp.*costhetaSp)-(rindexObj.*costhetaObj));
  % Calculate the optical path difference due to spherical aberrations
  OPDSA = exp(phisa);
  
  for zk = 1:nslices
    OPDDefocus = exp(defocus(1, zk).*phid);
    pupilDefocus = pupil.*OPDDefocus;
    pupilSA = pupilDefocus.*OPDSA;
    % Calculate the coherent PSF by using inverse Fourier Transform
    psf(:, :, zk) = ifft2(pupilSA);
    % Calculate the incoherent PSF from the coherent PSF
    psf(:, :, zk) = fftshift(abs(psf(:, :, zk)).^2);
  end
%   psf = sqrt(psf);
end