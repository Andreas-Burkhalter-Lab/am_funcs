function psf = psf_preset(instrument_params,sheet_params,lambda,ds,nxy,ns)
% psf_preset: OCPI microscopy point-spread functions
%
% This calculates the theoretical point-spread function for different
% OCPI microscopes at WashU.
%
% Syntax:
%   psf = psf_preset(instrument_params,sheet_params,lambda,ds,nxy,ns)
% where
%   instrument_params is a structure of the format produced by ocpi_params;
%   sheet_params is a structure specifying the parameters of the light
%     sheet, with fields
%       sigma: the sigma for the Gaussian beam specifying the light
%         sheet (in microns) (note sigma = FWHM/2.355).
%       shift (default 0): the amount of defocus, in microns, of the center
%         of the light sheet (the sign is meaningful)
%   lambda is the emission wavelength in nm (e.g., 509 for GFP, etc; as a
%     substitute, use the excitation wavelength, e.g., 488)
%   ds is the displacement along the scanning axis between frames (microns)
%   nxy is the number of pixels along the x- and y-axes
%   ns is the number of frames along the scanning axis
% and
%   psf is a 3-array of the psf(x, y, z). For microscopes with the
%     scan axis not aligned to the detection optical axis, the final
%     size of the array may not be [nxy nxy ns].
%
% See also: psf_model, ocpi_params.

% Copyright 2012 by Timothy E. Holy

  % Check that lambda is sensible
  if lambda < 350 || lambda > 1000
    error('lambda must be supplied in nm')
  end
  rindexObj = opt_refrindx('water',lambda);
  lambda = lambda/1000;  % convert to um, like the rest of the spatial quantities
  
  instrument_lookup = {'OCPI1','HS-OCPI'};
  default_sigma = [-1, 2.47];   % OCPI1 not currently calibrated; HS-OCPI for 488
  sheet_params = default(sheet_params,'sigma',default_sigma(strcmp(instrument_params.instrument,instrument_lookup)),'shift',0);
  
  %% Oversample along z and then sum
  dxy = instrument_params.pixel_xy_microns;
  dz = ds*cos(instrument_params.angle_scan_optax); % spacing along the z-axis
  r = max(1,ceil(dz/dxy));  % oversampling ratio
  dz = dz/r;
  nz = ns*r;
  
  %% Calculate the psf in orthogonal units
  % The detection psf
  psf = psf_model(lambda, instrument_params.numAper, rindexObj, dxy, dz, nxy, nz);
  % Illumination PSF
  centerz = (nz+1)/2;
  if nargin > 7
    centerz = centerz + sheet_params.shift/dz;
  end
  z = ((1:nz)-centerz)*dz/sheet_params.sigma;
  g = exp(-z.^2/2);
  psf = bsxfun(@times,psf,reshape(g,[1 1 nz]));
  
  %% Correct for non-orthogonal scan geometry
  if instrument_params.angle_scan_optax ~= 0
    T = maketform('affine',[1 -(ds/r)/dxy*sin(instrument_params.angle_scan_optax); 0 1; 0 0]);  % note (2,2) term should _not_ be cos(theta) (cancels with ds/dz)
    Tf = fliptform(T);
    ppsf = permute(psf,[2 3 1]);
    ppsf = imtransform(ppsf,Tf);
    psf = ipermute(ppsf,[2 3 1]);
  end
  
  %% Undo oversampling
  sz = size(psf);
  psf = reshape(psf,[sz(1:2) r sz(3)/r]);
  psf = squeeze(sum(psf,3));
  
  % Normalize
  psf = psf / sum(psf(:));
  