function [out,phiv] = biharmonic_phi(varargin)
% Syntax:
%   internalData = biharmonic_phi(pupildata,act_xy) % first call
%   [phi,phiv] = biharmonic_phi(params,v,internalData)  % subsequent calls
%
% params has the following fields:
%   R: radius at which the membrane is clamped, in grid-coordinates
%   v2f: the conversion of voltage to force
%   Zdefoc (optional, default 0): the coefficient of the sample's defocus Zernike (Z(2,0))
%   scale: a 2-vector, applied first to each grid coordinate in a
%     grid->phase transformation
%   rot: the angle of the rotation
%   shift: a 2-vector
%
% phi is the aberration (including membrane shape + defocus)
% phiv is the membrane shape (without any contribution from defocus)
  
  if (nargin == 2)
    pupildata = varargin{1};
    act_xy = varargin{2};
    inpupil = pupildata.H0 > 0;
    rho = pupildata.rho(inpupil);
    theta = pupildata.theta(inpupil);
    internalData.X = [rho.*cos(theta),rho.*sin(theta)];
    Z4 = zernike_values(pupildata.rho,pupildata.theta,4);
    internalData.Z4 = Z4(inpupil);
    internalData.inpupil = inpupil;
    internalData.imsz = size(pupildata.rho);
    internalData.act_xy = act_xy;
    out = internalData;
    return
  end
  params = varargin{1};
  v = varargin{2};
  internalData = varargin{3};
  
  ct = cos(params.rot); st = sin(params.rot);
  A = [ct -st; st ct] * diag(params.scale,0);
  T = maketform('affine',[A; params.shift]);
  
  % Convert pupil coordinates to coordinates within the unit circle (with
  % respect to clamp coordinates)
  XT = tforminv(T,internalData.X)/params.R;
  XT2 = sum(XT.^2,2);
  if any(XT2 > 1)
    error('Some coordinates are outside the clamped region');
  end
  x0 = internalData.act_xy / params.R;
  
  % Calculate the membrane shape as a superposition of pointsource terms
  % at the different grid points
  f = v * params.v2f;
  nzIndex = find(f ~= 0);
  u = zeros(size(XT,1),1);
  for iter = nzIndex
    u = u + f(iter)*biharmonic_pointsource(XT,x0(iter,[2 1]));
  end
  
  % Pack this into phi
  phiv = zeros(internalData.imsz);
  phiv(internalData.inpupil) = u;
  out = phiv;
  if isfield(params,'Zdefoc')
    out(internalData.inpupil) = phiv(internalData.inpupil) + ...
      params.Zdefoc * internalData.Z4;
  end
