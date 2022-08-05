function sfunc = phi_parametrizations(config,pupildata)
% PHI_PARAMETRIZATIONS: a general "switchyard" for different phase models
% The general syntax is
%   sfunc = phi_parametrizations(config)
% OR
%   sfunc = phi_parametrizations(config,pupildata)
% where sfunc provides param2phi and gphi2gparam as required by calcphi2d.
% 'config' is a structure with a field titled 'mode' containing a string
% which specifies the particular type of phase model. Some examples are:
%   'register pair': for rigid translation of two images in terms of
%     first-order Zernike "aberrations"
%   'register pair defocus': for rigid translation + focus adjustment
%   'registered vraw craw unab-ab': a series of unaberrated and aberrated
%     images, where the aberrated path has a constant aberration (c) plus a
%     voltage-dependent (v) aberration. Both aberrations are represented as
%     raw phases. This also allows you to supply results from a previous
%     registration.
% and so on.
%
% It is intended that you can deposit more phase models in this function.
% The only real way to understand this function is to dive into the code!
%
% See also: CALCPHI2D.

% Copyright 2009 by Timothy E. Holy

  if (nargin > 1)
    config.H0 = pupildata.H0;
    config.rho = pupildata.rho;
    config.theta = pupildata.theta;
  end

  % Note: in cases below, "pair" has been generalized to include a first
  % unaberrated image and any number of aberrated images
  switch config.mode
    case 'register pair'
      % Use this if you want to register a pair of images in terms of
      % coefficients of the two first-order Zernike basis functions
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      sfunc.param2phi = @(p) pp_register_pair_expand(p,zern1);
      sfunc.gphi2gparam = @(grad) pp_register_pair_contract(grad,zern1);

    case 'register defocus pair'
      % Use this if you want to register a pair of images in terms of
      % coefficients of the two first-order Zernike basis functions,
      % plus allow there to be some defocus
      zern1 = zernike_values(config.rho,config.theta,[1 2 4]);
      sfunc.param2phi = @(p) pp_register_pair_expand(p,zern1);
      sfunc.gphi2gparam = @(grad) pp_register_pair_contract(grad,zern1);

    case 'Zn pair'
      % Use this if you want to "register" a pair of images in terms of
      % coefficients of an arbitrary set Zernike polynomials
      zern1 = zernike_values(config.rho,config.theta,config.Zindex);
      sfunc.param2phi = @(p) pp_register_pair_expand(p,zern1);
      sfunc.gphi2gparam = @(grad) pp_register_pair_contract(grad,zern1);

    case 'registered second zernike pair'
      % Use this if you don't need to tweak the registration but want to fit
      % second order Zernike polynomials
      zern1 = zernike_values(config.rho,config.theta,1:5);
      sfunc.param2phi = @(p) pp_registered_pair_expand(p,zern1);
      sfunc.gphi2gparam = @(grad) pp_registered_pair_contract(grad,zern1);
       
    case 'Zn block'
      % Use this if you want to compute Zernikes for a block of images with
      % a single underlying object
      zern1 = zernike_values(config.rho,config.theta,config.Zindex);
      sfunc.param2phi = @(p) pp_register_block_expand(p,zern1);
      sfunc.gphi2gparam = @(grad) pp_register_block_contract(grad,zern1);

    case 'registered vdep raw singlechannel'
      % This is a case where you assume the registration is already done (&
      % parametrized in terms of Z1 coefficients), and you just want to
      % find the voltage-dependent raw phase
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      if ~isfield(config,'Zshifts')
        config.Zshifts = zeros(length(config.v),2);
      end
      sfunc.param2phi = @(p) pp_register_vraw_expand(p,config.Zshifts,zern1,config.v);
      sfunc.gphi2gparam = @(grad) pp_register_vraw_contract(grad,config.v);

    case 'register vraw singlechannel'
      % Do both registration and a voltage-dependent raw phase for a
      % single imaging channel.
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      sfunc.param2phi = @(p) pp_doregister_vraw_expand(p,zern1,config.v,config.Zmultiplier);
      sfunc.gphi2gparam = @(grad) pp_doregister_vraw_contract(grad,zern1,config.v,config.baseIndex,config.Zmultiplier);

    case 'vgaussian singlechannel'
      % Voltage-dependent gaussian aberration
      isInPupil = config.rho <= 1;
      rho = config.rho(isInPupil(:));
      theta = config.theta(isInPupil(:));
      X = [rho .* cos(theta), rho .* sin(theta)];
      [p2f,gf2gp] = chainrule_gaussian(X,isInPupil);
      sfunc.param2phi = @(p) pp_vgaussian_expand(p,config.v,p2f);
      sfunc.gphi2gparam = @(grad) pp_vgaussian_contract(grad,config.v,gf2gp);
      sfunc.phicoef = @(p) p2f(p);

    case 'register vgaussian singlechannel'
      % Do both registration and a voltage-dependent raw phase for a
      % single imaging channel.
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      isInPupil = config.rho <= 1;
      rho = config.rho(isInPupil(:));
      theta = config.theta(isInPupil(:));
      X = [rho .* cos(theta), rho .* sin(theta)];
      [p2f,gf2gp] = chainrule_gaussian(X,isInPupil);
      sfunc.param2phi = @(p) pp_doregister_vgaussian_expand(p,zern1,config.v,config.Zmultiplier,p2f);
      sfunc.gphi2gparam = @(grad) pp_doregister_vgaussian_contract(grad,zern1,config.v,config.baseIndex,config.Zmultiplier,gf2gp);
      sfunc.phicoef = @(p) p2f(p(2*length(config.v)+1:end));

    case 'registered vdep Zernike singlechannel'
      % This is a case where you assume the registration is already done (&
      % parametrized in terms of Z1 coefficients), but unlike the above the
      % voltage-dependent phase is parametrized in terms of Zernike basis
      % functions
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      zernv = zernike_values(config.rho,config.theta, ...
        config.Zindex);
      sfunc.param2phi = @(p) pp_register_vzern_expand(p,config.Zshifts,zern1,zernv,config.v);
      sfunc.gphi2gparam = @(grad) pp_register_vzern_contract(grad,zernv,config.v);
      sfunc.phicoef = @(p) pp_zern2phi(p,zernv);

     case 'registered vraw unab-ab'
      % This is a case where you assume the registration is already done (&
      % parametrized in terms of Z1 coefficients), and you're looking at
      % both an "unaberrated" channel and an aberrated channel. The
      % aberrated channel has a voltage-dependent aberration component,
      % parametrized as raw.
      % The two channels can also differ by a static registration
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      sfunc.param2phi = @(p) pp_register_vraw2chan_expand(p,config.Zshifts,zern1,config.unabIndex,config.abIndex,config.v,config.Zstatic);
      sfunc.gphi2gparam = @(grad) pp_register_vraw2chan_contract(grad,config.abIndex,config.v);
     
      
    case 'registered vraw craw unab-ab'
      % This is a case where you assume the registration is already done (&
      % parametrized in terms of Z1 coefficients), and you're looking at
      % both an "unaberrated" channel and an aberrated channel. The
      % aberrated channel has a "constant" component of its aberration
      % and a voltage-dependent component (both parametrized raw).
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      sfunc.param2phi = @(p) pp_register_vrawcraw_expand(p,config.Zshifts,zern1,config.unabIndex,config.abIndex,config.v);
      sfunc.gphi2gparam = @(grad) pp_register_vrawcraw_contract(grad,config.abIndex,config.v);

    case 'registered vraw czern unab-ab'
      % This is a case where you assume the registration is already done (&
      % parametrized in terms of Z1 coefficients), and you're looking at
      % both an "unaberrated" channel and an aberrated channel. The
      % aberrated channel has a "constant" component of its aberration
      % (parametrized in terms of Zernikes) and a voltage-dependent
      % component (parametrized raw).
      % First coordinates are the Zernike coefficients, next are raw
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      zernab = zernike_values(config.rho,config.theta, ...
        config.Zindex);
      config = default(config,'Zmultiplier',1);
      zernab = config.Zmultiplier * zernab;
      sfunc.param2phi = @(p) pp_register_vrawczern_expand(p,config.Zshifts,zern1,zernab,config.unabIndex,config.abIndex,config.v,config.framesize);
      sfunc.gphi2gparam = @(grad) pp_register_vrawczern_contract(grad,zernab,config.abIndex,config.v);
      sfunc.phicoef = @(p) pp_zern2phi(p(1:length(config.Zindex)),zernab);
      
    case 'registered vzern czern unab-ab'
      % This is a case where you assume the registration is already done (&
      % parametrized in terms of Z1 coefficients), and you're looking at
      % both an "unaberrated" channel and an aberrated channel. The
      % aberrated channel has a "constant" component of its aberration
      % and a voltage-dependent component, both parametrized in terms of
      % Zernike coefficients.
      % First coefficients are for the static, next are the vdep
      zern1 = zernike_values(config.rho,config.theta,[1 2]);
      zernab = zernike_values(config.rho,config.theta, ...
        config.Zindex);
      config = default(config,'Zmultiplier',1);
      zernab = config.Zmultiplier * zernab;
      sfunc.param2phi = @(p) pp_register_vzernczern_expand(p,config.Zshifts,zern1,zernab,config.unabIndex,config.abIndex,config.v);
      sfunc.gphi2gparam = @(grad) pp_register_vzernczern_contract(grad,zernab,config.abIndex,config.v);
      sfunc.phicoef = @(p) pp_zern2phi(reshape(p,length(config.Zindex),2),zernab);
      
      case 'registered shift-gaussian unab-ab'
          % Registered (fixed) and an aberration that is a shift + a single gaussian hump
          % with amplitude proportional to voltage. Order of parameters:
          % Z1 (shift parameters), amplitude, meanx, meany, omega (inverse std dev)
          % OR Z1, amplitude, meanx, meany, omegaxx, omegaxy, omegayy
          zern1 = zernike_values(config.rho,config.theta,[1 2]);
          isInPupil = config.rho <= 1;
          rho = config.rho(isInPupil(:));
          theta = config.theta(isInPupil(:));
          X = [rho .* cos(theta), rho .* sin(theta)];
          [p2f,gf2gp,gf2gp_safe] = chainrule_gaussian(X,isInPupil);
          sfunc.param2phi = @(p) pp_register_gaussian_expand(p,config.Zshifts,zern1,config.unabIndex,config.abIndex,config.v,p2f);
          sfunc.gphi2gparam = @(grad) pp_register_gaussian_contract(grad,zern1,config.abIndex,config.v,gf2gp);
          sfunc.phicoef = @(p) p2f(p(3:end));

      case 'unab-ab'
          % Two images are input, the first is "unaberrated". So it only
          % contains the flattening parameters. The second is "aberrated".
          % It contains both the fixed flatenning parameters, and the
          % to-be-calculated unknown phase aberration.
          
          zern1 = zernike_values(config.rho,config.theta,config.Zindex);
          sfunc.param2phi = @(p) pp_pair_expand(p,zern1);
          sfunc.gphi2gparam = @(grad) pp_pair_contract(grad,zern1);

      otherwise
          error('mode not recognized');
  end
end
    
%% Linear combinations of static phase "templates" (e.g., Zernikes)
% Coefs are indexed by (diversityIndex,phiBasisIndex)
function phi = pp_linear_static_expand(coefs,phiBasisFix)
  sz = size3(phiBasisFix);
  K = size(coefs,1);
  pf = reshape(phiBasisFix,[sz(1)*sz(2) sz(3)]);
  phi = pf * coefs';
  phi = reshape(phi,[sz(1:2) K]);
end
  
function gradP = pp_linear_static_contract(grad,phiBasisFix)
  gradsz = size3(grad);
  K = gradsz(3);
  phisz = size3(phiBasisFix);
  N = phisz(3);
  gradP = zeros(K,N);
  for divIndex = 1:K
    for basisIndex = 1:N
      tmp = grad(:,:,divIndex) .* phiBasisFix(:,:,basisIndex);
      gradP(divIndex,basisIndex) = sum(tmp(:));
    end
  end
end

%% Some special cases of the former
function phi = pp_register_pair_expand(p,zern1)
  % "Pair" has been generalized to include a first unaberrated image and
  % any number of aberrated images
  if isvector(p)
    coefs = [zeros(1,length(p)); p(:)'];
  else
    coefs = [zeros(1,size(p,2)); p];
  end
  phi = pp_linear_static_expand(coefs,zern1);
end

function gradP = pp_register_pair_contract(grad,zern1)
  gradP = pp_linear_static_contract(grad,zern1);
  gradP = gradP(2:end,:);
end

%% Block functions are like pair except we don't assume first is unaberrated
function phi = pp_register_block_expand(p,zern1)
  phi = pp_linear_static_expand(p,zern1);
end

function gradP = pp_register_block_contract(grad,zern1)
  gradP = pp_linear_static_contract(grad,zern1);
end
  
%% Registration already done, fit higher-order coefficients
function phi = pp_registered_pair_expand(p,zern1)
  coefs = [zeros(1,length(p)); p(:)'];
  phi = pp_linear_static_expand(coefs,zern1); % OK for expand
end

function gradP = pp_registered_pair_contract(grad,zern1)
  gradP = pp_linear_static_contract(grad,zern1);
  gradP = gradP(2,:);
  gradP(1:2) = 0;   % skip the first two, which correspond to registration
end
 
%% Voltage-dependent raw phi, registration already completed
% Neither the registration parameters nor v are considered variable (just
% phiraw)
function phi = pp_register_vraw_expand(raw,coefs,zern1,v)
  %K = length(v);
  phiBasis = zern1;
  %sz = size3(phiBasis);
  phiBasis(:,:,end+1) = raw;
  coefs(:,3) = v(:);
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_register_vraw_contract(grad,v)
  sz = size3(grad);
  K = sz(3);
  gradP = zeros(sz(1:2));
  for k = 1:K
    gradP = gradP + v(k) * grad(:,:,k);
  end
end

%% Voltage-dependent raw phi, also do registration. Single-channel
function phi = pp_doregister_vraw_expand(p,zern1,v,multiplier)
  K = length(v);
  coefs = multiplier*reshape(p(1:2*K),K,2);
  phiBasis = zern1;
  sz = size3(phiBasis);
  phiBasis(:,:,end+1) = reshape(p(2*K+1:end),sz(1:2));
  coefs(:,3) = v(:);
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_doregister_vraw_contract(grad,zern1,v,baseIndex,multiplier)
  sz = size3(grad);
  K = sz(3);
  % Raw phase
  gradRaw = zeros(sz(1:2));
  for k = 1:K
    gradRaw = gradRaw + v(k) * grad(:,:,k);
  end
  % Registration
  gradZ = pp_linear_static_contract(grad,zern1); % project to zern basis
  gradZ = gradZ*multiplier;
  gradZ(baseIndex,:) = 0;
  gradP = [gradZ(:); gradRaw(:)];
end


%% Voltage-dependent gaussian phi, also do registration. Single-channel
function phi = pp_doregister_vgaussian_expand(p,zern1,v,multiplier,p2f)
  K = length(v);
  coefs = multiplier*reshape(p(1:2*K),K,2);
  phiBasis = zern1;
  sz = size3(phiBasis);
  gaussian_params = p(2*K+1:end);
  phiBasis(:,:,end+1) = p2f(gaussian_params);
  coefs(:,3) = v(:);
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_doregister_vgaussian_contract(grad,zern1,v,baseIndex,multiplier,gf2gp)
  sz = size3(grad);
  K = sz(3);
  % Raw phase
  gradRaw = zeros(sz(1:2));
  for k = 1:K
    gradRaw = gradRaw + v(k) * grad(:,:,k);
  end
  gradG = gf2gp(gradRaw);
  % Registration
  gradZ = pp_linear_static_contract(grad,zern1); % project to zern basis
  gradZ = gradZ*multiplier;
  gradZ(baseIndex,:) = 0;
  gradP = [gradZ(:); gradG(:)];
end


%% Voltage-dependent gaussian phi, single-channel
function phi = pp_vgaussian_expand(gaussian_params,v,p2f)
  phiBasis = p2f(gaussian_params);
  coefs = v(:);
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_vgaussian_contract(grad,v,gf2gp)
  sz = size3(grad);
  K = sz(3);
  % Raw phase
  gradRaw = zeros(sz(1:2));
  for k = 1:K
    gradRaw = gradRaw + v(k) * grad(:,:,k);
  end
  gradP = gf2gp(gradRaw);
end


%% Registration plus voltage-dependent Zernike-basis phi
% Only the voltage-dependent part is considered variable
function phi = pp_register_vzern_expand(p,coefs,zern1,zernv,v)
  %K = length(v);
  szzv = size3(zernv);
  phiBasis = zern1;
  phiBasis(:,:,end+1:end+szzv(3)) = zernv;
  coefs = [coefs, v(:) * p(:)'];
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradZ = pp_register_vzern_contract(grad,zernv,v)
  gradP = pp_register_vraw_contract(grad,v); % extracts v-dependent raw phase
  gradZ = pp_linear_static_contract(gradP,zernv); % project to zern basis
%  szzv = size3(zernv);
%  N = szzv(3);
%  gradZ = zeros(N,1);
%  for indx = 1:N
%    tmp = gradP .* zernv(:,:,indx);
%    gradZ(indx) = sum(tmp(:));
%  end
end

%% Aberrated+unaberrated, registration, constant raw aberration
%  plus a voltage-dependent raw aberration
function phi = pp_register_vrawcraw_expand(p,coefsI,zern1, ...
					    unabIndex,abIndex,v)
  % Handle the registration (common to both channels)
  phiBasis = zern1;
  coefs(unabIndex,:) = coefsI;
  coefs(abIndex,:) = coefsI;
  % Include the raw components
  Nab = length(abIndex);
  coefs(abIndex,end+1:end+2) = [ones(Nab,1) v(:)];
  phiBasis(:,:,end+1:end+2) = p;
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_register_vrawcraw_contract(grad,abIndex,v)
  gradV = pp_register_vraw_contract(grad(:,:,abIndex),v); % this extracts v-dependent phase
  gradStatic = pp_register_vraw_contract(grad(:,:,abIndex),ones(1,length(abIndex))); % this extracts static phase
  gradP = gradStatic;
  gradP(:,:,2) = gradV;
end

%% Aberrated+unaberrated, registration, voltage-dependent raw aberration
function phi = pp_register_vraw2chan_expand(p,coefsI,zern1, ...
					    unabIndex,abIndex,v,Zstatic)
  % Handle the registration (common to both channels)
  phiBasis = zern1;
  coefs(unabIndex,:) = coefsI;
  coefs(abIndex,:) = coefsI + repmat(Zstatic,length(abIndex),1);
  % Include the raw voltage-dependent component
  coefs(abIndex,end+1) = v(:);
  phiBasis(:,:,end+1) = p;
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_register_vraw2chan_contract(grad,abIndex,v)
  gradP = pp_register_vraw_contract(grad(:,:,abIndex),v); % this extracts v-dependent phase
end


%% Aberrated+unaberrated, registration, constant Zernike-basis aberration
%  plus a voltage-dependent raw aberration
function phi = pp_register_vrawczern_expand(p,coefsI,zern1,zernab, ...
					    unabIndex,abIndex,v,framesize)
  % Handle the registration (common to both channels)
  phiBasis = zern1;
  coefs(unabIndex,:) = coefsI;
  coefs(abIndex,:) = coefsI;
  % Include the zernike aberration components (constant, ie voltage-indep)
  szzv = size3(zernab);
  Nz = szzv(3);
  Nab = length(abIndex);
  phiBasis(:,:,end+1:end+Nz) = zernab;
  coefs(abIndex,end+1:end+Nz) = repmat(reshape(p(1:Nz),1,Nz),Nab,1);
  % Include the raw components
  raw = reshape(p(Nz+1:end),framesize);
  phiBasis(:,:,end+1) = raw;
  coefs(abIndex,end+1) = v(:);
  phi = pp_linear_static_expand(coefs,phiBasis); % OK for expand
end

function gradP = pp_register_vrawczern_contract(grad,zernab,abIndex,v)
  gradRaw = pp_register_vraw_contract(grad(:,:,abIndex),v); % this extracts v-dependent phase
  gradStatic = pp_register_vraw_contract(grad(:,:,abIndex),ones(1,length(abIndex))); % this extracts static phase
  gradZ = pp_linear_static_contract(gradStatic,zernab); % project to zern basis
  gradP = [gradZ(:); gradRaw(:)];
end

%% Aberrated+unaberrated, registration, constant Zernike-basis aberration
%  plus a voltage-dependent Zernike aberration
function phi = pp_register_vzernczern_expand(p,coefsI,zern1,zernab, ...
					    unabIndex,abIndex,v)
  % Handle the registration (common to both channels)
  phiBasis = zern1;
  coefs(unabIndex,:) = coefsI;
  coefs(abIndex,:) = coefsI;
  % Include the zernike aberration components
  szzv = size3(zernab);
  Nz = szzv(3);
  Nab = length(abIndex);
  phiBasis(:,:,end+1:end+Nz) = zernab;
  % Do the constant + voltage-dependent aberration
  coefs(abIndex,end+1:end+Nz) = repmat(reshape(p(1:Nz),1,Nz),Nab,1) + ...
    v(:)*reshape(p((1:Nz)+Nz),1,Nz);
  phi = pp_linear_static_expand(coefs,phiBasis);
end

function gradP = pp_register_vzernczern_contract(grad,zernab,abIndex,v)
  gradVdep = pp_register_vraw_contract(grad(:,:,abIndex),v); % this extracts v-dependent phase
  gradStatic = pp_register_vraw_contract(grad(:,:,abIndex),ones(1,length(abIndex))); % this extracts static phase
  gradVZ = pp_linear_static_contract(gradVdep,zernab); % project v-dep to zern basis
  gradSZ = pp_linear_static_contract(gradStatic,zernab); % project static to zern basis
  gradP = [gradSZ(:); gradVZ(:)];
end

%% Gaussian aberration + static offset
function phi = pp_register_gaussian_expand(p,coefsI,zern1,unabIndex,abIndex,v,p2f)
  % This is a nested function so that we don't have to compute gaussian
  % parameters twice (i.e., again in contract)
  % Handle the registration (common to both channels)
  phiBasis = zern1;
  coefs(unabIndex,:) = coefsI;
  Z1 = p(1:2);
  coefs(abIndex,:) = coefsI + repmat(Z1(:)',size(coefsI,1),1); % adjust with the offset
  % Create the gaussian
  phiGaussian = p2f(p(3:end));
  phiBasis(:,:,end+1) = phiGaussian;
  % Load the voltage-dependent coefficients
  coefs(abIndex,end+1) = v(:);
  phi = pp_linear_static_expand(coefs,phiBasis);
end

function gradP = pp_register_gaussian_contract(grad,zern1,abIndex,v,gf2gp)
  % Do the static-shift terms
  gradZ = zeros(2,1);
  for k = abIndex
    for cindx = 1:2
      tmp = grad(:,:,k) .* zern1(:,:,cindx);
      gradZ(cindx) = gradZ(cindx) + sum(tmp(:));
    end
  end
  % Extract the voltage-dependent part
  framesz = size(grad);
  K = length(abIndex);
  gradV = zeros(framesz(1:2));
  for k = 1:K
    gradV = gradV + v(k) * grad(:,:,abIndex(k));
  end
  % Extract the gaussian parameters
  gradG = gf2gp(gradV);
  %gradG = gf2gp(gradV,p(3:end));
  % Combine them into a single output
  gradP = [0*gradZ; gradG];
end

%% unaberrated + aberrated image
function phi = pp_pair_expand(p,zern1)
  % "Pair" has been generalized to include a first unaberrated image and
  % any number of aberrated images

  unab_coefs = p(1,:); % fixed aberration
  coefs = unab_coefs;
  coefs(2,:) = p(2,:); % p2's starting guess already contains the fixed aberration
  %coefs(2,:) = p(2,:) + unab_coefs;
  %coefs = [zeros(1,size(p,2)); p(2,:)]; % setting first images zernikes to zero

  phi = pp_linear_static_expand(coefs,zern1);
end

function gradP = pp_pair_contract(grad,zern1)
  gradP = pp_linear_static_contract(grad,zern1);
  %gradP = gradP(2:end,:);
  gradP(1,:) = zeros(1,size(gradP,2));
end



%% Utility functions
function phi = pp_zern2phi(p,zernv)
  if isvector(p) && size(p,1) == 1
    warning('phi:zernike','p may be of the wrong shape');
  end
  sz = size3(zernv);
  N = sz(3);
  K = size(p,2);
  phi = zeros([sz(1:2) K]);
  for k = 1:K
    for n = 1:N
      phi(:,:,k) = phi(:,:,k) + p(n,k) * zernv(:,:,n);
    end
  end
end

function sz = size3(A)
  sz = size(A);
  sz(end+1:3) = 1;
  if (length(sz) > 3)
    error('Input was more than 3-dimensional');
  end
end
