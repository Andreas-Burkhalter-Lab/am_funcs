function varargout = biharmonic_penalty(param,v,pupildata,im)
% Syntax:
%   phik = biharmonic_penalty(param,v,pupildata)
%   err = biharmonic_penalty(param,v,pupildata,im)
%   [err,gparam] = biharmonic_penalty(param,im,v,pupildata)
% param has the following fields:
%   m,zeta: vectors (or scalars) of slopes and intercepts for the
%     voltage-to-displacement relation. If supplied as a vector, there
%     should be one element for each actuator.
%   a: a 3-vector or 4-vector. The 3-vector gets unpacked as
%         A = [a(1) a(2); a(2) a(3)],
%     a symmetric matrix which describes the pupil coordinate rotation and
%     scaling for its projection onto the DM. The 4-vector is unpacked as
%        reshape(a,[2 2]) = [a(1) a(3); a(2) a(4)].
%   xi0: a 2-vector describing the offset (translation) component of the
%     pupil.  In total, pupil coordinates u are related to mirror
%     coordinates xi by
%         xi = A * u + xi0
%   R: radius at which the membrane is clamped, in mirror coordinates xi
%   act_xi: mirror coordinates of the actuators, in mirror coordinates
%     xi.  This should be n_actuators-by-2.
%   alpha (default 0): coefficient of defocus
% im: ny-by-nx-by-K image "stack" of K images
% v: a K-by-n_actuators matrix, specifying the setting of all voltages
%   corresponding to each image
%
% gparam is a structure with the same fields as param, containing the
%   gradient of the penalty with respect to each coordinate (exception: no
%   act_xi field).

% Copyright 2009 by Timothy E. Holy
  
  %% Parse the inputs and do the forward computations needed for phase calculations
  inpupil = pupildata.H0>0;
  theta = pupildata.theta(inpupil);
  u = repmat(pupildata.rho(inpupil),1,2) .* [cos(theta),sin(theta)];
  n_pts = size(u,1);
  if (length(param.a) == 3)
    A = reshape(param.a([1 2 2 3]),2,2);
  else
    A = reshape(param.a,2,2);
  end
  xi = u*A + repmat(param.xi0(:)',n_pts,1);
  R = param.R;
  [K,n_actuators] = size(v);
  zeta = param.zeta;
  if isscalar(zeta)
    zeta = repmat(zeta,1,n_actuators);
  end
  bflag = any([v;zeta] ~= 0,1);
  b = zeros(n_pts,n_actuators);     % shape due to each actuator
  gub = zeros(n_pts,2,n_actuators); % gradient of b with respect to u
  gRb = zeros(n_pts,n_actuators);   % gradient of b with respect to R
  for actIndex = 1:n_actuators
    if bflag(actIndex)
      [b(:,actIndex),gub(:,:,actIndex),gRb(:,actIndex)] = biharmonic_pointsource(xi,param.act_xi(actIndex,:),R);
    end
  end
  m = param.m;
  if isscalar(m)
    m = repmat(m,1,n_actuators);
  end
  alpha = 0;
  if isfield(param,'alpha')
    alpha = param.alpha;
  end

  imsz = size(pupildata.H0);
  imsz(end+1:2) = 1;
  imsz(3) = K;
  
  %% Compute the phase associated with each image
  phik = zeros(imsz);
  Z4 = (2*pupildata.rho.^2 - 1).*pupildata.H0;
  for k = 1:K
    psi = zeros(n_pts,1);
    for actIndex = find(bflag)
      psi = psi + (m(actIndex)*v(k,actIndex)+zeta(actIndex))*b(:,actIndex);
    end
    phi = zeros(imsz(1:2));
    phi(inpupil) = psi;
    phik(:,:,k) = phi + alpha*Z4;
  end
  
  %% Check to see if we need an early return
  if (nargin < 4)
    varargout = {phik};
    return
  end
  
  if (imsz(3) ~= K)
    error('Size of im and v do not agree');
  end
  
  
  %% Calculate value and gradient with respect to phi
  varargout = cell(1,nargout);
  [varargout{:}] = pdpenalty(phik,im,pupildata.H0);
  
  %% Pack the gradient with respect to the input parameters
  if length(varargout) > 1
    gphik = varargout{2};  % save the gradient with respect to phi
    gparam = rmfield(param,'act_xi');  % use the same structure format
    % Pull out just the (nonzero) pixels that are in the pupil
    gphik_inpupil = zeros(n_pts,K);
    for k = 1:K
      gphitmp = gphik(:,:,k);
      gphik_inpupil(:,k) = gphitmp(inpupil);
    end
    % Gradient with respect to m and zeta
    gparam.m = zeros(1,n_actuators);
    gparam.zeta = zeros(1,n_actuators);
    for actIndex = find(bflag)
      for k = 1:K
        tmp = sum(b(:,actIndex).*gphik_inpupil(:,k));
        gparam.m(actIndex) = gparam.m(actIndex)+v(k,actIndex)*tmp;
        gparam.zeta(actIndex) = gparam.zeta(actIndex) + tmp;
      end
    end
    if isscalar(param.m)
      gparam.m = sum(gparam.m);
    end
    if isscalar(param.zeta)
      gparam.zeta = sum(gparam.zeta);
    end
    % Compute gradients with respect to the rigid registration parameters
    % and R
    gparam.a = zeros(1,length(param.a));
    gparam.xi0 = zeros(1,2);
    gparam.R = 0;
    for k = 1:K
      % Compute grad_u psi
      gupsi = zeros(n_pts,2);
      gRpsi = zeros(n_pts,1);
      for actIndex = find(bflag)
        tmp = m(actIndex)*v(k,actIndex) + zeta(actIndex);
        gupsi = gupsi + tmp*gub(:,:,actIndex);
        gRpsi = gRpsi + tmp*gRb(:,actIndex);
      end
      % Compute individual gradient terms
      if (length(param.a) == 3)
        % 3-parameter A
        gparam.a(1) = gparam.a(1) + sum(gupsi(:,1) .* u(:,1) .* ...
          gphik_inpupil(:,k));
        gparam.a(2) = gparam.a(2) + sum(sum(gupsi .* u(:,[2 1]),2) .* ...
          gphik_inpupil(:,k));
        gparam.a(3) = gparam.a(3) + sum(gupsi(:,2) .* u(:,2) .* ...
          gphik_inpupil(:,k));
      else
        % 4-parameter A
        thisg = zeros(2,2);
        for i1 = 1:2
          for i2 = 1:2
            thisg(i1,i2) = sum(gupsi(:,i1) .* u(:,i2) .* gphik_inpupil(:,k));
          end
        end
        gparam.a = gparam.a + reshape(thisg',size(gparam.a));
      end
      gparam.xi0 = gparam.xi0 + sum(gupsi .* repmat(gphik_inpupil(:,k),1,2),1);
      gparam.R = gparam.R + sum(gRpsi .* gphik_inpupil(:,k));
    end
    % Compute gradient with respect to alpha (defocus)
    if isfield(param,'alpha')
      gparam.alpha = 0;
      for k = 1:K
        tmp = gphik(:,:,k) .* Z4;
        gparam.alpha = gparam.alpha + sum(tmp(:));
      end
    end
    
    varargout{2} = gparam;
  end
