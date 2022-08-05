function [x1,s1] = raytrace_ideal_lens(x0,s0,lens)
% RAYTRACE_IDEAL_LENS: ray propagation through well-corrected lenses
%
% This function allows you to model the focusing properties of
% well-corrected lens systems, ones which may be assumed to be a good
% approximation of an "ideal" lens.  Such lenses satisfy a condition
% called the sine condition, and generalizations of this condition allow
% for the tracing of rays even without explicit information on the design
% and construction of the lens.   This does not require the paraxial
% approximation; it is valid even for systems with large numerical aperture.
%  
% This function is based on the following paper:
%    Holy, T.E. "Ray-tracing through ideal lens systems." Some journal,
%    2008.
% This paper should be consulted for detailed information.
%
% The optic axis is taken to be along the z-axis.  We define "right" (R)
% as the direction of increasing z, and "left" (L) as the direction of
% decreasing z.
%
% Syntax:
%   [x1,s1] = raytrace_ideal_lens(x0,s0,lens)
% where
%   x0 is a 3-by-N matrix, where each column contains the initial
%     position of a single ray as an [x;y;z] triple (again, the optic axis
%     is the z-axis), one for each of N rays;
%   s0 is a 3-by-N matrix, where each column is a unit vector giving the
%     direction of propagation at the corresponding x0 (note that
%     propagation can be either to the right or to the left, depending on
%     the sign of s0(3,i), and that this is used to define whether the
%     ray is entering the optical system from the left or from the right);
%
%   lens is a structure which may assume one of two forms.  Both forms
%   have the following information:
%     nL,nR: the refractive indices of the left and right spaces.  Note
%       that in cases where you are tracing rays with different
%       wavelengths, a chromatically-corrected lens can be treated as if
%       the refractive indices and lens properties are constant for all
%       wavelengths, so these are scalar quantities.
%     zL,zR: the z-coordinate(s) of the left (zL) and right (zR) nodal
%       points.  (In rough terms, these specify the position of the left
%       and right faces of the lens system along the optic axis; more
%       precisely, the nodal points are typically inside the lens.) The
%       position of the nodes typically varies with ray angle, but here we
%       specify only the position for the limit theta->0 (i.e., in the limit
%       of the axial ray).
%     fL,fR: the displacement from the axial nodal point (theta->0) to the
%       axial point on the corresponding focal surface. (Note sign
%       convention: for a convergent (imaging) lens, fL < 0 and fR > 0.)
%       For an infinite-conjugate system, one or the other of these (but
%       not both) can be +/-Inf.
%
%    For the "simple" form, these are all the fields you need.  (The simple
%    form is appropriate if you are interested only in an abstract model of
%    an ideal lens.)  For the "general" form (appropriate if you are
%    modeling a particular lens for which certain performance
%    specifications are available), several additional fields are
%    needed, which are provided by PREPARE_IDEAL_LENS.  See that function
%    for help.  (You can also supply them by hand, if you prefer.)
%
%    Additional optional fields:
%      start_on_focal_surface (default false): if true, the ray-tracing
%        algorithm will assume the starting positions x0 lie on the object
%        focal surface. If false, this function will first trace the rays
%        to their intersection with the focal surface.
%
% On output,
%   x1 is 3-by-N matrix, where each column specifies a position on the
%     output ray (x1 will be on the output focal surface, or, for an
%     infinite-conjugate system, on the nodal plane);
%   s1 is a 3-by-N matrix of ray propagation directions, each column for
%     a single ray.
%
% By default, the lens will be tested to see if it follows typical sign
% conventions, and a warning issued if not.  You can suppress these
% warnings using
%    warning('off','ideallens:checksigns')
%
% See also: PREPARE_IDEAL_LENS.
  
% Copyright 2008 by Timothy E. Holy
  
  %% Input validation
  % Check the lens structure and set the mode
  general = false;
  common_fields = {'zL','zR','nL','nR','fL','fR'};
  required_general_fields = {'z_N','z_Np','zeta','rho2theta','theta2rho','z_F_rho','theta_max','rho_max','tangent_rhocoef','tangent_zcoef'};
  other_general_fields = {'z_F_rhop'};
  all_general_fields = [required_general_fields other_general_fields];
  if (nargin < 3 || ~isstruct(lens))
    error('You must supply a lens structure');
  end
  fields = fieldnames(lens);
  if (length(intersect(fields,common_fields)) < length(common_fields))
    error(['You must supply all of the common fields: ' ...
      sprintf('%s, ', common_fields{1:end-1}) ...
      sprintf('and %s.',common_fields{end})]);
  end
  if (length(intersect(fields,required_general_fields)) == length(required_general_fields))
    general = true;
  elseif (length(intersect(fields,common_fields)) == length(common_fields))
    if (length(intersect(fields,all_general_fields)) > 1)
      error('Some "general" mode fields are present, but not enough to use. Quitting out of sheer confusion.');
    end
  else
    error(['You must supply all of the "simple" fields, and if desired' ...
      ' the "general" fields']);
  end
  if ~isfield(lens,'start_on_focal_surface')
    lens.start_on_focal_surface = false;
  end

  % Check x0 and s0
  [d,N] = size(x0);
  if (d ~= 3 || ~isequal(size(x0),size(s0)))
    error('You must supply 3-by-N matrices for both x0 and s0');
  end

  s0norm = sum(s0.^2,1);
  if any(abs(s0norm - 1) > 100*eps)
    warning('ideallens:normalization',...
      's0 did not consist of unit vectors. Proceeding anyway after normalizing');
    s0 = s0 ./ repmat(sqrt(s0norm),3,1);
  end

  % Basic sign convention checking (more later) & parameter extraction
  z = [lens.zL lens.zR];
  if (lens.zL > lens.zR)
    warning('ideallens:checksigns','zL and zR may be reversed');
  end
  n = [lens.nL lens.nR];
  f = [lens.fL lens.fR];
  if (f(1) > 0 || f(2) < 0)
    warning('ideallens:checksigns',...
      ['The focal lengths violate the sign convention for forming' ...
      ' a real image, is this intended?']);
  end
  if all(isinf(f))
    error('fL and fR are both infinite (not implemented)');
  end

  %% For each ray, determine whether it is entering from the left or right
  entry_sign = sign(s0(3,:));
  if any(entry_sign == 0)
    error(['All rays must propagate either to the right or left, but at' ...
      ' least one ray has s0(3) == 0']);
  end
  % For simplicity, insist that all rays are propagating in one
  % direction. (Call this function twice if you need to have both directions.)
  if any(entry_sign ~= entry_sign(1))
    error('All rays must propagate either to the right or to the left (not both)');
  end
  entry_sign = entry_sign(1);
  exit_index = (entry_sign + 3)/2;  % 1 (left) or 2 (right)
  entry_index = 3-exit_index;        % for the exiting (emitted) ray
  order_index = [entry_index exit_index];

  zero1 = zeros(1,N);
  zero2 = zeros(2,N);
  one1 = ones(1,N);
  %% Perform calculations at each surface of the lens; many of these are
  %% common, so can be done in a loop
  for i = 1:2                % i=1 is entry face, i=2 is exit face
    index = order_index(i);  % index=1 is left face, index=2 is right face
    if ~isinf(f(index))
      %% Finite conjugate
      if (i == 1)
        % We're doing the input side, trace to focal surface
        if ~lens.start_on_focal_surface
          if (general && ~isempty(lens.z_F_rho{index}))
            t = trace_to_surface(x0,s0,z(index)+f(index),index);
          else
            t = (z(index)+f(index) - x0(3,:)) ./ s0(3,:);  % to the focal plane
          end
          x0 = x0 + repmat(t,3,1).*s0;
        end
        % Calculate the radius from the position
        rho = sqrt(sum(x0(1:2,:).^2,1));
        % Calculate theta from rho
        theta = lens.rho2theta{index}(rho);
        sintheta = sin(theta);
        sinthetafix = sintheta;
        thz = (theta == 0);
        sinthetafix(thz) = 1;  % avoid some NaNs
      else
        % Calculate rho from theta
        rho = lens.theta2rho{index}(theta);
      end
      % Calculate displacements between nodal point & focal surface point
      zN = lens.z_N{index}(theta);
      zeta = lens.zeta{index}(theta);
      xN = [zero2; z(index) + zN];  % location of nodal point
      d = f(index)*(1+zeta) - zN; % displacement along z to point on focal surface
      ell(i,:) = d ./ cos(theta); % (signed) length to point on focal surface
      if (i == 1)
        % Calculate the nodal ray
        sN = (x0 - xN) ./ repmat(ell(i,:),3,1);
        % Calculate the cylindrical coordinates unit vectors
        rho_hat = sN(1:2,:) ./ repmat(-sinthetafix,2,1);
        %rho_hat = x0(1:2,:) ./ repmat(rho.*(-sign(d)),2,1);
        if any(thz)
          % theta=0 requires special treatment
          s0rs = s0(1:2,thz);
          s0n = sqrt(sum(s0rs.^2,1));
          s0rs = s0rs ./ repmat(s0n,2,1);
          s0rs(:,s0n == 0) = repmat([1;0],1,sum(s0n == 0)); % arbitrary
          rho_hat(:,thz) = s0rs;
        end
        rho_hat3 = [rho_hat; zero1];
        ephi = [-rho_hat(2,:); rho_hat(1,:); zero1];
      else
        % Calculate the output position
        x1 = repmat(ell(i,:),3,1) .* sN + xN;
      end
      % Calculate the tangent vectors
      rhocoef = lens.tangent_rhocoef{index}(theta);
      zcoef = lens.tangent_zcoef{index}(theta);
      t1 = repmat(rhocoef,3,1) .* rho_hat3 + ...
        [zero2; zcoef];
      t2 = - ephi;
      if (i == 1)
        % Calculate the dot product with the tangent vectors
        ds = s0 - sN;
        proj = repmat(n(index)*ell(i,:),2,1) .* [sum(ds .* t1,1); sum(ds .* t2,1)];
      else
        % Determine the output ray direction from the projections onto
        % the tangent vectors
        proj0 = proj;
        proj = proj ./ repmat(n(index)*ell(i,:),2,1);
        proj = proj + [sum(sN.*t1,1); sum(sN.*t2,1)];
        t1len = sqrt(sum(t1.^2,1));
        t1hat = t1 ./ repmat(t1len,3,1);
        t2hat = - ephi;
        % Projections onto unit vectors
        projhat = [proj(1,:)./t1len; proj(2,:)];
        t1perp = repmat(zcoef./t1len,3,1).*rho_hat3 + ...
          [zero2; -rhocoef./t1len];
        % Reassemble from unit vectors
        s1 = repmat(projhat(1,:),3,1) .* t1hat + ...
          repmat(projhat(2,:),3,1) .* t2hat + ...
          repmat(entry_sign*sqrt(1-sum(projhat.^2,1)),3,1) .* t1perp;
      end
    else
      %% Infinite conjugate
      if (i == 1)
        sN = s0;
        theta = acos(abs(sN(3,:)));
        sintheta = sin(theta);
        sinthetafix = sintheta;
        thz = (theta == 0);
        sinthetafix(thz) = 1;
        rho_hat = s0(1:2,:) ./ repmat(-sinthetafix,2,1);
        %rho_hat(:,thz) = repmat([1;0],1,sum(thz));  % pick something arbitrary
        if any(thz)
          % theta=0 requires special treatment
          x0rs = x0(1:2,thz);
          x0n = sqrt(sum(x0rs.^2,1));
          x0rs = x0rs ./ repmat(x0n,2,1);
          x0rs(:,x0n == 0) = repmat([1;0],1,sum(x0n == 0)); % arbitrary
          rho_hat(:,thz) = x0rs;
        end
        rho_hat3 = [rho_hat; zero1];
        ephi = [-rho_hat(2,:); rho_hat(1,:); zero1];
      else
        s1 = sN;
      end
      zN = lens.z_N{index}(theta);
      xN = [zero2; z(index) + zN];  % location of nodal point
      if (i == 1)
        % Trace to plane containing nodal point, perp to sN
        dx = xN - x0;
        t = sum(dx .* sN,1);
        dx0 = -dx + repmat(t,3,1) .* sN; % x_traced - xN (y-xN)
      end
      if (i == 1)
        % Calculate the projection with the tangent vectors
        proj = n(index) * [sum(dx0.*rho_hat3) .* sec(theta); ...
          sum(dx0.*ephi)];
      else
        % Determine the output position in the nodal plane
        proj = proj ./ (n(index)* [sec(theta); one1]);
        x1 = repmat(proj(1,:),2,1).*rho_hat + ...
          repmat(proj(2,:),2,1).*ephi(1:2,:);
        x1(3,:) = - sum(x1 .* sN(1:2,:),1) ./ sN(3,:) + xN(3,:);
      end
    end
  end
  % To compensate for roundoff errors, we normalize the length of s1
  s1 = s1 ./ repmat(sqrt(sum(s1.^2,1)),3,1);
  return;


  
  
  
  %% Do the raytracing. We handle the "simple" and "general" cases together
  %% (for the most part).  This causes a slight performance penalty on
  %% simple cases, but from a debugging perspective having common code
  %% is well worth it.
  if ~isinf(f(1))
    %% Rays are entering from a finite-conjugate side
    % Trace rays to the focal surface
    if ~lens.start_on_focal_surface
      if (general && ~isempty(lens.z_F_rho{1}))
        t = trace_to_surface(x0,s0,z(1)+f(1));
      else
        t = (z(1)+f(1) - x0(3,:)) ./ s0(3,:);  % to the focal plane
      end
    end
    x0 = x0 + repmat(t,3,1).*s0;
    % Calculate theta, d0, and ell0 for each object point
    r0 = sqrt(sum(x0(1:2,:).^2,1)); % distance from optic axis
    if general
      theta = lens.r2theta{1}(r0);
      d0 = f(1) - lens.z_N{1}(theta) + lens.zeta{1}(theta);
      ell0 = d0 ./ cos(theta);
    else
      theta = atan(r0 ./ abs(f(1)));
      ell0 = repmat(f(1),size(theta));
      %d0 = ell0 .* cos(theta);
    end
    % Calculate the perpendicular component of the nodal ray for each
    % object point
    if general
      sN = [zero2; z(1) + lens.z_N{1}(theta)] - x0;
      sN = sN ./ repmat(abs(ell0),3,1);
      % Surface normal
      g = [repmat(lens.z_F_rhop{1}(r0) ./ r0,2,1) .* x0(1:2,:); -ones(1,N)];
      normal0 = g ./ sqrt(sum(g.^2,1));
      % Get the component of the nodal ray that is perpendicular to the
      % focal surface normal
      sNperp = sN - repmat(sum(sN .* normal0,1),3,1) .* normal0;
    else
      sNperp = x0(1:2,:) ./ repmat(ell0*entry_sign,2,1);
      sNz = sqrt(1-sum(sNperp.^2,1)) * entry_sign;
      %sNperp(3,:) = 0;
    end
    % Prepare for later calculations
    F0 = ell0;
    if general
      F0 = F0 + (-lens.z_Np{1}(theta) + lens.zetap{1}(theta)) .* sin(theta);
    end
  else
    %% Rays are entering from an infinite-conjugate side
    % Calculate the nodal ray
    sNperp = s0(1:2,:);  % even for the general case, perp is wrt z-axis
    % Calculate theta
    theta = asin(sqrt(sum(sNperp.^2,1)));
    % Trace source ray to (tilted) nodal plane
    t = z(1) * s0(3,:) - sum(x0 .* s0,1);
    if general
      t = t + lens.z_N{1}(theta) .* s0(3,:);
    end
    x0 = x0 + repmat(t,3,1) .* s0;
  end
  if general
    if any(theta > lens.theta_max)
      warning('ideallens:extrapolation','Warning: theta exceeded the supplied range, extrapolation required (dangerous)');
    end
    if any(r0 > lens.r_max(1))
      warning('ideallens:extrapolation','Warning: r0 exceeded the supplied range, extrapolation required (dangerous)');
    end
  end
  if ~isinf(f(2))
    %% The output space is at finite conjugate
    % Calculate image point position
    if general
      zF1 = lens.zeta{2}(theta);
      zN1 = lens.z_N{2}(theta);
      d1 = f(2) - zN1 + zF1;
      ell1 = d1./cos(theta);
    else
      zF1 = zeros(size(theta));
      ell1 = repmat(f(2),size(theta));
      %d1 = ell1 .* cos(theta);
    end
    x1perp = repmat(ell1*entry_sign,2,1) .* sNperp;
    if general
      % Check that we are not outside of known territory
      r1 = sqrt(sum(x1perp.^2,1));
      if any(r1 > lens.r_max(2))
        warning('ideallens:extrapolation','Warning: r1 exceeded the supplied range, extrapolation required (dangerous)');
      end
    end
    x1 = [x1perp; z(2) + f(2) + zF1];
    % Calculate the output ray direction
    if general
      F1 = ell1 + (-lens.z_Np{2}(theta) + lens.zetap{2}(theta)).*sin(theta ...
	);
    else
      F1 = ell1;
    end
    if ~isinf(f(1))
      Mprime = F1*n(2)./(F0*n(1));
      s1perp = (s0(1:2,:) - sNperp)./repmat(Mprime,2,1) + sNperp;  % FIXME
    else
      coef = (n(2)/n(1)*entry_sign)*F1;
      s1perp = sNperp - x0(1:2,:)./repmat(coef,2,1);
    end
  else
    % The output space is at infinite conjugate.  Calculate the
    % position in the nodal plane.
    coef = (n(1)/n(2)*entry_sign) * F0;
    x1perp = repmat(coef,2,1) .* (sNperp - s0(1:2,:));
    x1 = [x1perp; z(2) - sum(x1perp.*sNperp,1) ./ sNz];
    % Calculate the output ray direction
    s1perp = sNperp;
  end
  s1z = sqrt(1-sum(s1perp.^2,1))*entry_sign;
  s1 = [s1perp; s1z];

  function t = trace_to_surface(x,s,z,index)
    tol = 1e-6;
    % Check to see if the current position is on the surface
    rho = sqrt(sum(x(1:2,:).^2,1));
    zdiff = z + lens.z_F_rho{index}(rho) - x(3,:);
    if all(abs(zdiff) < tol)
      t = zeros(1,N);
      return
    end
    % Check whether starting position or the point in the plane containing
    % the axial point is the best starting point
    t = (z - x(3,:)) ./ s(3,:);
    xc = x+repmat(t,3,1).*s;
    rhoc = sqrt(sum(xc(1:2,:).^2,1));
    zdiffc = z + lens.z_F_rho{index}(rhoc) - xc(3,:);
    zdiffboth = [zdiff; zdiffc];
    [tmp,bestSub] = min(abs(zdiffboth),[],1);
    bestIndex = sub2ind([2 N],bestSub,1:N);
    tboth = [zeros(1,N); t];
    t = tboth(bestIndex);
    zdiff = zdiffboth(bestIndex);
    rboth = [rho; rhoc];
    rho = rboth(bestIndex);
    xc = x + repmat(t,3,1).*s;
    % Do a Newton-type optimization
    converged = false;
    iter = 0;
    itermax = 10;
    while (~converged && iter < itermax)
      % Compute dr/dt
      rnz = rho > 0;   % have to handle rho=0 separately
      J(rnz) = sum(xc(1:2,rnz) .* s(1:2,rnz),1) ./ rho(rnz);
      J(~rnz) = sqrt(sum(s(1:2,~rnz).^2,1));
      J = -J.*lens.z_F_rhop{index}(rho) + s(3,:);
      dt = zdiff ./ J;
      t = t + dt;
      xc = x + repmat(t,3,1).*s;
      rho = sqrt(sum(xc(1:2,:).^2,1));
      zdiff = z + lens.z_F_rho{index}(rho) - xc(3,:);
      converged = all(abs(dt) < tol);
      iter = iter+1;
    end
    if (iter > itermax)
      error('Propagation to surface failed');
    end
  end
end
