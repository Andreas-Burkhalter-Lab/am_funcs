function lens = prepare_ideal_lens(f,nodes,FL,FR)
% PREPARE_IDEAL_LENS: configure an ideal lens for raytracing
% This function creates the spline interpolants needed for raytracing by
% RAYTRACE_IDEAL_LENS.
%
% This function is based on the following paper:
%    Holy, T.E. "Ray-tracing through ideal lens systems." Some journal,
%    2008.
% It should be consulted for detailed information.
%
% Syntax:
%   lens = prepare_ideal_lens(f,nodes,FL,FR)
% where
%   f is a 2-vector [fL fR] containing the displacement from the axial
%     nodal point to the axial point on the corresponding focal surface.
%     fL is for the left node/focal surface, fR for the right node/focal
%     surface.  Note that the sign convention is that rightward is
%     positive, i.e., if the focal surface is to the left of its node,
%     then f < 0; if the surface is to the right, f > 0.  An infinite
%     conjugate is described as f = +/-inf.
%   nodes is a 3-by-k matrix describing the nodal behavior of the lens.
%     The first row is theta, which should be arranged in increasing order
%     and should start from 0.  The second row is the displacement of the
%     left node for each value of theta, following the same sign convention
%     as described above for f.  The displacement should be 0 for theta=0.
%     The third row is the displacement of the right node as a function
%     of theta.
%   FL and FR are 2-by-m matrices that describe the shape of the focal
%     surfaces (i.e., the field curvature).  The first row is rho, the
%     radial distance from the optic axis.  The second row is the
%     displacement of the focal surface along the optic axis as a
%     function of rho, measured from the plane containing the axial point
%     of the focal surface.  This too obeys the sign convention for f.  rho
%     should start at 0, and the corresponding displacement should also
%     be 0.
%     FL and FR do not have to be of the same size.  If the lens has a
%     flat field on one or both of the sides, you can supply the empty
%     matrix.  You can also leave it empty for a side operating at
%     infinite conjugate.
%
% On output, the structure "lens" will contain fields that are 1-by-2
% cell arrays of function handles.  In each such cell array, the first is
% for the left side, the second for the right side.  These fields are the
% following: 
%   z_N: the nodal positions z_N(theta)
%   z_F_rho: the focal surface displacement as a function of rho, z_F(rho)
%   zeta: the focal surface displacement as a function of theta, z_F(theta)
%   r2theta: a function converting radius to theta
%   theta2r: a function converting theta to radius
% Other fields are:
%   fL and fR (copied from f, see above)
%   theta_max: the maximum theta for which z_N can be calculated without
%     extrapolation;
%   rho_max: a 1-by-2-vector, [rhoLmax rhoRmax] containing the maximum rho
%     forwhich z_F_rho can be calculated without extrapolation.
%
% Before calling RAYTRACE_IDEAL_LENS, you must also set the nL, nR, zL,
% and zR fields of lens. (See RAYTRACE_IDEAL_LENS.)
%
% See also: RAYTRACE_IDEAL_LENS.
  
% Copyright 2008 by Timothy E. Holy
    
  lens.fL = f(1);
  lens.fR = f(2);
  %% Process the node information
  default_nodes = zeros(3,4); default_nodes(1,:) = linspace(0,pi/2,4);
  if isempty(nodes)
    nodes = default_nodes;
  end
  [sz,n_theta] = size(nodes);
  if (sz ~= 3)
    error('nodes must be 3-by-k');
  end
  lens.theta_max = nodes(1,end);
  if any(nodes(:,1) ~= 0)
    error('The nodal displacement for theta=0 should be 0');
  end
  if (n_theta < 4 && n_theta ~= 3)
    % There aren't enough points specified for spline interpolation, fill
    % in intermediate values
    warning('ideallens:interpolation',...
	    ['Not enough theta values specified to create splines, assuming' ...
	     ' constant or linear behavior']);
    if (n_theta == 1)
      % Must be constant & 0
      nodes = default_nodes;
    elseif (n_theta == 2)
      % Assume linear behavior
      nodes = [nodes(:,1) nodes(:,2)*([1 2 3]/3)];
    end
  elseif (n_theta == 3)
    error(['Not enough theta values specified to create splines, but with 3' ...
	   ' provided values linear behavior can''t be assumed: provide a 4th point']);
  end
  for index = 1:2
    lens.z_N{index} = ppcreate(nodes(1,:),nodes(index+1,:),'spline');
    lens.z_Np{index} = ppcreate(lens.z_N{index},'diff');
  end

  %% Process the information about field curvature
  F = {FL,FR};
  r_max = nan(1,2);
  for index = 1:2
    if ~isempty(F{index})
      r_max(index) = F{index}(1,end);
      lens.z_F_rho{index} = ppcreate(F{index}(1,:),F{index}(2,:),'spline');
      lens.z_F_rhop{index} = ppcreate(lens.z_F_rho{index},'diff');
    else
      lens.z_F_rho{index} = [];
      % If field curvature data are not supplied, we set it to have a
      % flat field. Do it as a spline so 'diff' works
      %rmax = tan(lens.theta_range(2)) * (f(index) + nodes(1+index,end));
      %lens.z_F_rho{index} = ppcreate(linspace(0,rmax,4),[0 0 0 0]);
    end
  end
  
  %% Set up functions for converting r0 <--> theta <--> r1
  % Also calculate the field curvature as a function of theta
  for index = 1:2
    if isinf(f(index))
      lens.rho2theta{index} = [];
      lens.theta2rho{index} = [];
      lens.zeta{index} = [];
      rho_max(index) = inf;
      continue
    end
    if isempty(F{index})
      % Has a flat field, the relationship can be determined directly
      theta = nodes(1,:);
      z_N = nodes(1+index,:);
      rho = tan(theta) .* abs(f(index) - z_N);
      z_F = zeros(size(theta));
      rho_max(index) = rho(end);
    else
      % Has a curved field, we have to iteratively solve the nonlinear
      % equation for theta given rho
      rho = F{index}(1,:);
      rho_max(index) = max(rho);
      z_F = F{index}(2,:);
      theta = atan(rho/abs(f(index)));  % starting guess
      d = f(index) - lens.z_N{index}(theta) + z_F;
      s = sign(d);
      converged = false;
      iter = 0;
      itermax = 20;
      % Iterate by Newton's method. We get global convergence by insuring
      % that each accepted step is no worse than the previous position.
      resid = rho./s - tan(theta).*d;
      while (~converged && iter < itermax)
        J = sec(theta).^2 .* d - tan(theta) .* lens.z_Np{index}(theta);
        dtheta = resid ./ J;
        converged = all(abs(dtheta) < 1e-6);
        iter = iter+1;
        thetanew = theta+dtheta;
        % Insure that theta doesn't go out of bounds
        thetanew(thetanew > pi/2) = pi/2;
        dtheta = thetanew-theta;
        d = f(index) - lens.z_N{index}(thetanew) + z_F;
        % Insure global convergence
        residnew = rho./s - tan(thetanew).*d;
        isworse = abs(residnew) > abs(resid);
        while any(isworse)
          dtheta(isworse) = dtheta(isworse)/2;
          thetanew = theta+dtheta;
          thetanew(thetanew > pi/2) = pi/2;
          d = f(index) - lens.z_N{index}(thetanew) + z_F;
          residnew = rho./s - tan(thetanew).*d;
          isworse = abs(residnew) > abs(resid);
        end
        theta = thetanew;
        resid = residnew;
      end
      if (iter > itermax)
        error('Calculation of theta from rho failed!');
      end
      if any(theta > lens.theta_max)
        warning('ideallens:extrapolation',...
          ['Insufficient theta range supplied for the nodal behavior;' ...
          ' extrapolation was required (dangerous)']);
      end
    end
    lens.rho2theta{index} = ppcreate(rho,theta,'spline');
    lens.theta2rho{index} = ppcreate(theta,rho,'spline');
    lens.zeta{index} = ppcreate(theta,z_F/f(index),'spline');
  end

  lens.rho_max = rho_max;

  %% Set up components of tangent vectors needed for raytracing
  for index = 1:2
    if ~isinf(f(index))
      zeta = lens.zeta{index}(theta);
      zetap_fcn = ppcreate(lens.zeta{index},'diff');
      zetap = zetap_fcn(theta);
      d = f(index)*(1+zeta) - lens.z_N{index}(theta);
      dp = -lens.z_Np{index}(theta) + f(index)*zetap;
      lens.tangent_rhocoef{index} = ppcreate(theta,...
        -(dp./d .* sin(theta) + sec(theta)),'spline');
      lens.tangent_zcoef{index} = ppcreate(theta,...
        f(index)*zetap .* cos(theta)./ d, 'spline');
    end
  end
