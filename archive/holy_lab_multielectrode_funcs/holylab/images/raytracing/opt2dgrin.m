function [t,ro] = opt2dgrin(params,r)
% OPT2DGRIN: define a grin lens for ray tracing
% Syntax:
%   hline = opt2dgrin(params)
%   [t,ro] = opt2dgrin(params,ray)
% The first syntax plots the grin on the current axis. The second
% traces a ray.
%
% params is the parameters structure describing the planar surface:
%   xc is the location of the center of the first face
%   orient is the optic axis, a vector pointing from the center of the
%     first face to the center of the second face (the length is equal to
%     the length of the GRIN lens)
%   diam is the diameter of the lens
%   n0 is index of refraction along the axis (if available, customize for
%     your wavelength)
%   g is the gradient constant (sometimes called sqrt(A))
%   R1 and R2 are the radius of curvatures of the first and second faces,
%     respectively (use Inf if they are plano) NOTE: currently R2 must be
%     Inf.
%   mat1 is the immersion medium on the first face (e.g., 'saline')
%   mat2 is the immersion medium on the second face (e.g., 'air')
% r is the input ray (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro is the output ray after intersection.
%
% See also: RAY, OPT2DLINE, OPT2DCIRCLE, OPT2DASPHERIC.

% Copyright 2006 by Timothy E. Holy
  
% Debugging note: if you set g tiny, and make R1 finite, you get a
% plano-convex lens. That makes it easy to compare its output to pcx.

  if ~isinf(params.R2)
    error('params.R2 must be inf, currently')
  end
  grinlength = norm(params.orient);
  rotmtrx = [params.orient(1) -params.orient(2); params.orient(2) params.orient(1)]/grinlength;
  if (nargin == 1)
    % Plot GRIN outline
    y = linspace(-params.diam/2,params.diam/2,100);
    x1 = grinx(y,params.R1);
    x2 = grinlength+grinx(y,params.R2);
    L1 = rotmtrx*[x1; y] + repmat(params.xc',1,length(y));
    L2 = rotmtrx*[x2; y] + repmat(params.xc',1,length(y));
    t(1) = line(L1(1,:),L1(2,:),'Color','k');
    t(2) = line(L2(1,:),L2(2,:),'Color','k');
  else
    n_rays = length(r);
    ro = r;
    validFlag = true(1,n_rays);
    rot_e = rotmtrx'*[r.e];  % rotated back to horizontal coordinate frame
    dx = [r.x0] - repmat(params.xc,1,n_rays);
    if ~isinf(params.R1)
      b = sum(dx.*rot_e) - params.R1*rot_e(1);
      c = sum(dx.^2) - 2*params.R1*dx(1);
      sqrtarg = b.^2-c;
      validFlag = sqrtarg >= 0;
      t = -b - sqrt(sqrtarg);
      if (t < 0)
        t = -b + sqrt(sqrtarg);
      end
      if ~all(validFlag)
        % No intersection
        t(~validFlag) = NaN;
        [ro(~validFlag).valid] = false;
      end
    else
      % R1 is infinite, so it's just the intersection with the plane
      t = abs(dx(1,:))./rot_e(1,:);
    end
    % Find out if the ray strikes outside the diameter of the lens
    validFlag = validFlag | (abs(dx(2,:) + t.*rot_e(2,:)) < params.diam/2);
    if ~all(validFlag)
      t(~validFlag) = NaN;
      [ro(~validFlag).valid] = false;
      if ~any(validFlag)
        return
      end
    end
    % OK, we know it hit. Calculate the position & normal, relative to xc
    xp = dx + repmat(t,2,1).*rot_e;
    % Calculate the index of refraction at that point
    ngrini = params.n0*(1-(xp(2,:)*params.g).^2/2);
    % Calculate the immersion indices for each unique wavelength
    w = [r.w];
    [uw,tmp,wIndex] = unique(w);
    n1 = opt_refrindx(params.mat1,uw);
    n2 = opt_refrindx(params.mat2,uw);
    % Use Snell's law to calculate the angle of refraction
    if isinf(params.R1)
      normal_all = repmat([-1;0],1,n_rays);
    else
      normal_all = repmat([params.R1;0],1,n_rays) - xp;
      normal_all = -normal_all./repmat(norm(normal_all),2,1);
    end
    for rayIndex = 1:n_rays
      if ro(rayIndex).valid
        normal = normal_all(:,rayIndex);
        e = rot_e(:,rayIndex);
        P = perpproj(normal);
        nratio = n1(wIndex(rayIndex))/ngrini(rayIndex);
        e2 = nratio*P*e;  % Snell's law
        e2 = e2 + sign(sum(e .* normal)) * sqrt(1-sum(e2.*e2))*normal;
        % Calculate the parameters of the sinsoid A sin(gx) + B cos(gx)
        gx = params.g * xp(1,rayIndex);
        AB = inv([sin(gx) cos(gx); cos(gx) -sin(gx)]) * [xp(2,rayIndex); e2(2)/(params.g*e2(1))];
        % Will this vector propagate? See if the peak of the sinusoid is
        % outside the diameter of the GRIN
        gxexit = params.g * grinlength;
        if (all(AB) == 0)
          zmax = 0;
        else
          zmax = atan(AB(1)/AB(2));
        end
        if (zmax < 0)
          zmax = zmax+pi;
        end
        zmax = min(zmax,gxexit);
        ymax = AB(1)*sin(zmax) + AB(2)*cos(zmax);
        if (abs(ymax) > params.diam/2)
          t(rayIndex) = NaN;
          ro(rayIndex).valid = false;
        end
        % Now propagate this vector until it encounters the second surface
        % fixme: this is where requiring that R2=inf comes in handy
        yexit = AB(1)*sin(gxexit) + AB(2)*cos(gxexit);
%     xstart = xp + params.xc;
%     xend = [grinlength,yexit] + params.xc;
%     line([xstart(1) xend(1)],[xstart(2) xend(2)],'Color','m')
        dydxexit = params.g * (AB(1)*cos(gxexit) - AB(2)*sin(gxexit));
        exexit = 1/sqrt(1+dydxexit^2);
        eexit = [exexit; dydxexit * exexit];
        % Finally, Snell's law again
        xpo = [grinlength; yexit];
        ngrino = params.n0*(1-(xpo(2)*params.g)^2/2);
        normal = [-1;0];
        P = perpproj(normal);
        nratio = ngrino/n2(wIndex(rayIndex));
        e2 = nratio*P*eexit(:);  % Snell's law
        if (norm(e2) <= 1)
          e2 = e2 + sign(sum(eexit .* normal)) * sqrt(1-sum(e2.*e2))*normal;
        else
          % Total internal reflection off the back surface
          e2 = eexit - 2*sum(eexit.*normal)*normal;
        end
        % Output ray---remember to rotate back
        ro(rayIndex).x0 = (rotmtrx*xpo) + params.xc;
        ro(rayIndex).e = rotmtrx*e2;
      end
    end
  end
  
function x = grinx(y,R)
  if ~isinf(R)
    x = abs(R)-sqrt(R^2-y.^2);
    x = x*sign(R);
  else
    x = 0*y;
  end
