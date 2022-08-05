function [e0,e1,x0tildep,x1tildep] = points2ray_sinelens(dx0,dx1,f)
% POINTS2RAY_SINELENS: solve for ray parameters given initial & final points
%
% This function finds the ray that passes through two points, an object
% point and an image point, for a lens that satisfies the sine condition.
% The results are reported in terms of the propagation directions in the
% object and image field. The lens's optic axis is oriented along the
% z-axis, where coordinates are parametrized as [x y z].
%
% Syntax:
%   [e0,e1] = points2ray_sinelens(dx0,dx1,f)
% where
%   dx0 is the separation of the object point, x0, from the axial point xf
%     in the focal plane, i.e., dx0 = x0-xf;
%   dx1 is the separation of the final point, x1, from the axial point xl in
%     the objective lens plane, i.e., dx1 = x1 - xl (i.e. xl is the
%     center of the back aperture of the lens);
%   f is a two vector giving the [front back] focal lengths of the
%     objective.  For a convergent lens, f(1) < 0 and f(2) > 0. 
% and
%   e0 is the 3-vector of the propagation direction before the objective
%     (NOTE: will be NaNs if no such ray exists!);
%   e1 is the 3-vector of the propagation direction after the objective.
%
% Two additional outputs are also available:
%   [e0,e1,x0tildep,x1tildep] = points2ray_sinelens(dx0,dx1,f)
% where
%   x0tildep is the two-vector of the ray's strike point in the object
%     focal plane (i.e., relative to xf);
%   x1tildep is the two vector of the ray's strike point at the exit of the
%     objective (i.e., relative to xl);
%
% The sine lens is modeled according to de Meijere & Velzel (1987),
% "Linear ray-propagation models in geometrical optics," J. Opt. Soc. Am.
% A, Vol. 4, 2162-2165. The front focal plane and back aperture were picked
% as the two planes to be imaged in a manner consistent with the sine
% condition (and rays from points not on these planes are traced by first
% propagating to these two planes). The sign convention for the lens
% follows their sign convention, which is opposite that of Born & Wolf.
%
% See also my (unpublished) writeup on the topic.

% Copyright 2008 by Timothy E. Holy

  dx0 = dx0(:);
  dx1 = dx1(:);
  if (dx1(3) == 0)
    % Special code path for what should be the most important case: tracing
    % to the back aperture of the lens. This version has been quite
    % carefully tuned, and should be fast and robust.
    x1tildep = dx1(1:2);
    if (dx0(3) == 0)
      % Trivial case (analytically solvable)
      x0tildep = dx0(1:2);
      e0tilde = (x0tildep-x1tildep)/f(1);
      e1tilde = -x0tildep/f(2);
      e0t2 = e0tilde'*e0tilde;
      e1t2 = e1tilde'*e1tilde;
      if (e0t2 > 1 || e1t2 > 1)
        e0 = nan(3,1);
        e1 = nan(3,1);
      else
        e03 = sqrt(1-e0t2);
        e0 = [e0tilde; e03];
        e13 = sqrt(1-e1t2);
        e1 = [e1tilde; e13];
      end
      %p2rsl_finish_e;
      return;
    end
    dx = x1tildep-dx0(1:2);
    gamma = sqrt((dx'*dx))/abs(f(1));
    delta = dx0(3)/f(1);  % defocus compared with focal length of lens
    if (delta < 0)
      no_solution = false;
      if (gamma > 1 || delta <= -1)
        no_solution = true;  % easy case!
      else
        umax = sqrt(1-(-delta)^(2/3));  % value of u at the maximum
        fmax = umax + delta*umax/sqrt(1-umax^2); % value of function at max
        if (fmax < gamma)
          no_solution = true;
        end
      end
      if no_solution
        e0 = nan(3,1);
        e1 = nan(3,1);
        x0tildep = dx0(1:2);
        return;
      end
      % Pick an initial guess by linear interpolation from 0 to the maximum
      % (we want the leftmost root)
      u = umax * (gamma/fmax);
    else
      umax = 1;
      % Initial guess is made by a 3-part fit
      if (gamma < 1)
        u = gamma/(1+delta);
      elseif (gamma < 1 + delta/sqrt((1+delta)^2 - 1))
        u = 1/(1+delta);
      else
        u = (gamma-1)/sqrt((gamma-1)^2 + delta^2);
      end
    end
    % Now polish this guess by Newton's method
    du = Inf;
    du_thresh = 1e-8*gamma;
    while (abs(du) > du_thresh)
      tmp = sqrt(1-u^2);
      r = u + delta*u/tmp - gamma;  % residual
      J = 1 + delta/tmp^3; % Jacobian
      du = -r/J;
      % Keep u between 0 and umax
      unew = u+du;
      if (unew < 0)
        u = u/2;
      elseif (unew >= umax)
        u = (u+umax)/2;
      else
        u = unew;
      end
    end
    alpha = u/gamma; % if delta < 0, we expect alpha > 1, otherwise 0<alpha<=1
    x0tildep = alpha*dx0(1:2) + (1-alpha)*x1tildep;
    e0tilde = (x0tildep-x1tildep)/f(1);
    e1tilde = -x0tildep/f(2);
    e0t2 = e0tilde'*e0tilde;
    e1t2 = e1tilde'*e1tilde;
    if (e0t2 > 1 || e1t2 > 1)
      e0 = nan(3,1);
      e1 = nan(3,1);
    else
      e03 = sqrt(1-e0t2);
      e0 = [e0tilde; e03];
      e13 = sqrt(1-e1t2);
      e1 = [e1tilde; e13];
    end
    %p2rsl_finish_e;
    return;
  end
  % Set up initial guesses: for x0, just snip off 3rd coordinate (project along
  % z-axis), on the assumption we're pretty close to the focal plane and
  % that the result of any propagation to the focal plane will not cause
  % much of a shift
  x0tildepnew = dx0(1:2);
  % For x1, use this guess for x0 to calculate the propagation direction,
  % and then back-propagate to the back aperture of the objective
  e1tilde = -x0tildepnew/f(2);
  x1tildepnew = dx1(1:2) - (dx1(3)/sqrt(1-e1tilde'*e1tilde)) * e1tilde;
  % Initialize for an iterative solution
  dxtotal = inf;
  dxthresh = 1e-8 * sum(abs(f));
  failed = false;
  f2 = f.^2;
  one = eye(2,2);
  cheat = 0;
  cheat_thresh = 3;
  while (dxtotal > dxthresh && cheat < cheat_thresh)
    x0_2 = x0tildepnew'*x0tildepnew;
    dxtp = x1tildepnew-x0tildepnew;
    dxtp_2 = dxtp'*dxtp;
    if (dxtp_2 > f2(1))
      dxtp_2 = 0.99*f2(1);
      cheat = cheat+1;
    end
    if (x0_2 > f2(2) || dxtp_2 > f2(1))
      % Requires rays that cannot have unit norm
      failed = true;
      break;
    end
    x0tildep = x0tildepnew;
    x1tildep = x1tildepnew;
    a = dx1(3)/sqrt(1-x0_2/f2(2))/f(2);
    b = dx0(3)/sqrt(1-dxtp_2/f2(1))/f(1);
    Ma = ((a/(f2(2)-x0_2))*x0tildep)*x0tildep';
    Mb = ((b/(f2(1)-dxtp_2))*dxtp)*dxtp';
    A = [a*one + Ma, -one; ...
      -(b*one+Mb) - one, b*one+Mb];
    r = [x1tildep - dx1(1:2) - a*x0tildep;...
      x0tildep - dx0(1:2) - b*dxtp];
    dx = A\r;
    dxtotal = sum(abs(dx));
    x0tildepnew = x0tildep + dx(1:2);
    x1tildepnew = x1tildep + dx(3:4);
  end
  if (failed || cheat == cheat_thresh)
    %warning('Failed to converge');
    e0 = nan(3,1);
    e1 = nan(3,1);
%     x0tildep = nan(2,1);
%     x1tildep = nan(2,1);
  else
    x0tildep = x0tildepnew;
    x1tildep = x1tildepnew;
    e0tilde = - (x1tildep - x0tildep)/f(1);
    e1tilde = - x0tildep/f(2);
    p2rsl_finish_e;
  end
  return
  
  % Old, slow, fallback code in case the above doesn't converge
  % Pack into a single vector
  p0 = [dx0(1:2); dx1(1:2)];
  [p,fval,exitflag] = fsolve(@tis_fun,p0,struct('Display','none'));
  % They're mostly unpacked already, thanks to the nested function
  e0 = e0tilde;
  e0(3) = e03;
  e1 = e1tilde;
  e1(3) = e13;
  if (exitflag ~= 1 || e0'*e0 > 1.0001 || e1'*e1 > 1.0001)
    warning('Failed to converge');
    e0 = nan(size(e0));
    e1 = nan(size(e1));
    %x0tildep = nan(size(x0tildep));
    %x1tildep = nan(size(x1tildep));
  else
    fprintf('!');
  end
  
  function ret = tis_fun(pp)
    % This is the function that imposes the conditions
    % First, unpack the variables
    x0tildep = pp(1:2);
    x1tildep = pp(3:4);
    e0tilde = - (x1tildep - x0tildep)/f(1);
    e03 = abs(sqrt(1-e0tilde'*e0tilde));
    e1tilde = - x0tildep/f(2);
    e13 = abs(sqrt(1-e1tilde'*e1tilde));
    
    % The t's are the deviations from the desired solution
    t0 = dx0(1:2) - (dx0(3)/e03)*e0tilde - x0tildep;
    t1 = dx1(1:2) - (dx1(3)/e13)*e1tilde - x1tildep;
    ret = [t0; t1];
  end
  
  function p2rsl_finish_e
    % A nested function to finish the computation of e0, e1
    e0t2 = e0tilde'*e0tilde;
    e1t2 = e1tilde'*e1tilde;
    if (e0t2 > 1 || e1t2 > 1)
      e0 = nan(3,1);
      e1 = nan(3,1);
    else
      e03 = sqrt(1-e0t2);
      e0 = [e0tilde; e03];
      e13 = sqrt(1-e1t2);
      e1 = [e1tilde; e13];
    end
  end
end