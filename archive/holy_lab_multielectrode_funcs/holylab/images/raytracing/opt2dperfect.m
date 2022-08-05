function [t,ro] = opt2dperfect(params,r)
% this is a hack version
% y is not equal to yout
% changed sign of thetap ... both for no logical reason other than it
% sorta worked
% OPT2DPERFECT: paraxial approximation for all rays (even non-paraxial)
% Syntax:
%   hline = opt2dprojective(params)
%   [t,ro] = opt2dprojective(params,ray)
% The first syntax plots the "surface" on the current axis. The second
% traces a ray.
%
% Params is the parameters structure describing the projective
% transformation:
%   xc is the center of the transformation (i.e., position along the optic
%     axis), a 2-vector
%   normal is a unit vector defining the optical axis, a 2-vector
%   f is a 2-vector of focal lengths, [frontfocallength backfocallength]
%   plotsize is the length of the line used to represent this element when
%     plotted on the screen (this has no functional consequences, only for
%     visual representation; default 1)
% r is the input ray (see RAY)
% and
%   t is the distance along the ray to the point of intersection (NaN if no
%     intersection);
%   ro is the output ray after intersection.
%
% See also: RAY, OPT2DCIRCLE.

  if (nargin == 1)
    % Plot surface on screen
    if ~isfield(params,'plotsize')
        params.plotsize = 1;
    end
    parallel = [params.normal(2) -params.normal(1)];
    xl = params.xc + params.plotsize*parallel;
    xr = params.xc - params.plotsize*parallel;
    t = line([xl(1) xr(1)],[xl(2) xr(2)],'Color','k','LineWidth',2);
  else
      R = opt2p_rotmtrx(params);
      dx = r.x0 - params.xc;
      dx = dx*R;   % rotate so that optic axis is along x
      e = r.e*R;
      if (e(1) <= 0) % hack...rays only going to positive x will be considered
          % It's going the wrong way, it never intersects
          t = nan;
          ro = struct;
          return
      end
      t = -dx(1)/e(1);
      % The strike position is [0 y] (in our new coordinate system)
      y = dx(2)-dx(1)*e(2)/e(1);
      % Now compute the image of the source point x0.
      if any(isinf(params.f))
          error('We don''t build no stinkin telescopes here!');
      end
      ff = prod(params.f);
      z = dx(1) + params.f(1);  % measure distance from front focal plane
      if (z ~= 0)
          zp = ff/z;
          yp = y*params.f(1)/z;
          dxp = [zp+params.f(2) yp];  % this is the image point
          % Now compute the ray direction
          tantheta = e(2)/e(1);
          tanthetap = tantheta * (-z/params.f(2));
          l = sqrt(tanthetap^2+1);
          tanthetap = -tanthetap;
          eout = [1 tanthetap]/l;
          % Back-propagate this so it seems as if it's leaving the x=0
          % line
          ttmp = -dxp(1)/eout(1);
          yout = dxp(2) + ttmp*eout(2);
          %xout = [0 yout];
          xout = [0 y];
          
      else
          eout = [sign(e(1)) 0];
          xout = [0 y];
      end
      ro = r;
      ro.x0 = xout*R' + params.xc;
      ro.e = eout*R';
  end
  
function R = opt2p_rotmtrx(params)
  R = [params.normal; -params.normal(2) params.normal(1)]; % Rotation matrix
%   R = R * [0 1; -1 0]; % To compensate for the fact that y = f(x), but
                       % theta=0 corresponds to propagation in the x
                       % direction
