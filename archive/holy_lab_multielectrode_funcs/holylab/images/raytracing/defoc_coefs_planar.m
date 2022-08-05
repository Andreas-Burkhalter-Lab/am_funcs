function [Cz,Cy] = defoc_coefs_planar(alpha,tm,yrange)
% DEFOC_COEFS_PLANAR: defocus coefficients for planar "tissue"
% Syntax:
%   [Cz,Cy] = defoc_coefs_planar(alpha,theta_max)
% where
%   alpha is the angle of the tissue's normal with respect to the optic
%     axis;
%   theta_max is the angle of the margin ray, i.e., asin(NA/n), where n
%     is the index of refraction of the immersion medium.
% and
%   Cz and Cy are the defocus coefficients described in Tim's notes.
  
% Copyright 2008 by Timothy E. Holy
  
  ta = tan(alpha);
  if (nargin < 3)
    % We're using the whole aperture, so it's much more efficient to do the
    % integral in polar coordinates
    sc = sin(tm)*cos(tm);  % twice the integral of sin(t)^2 from 0 to tm
    zfun0 = @(t) (sin(tm)./cos(t) - tm)./sqrt(1-ta^2*tan(t).^2);
    Cz = quad(zfun0,0,tm);
    Cz = Cz / (tm*(tm+sc)/2 - sin(tm)^2);

    yfun0 = @(t) tan(t).^2./(1 - ta^2*tan(t).^2 + sqrt(1 - ta^2*tan(t).^2));
    Cy = ta * (4*quad(yfun0,0,tm)/(tm - sc));
  else
    % We're blocking part of the aperture, so we have to do it in cartesian
    % coordinates
    a = sin(tm);
    rngs = {-1,1,yrange(1),yrange(2)};
    b = - [dblquad(@(x,y) cfun(x,y) .* srcfun(x,y), rngs{:});
      dblquad(@(x,y) zfun(x,y) .* srcfun(x,y), rngs{:});
      dblquad(@(x,y) yfun(x,y) .* srcfun(x,y), rngs{:})];
    A = [dblquad(@(x,y) cfun(x,y).^2,rngs{:}) dblquad(@(x,y) cfun(x,y).*zfun(x,y),rngs{:}) dblquad(@(x,y) cfun(x,y).*yfun(x,y),rngs{:});
      0 dblquad(@(x,y) zfun(x,y).^2,rngs{:}) dblquad(@(x,y) zfun(x,y).*yfun(x,y),rngs{:});
      0 0 dblquad(@(x,y) yfun(x,y).^2,rngs{:})];
    A(2,1) = A(1,2);
    A(3,1) = A(1,3); A(3,2) = A(2,3);
    sol = A\b;
    %C0 = sol(1);
    Cz = sol(2);
    Cy = sol(3);
  end

    function r = srcfun(x,y)
      r = zeros(size(x));
      keep = x.^2 + y.^2 <= 1;
      if ~any(keep)
        return
      end
      xkeep = x(keep);
      r(keep) = 1 ./ (sqrt(1 - a^2 * (xkeep.^2+y.^2)) + (a*ta)*y);
    end
    function r = zfun(x,y)
      r = zeros(size(x));
      keep = x.^2 + y.^2 <= 1;
      if ~any(keep)
        return
      end
      xkeep = x(keep);
      r(keep) = sqrt(1 - a^2 * (xkeep.^2+y.^2));
    end
    function r = yfun(x,y)
      r = zeros(size(x));
      keep = x.^2 + y.^2 <= 1;
      if ~any(keep)
        return
      end
      r(keep) = a*y;
    end
    function r = cfun(x,y)
      r = zeros(size(x));
      keep = x.^2 + y.^2 <= 1;
      if ~any(keep)
        return
      end
      r(keep) = -1;
    end

end
  
