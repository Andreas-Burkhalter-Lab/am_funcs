function [x,fx] = brentmin(f,xbrack,fxbrack,options)
% BRENTMIN: 1-dimensional minimization by Brent's method
%
% This differs from FMINBND in two ways:
%  1. Information about the middle point in the triple is retained. Thus,
%     it may not throw away as much useful work as FMINBND.
%  2. This function allows you to set a tolerance on function value, if you
%     choose. The default convergence criterion is on the parameter value.
%
% Syntax:
%   [x,fx] = brentmin(f,xbrack,fxbrack)
%   [x,fx] = brentmin(f,xbrack,fxbrack,options)
% where
%   f is the function to be minimized
%   xbrack is a bracketing triple of points
%   fxbrack contains the function values at the bracketing triple
%   options may contain the following fields:
%     Display (default 'none'): set to 'iter' to see incremental results
%     TolFun: if present, also uses a fractional tolerance on function value,
%       not just on x
%     iter_max (default 100): the maximum number of iterations
% and
%   x is the optimum parameter value
%   fx is the function value at the minimum.
%
% See also: FMINBND.

% Original code: Philip D. Loewen
% 26 Feb 98 - original
% 30 Apr 01 - tuning
%  This is just typed from Numerical Recipes, modified for Matlab.
% 2010: wrap in its own function + a few modifications, Timothy E. Holy

% TODO: consider setting it so that it uses parabolic interpolation on the
% first iteration by default?

if (nargin < 4)
  options = struct;
end
options = default(options,'iter_max',100,'TolX',1e-6,'Display','none');
checkval = isfield(options,'TolFun');

% Initializations...
a = xbrack(1);  % Left bracket point
b = xbrack(3);  % Right bracket point
v = xbrack(2);  % Middle bracket point
x = v;   % Point with best fcn value so far
w = v;   % Point with second-best fcn value so far
         % (Later, v stores the previous value of w.)
e = 0.0; % Distance moved on the step before last.

fx = fxbrack(2);  % Remember fcn value at middle pt of bracket from before.
fv = fx;
fw = fx;

Zeps = eps^0.75;
gold   = (3.0+sqrt(5.0))/2.0;   % Golden ratio expansion factor

for iter=1:options.iter_max         % Main program loop.
  %fprintf('a %g, x %g, v %g, w %g, b %g\n',a,x,v,w,b);
  xm = 0.5*(a+b);  % Midpoint of current bracket
  tol1 = options.TolX*abs(x)+Zeps;
  tol2 = 2*tol1;

  % Test for done here.
  if ( abs(x-xm) <= (tol2-0.5*(b-a)) ) || (checkval && iter > 2 && (fw-fx)+(fv-fx) < options.TolFun*(abs(fv)+abs(fx)+abs(fw)))
    return
  end

  if ( abs(e) <= tol1),
    % Step before last was very small, so
    % take a Golden Section Step:
    if (x >= xm),
      e = a-x;
    else
      e = b-x;
    end;
    d = e/gold;
  else
    % Step before last was of decent size, so 
    % construct a trial parabolic fit.
    r = (x-w)*(fx-fv);
    q = (x-v)*(fx-fw);
    p = (x-v)*q - (x-w)*r;
    q = 2.0*(q-r);
    if (q>0.0), p = -p; end;
    q = abs(q);
    etemp = e;
    e = d;

    % Test viability of trial fit.
    if (abs(p)>= abs(0.5*q*etemp)) || (p <= q*(a-x)) || ( p>= q*(b-x)) 
      % Parabolic fit is poor, so take golden section step.
      if (x >= xm),
        e = a-x;
      else
        e = b-x;
      end;
      d = e/gold;
    else
      % Parabolic fit is OK, so use it.
      d = p/q;
      u = x+d;
      if (u-a < tol2) || (b-u < tol2), d = tol1*sign(xm-x); end;
    end;
  end;

  % Arrive here with increment  d  in hand.
  % Use d to form new x-value, insisting on moving at least tol1.
  if (abs(d) >= tol1),
    u = x+d;
  else
    u = x + tol1*sign(d);
  end;

  % Evaluate given function at point u
  % (This is the one function evaluation per iteration.)
  fu = f(u);
  if strcmp(options.Display,'iter')
    fprintf('u %g, fu %g\n',u,fu);
  end
  if (fu <= fx),
    % New evaluation point  u  is better than best-so-far  x
    if (u >= x), a=x; else b=x; end;   % Contract bracketing interval
    v=w; fv=fw;
    w=x; fw=fx;
    x=u; fx=fu;
  else
    % New evaluation point  u  is worse than best-so-far  x
    if (u < x), a=u; else b=u; end;
    if (fu <= fw) || (w==x),
      v = w; fv = fw;
      w = u; fw = fu;
    elseif(fu<=fv) || (v==x) || (v==w),
      v=u; fv=fu;
    end;
  end;

end;    % End of main program loop

error('Maximum iterations exceeded')
