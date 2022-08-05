function [x,fx] = linbrack(f,x,fx,options)
% LINBRACK: find a bracketing triple for 1-dimensional minimization
%
% Syntax:
%  [x,fx] = linbrack(f,x0)
%  [x,fx] = linbrack(f,x0,fx0)
%  [x,fx] = linbrack(f,x0,fx0,options)
% where
%   f is the 1-dimensional function to be minimized
%   x0 is a vector of starting points at which you want to evaluate the
%     function (must have at least two points)
%   fx0 (optional) is a partial or complete list of corresponding function
%     values (by supplying these, you can prevent re-evaluation of function
%     values you already know)
% and
%   x is an ordered bracketing triple [a b c]
%   fx contains the function values [f(a) f(b) f(c)]
%
% The bracketing triple will have finite function values, even if the
% function to be optimized sometimes returns Inf or NaN.
%
% Example: create and save the following test function:
%   function fx = testfun(x)
%     fx = (exp(-x)-12).^2;
%     fx(x < -3) = inf;  % just to make it harder
%
%   Then find a bracketing triple this way:
%   [x,fx] = linbrack(@testfun,[0 1]);
%
% See also: LINBRACKGD, FMINBND.

% Copyright 2010 by Timothy E. Holy

  if (length(x) < 2)
    error('Must supply at least two evaluation points');
  end
  if (nargin < 3)
    fx = [];
  end
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'iter_max',50);
  iter = 0;
  
  % Evaluate at all the user-specified points
  while (length(fx) < length(x))
    n = length(fx)+1;
    fx(n) = f(x(n));
  end
  
  % Find the lowest point tested so far
  [fa,index] = min(fx);
  a = x(index);
  % Find neighbors of this point
  [xs,sortOrder] = sort(x);
  index = find(xs == a);
  if (index > 1 && index < length(x))
    % We have already found a triple, quit
    index = index-1:index+1; 
    x = xs(index);
    fx = fx(sortOrder(index));
    if ~any(isinf(fx) | isnan(fx))
      return
    end
  elseif (index == 1)
    index = 1:min(3,length(xs));
  else
    index = max(1,length(xs)-2):length(xs);
  end
  x = xs(index);
  fx = fx(sortOrder(index));

  while (iter < options.iter_max)
    [x,fx] = get_triple(f,x,fx);
    if (fx(2) < fx(1) && fx(2) < fx(3))
      return
    end
    % Consider parabolic extrapolation, but only if we have evidence of
    % convexity
    dx = diff(x);
    if ~any(isinf(fx) | isnan(fx)) && fx(2) < (dx(1)*fx(3) + dx(2)*fx(1))/(dx(1)+dx(2))
      dfx = diff(fx);
      r = -dx(1)*dfx(2);
      q = -dx(2)*dfx(1);
      u = x(2) + sign(q-r)*(q*dx(2)+r*dx(1))/(2*max(abs(q-r),eps(x(2))));
      if (u-x(3) > x(2)-x(1) || u-x(1) < x(2)-x(3) || (u > x(1) && u < x(3)))
        fu = f(u);
        if (isnan(fu) || isinf(fu))
%           warning('lineminimization:brackparabolic','parabolic extrapolation during bracketing gave Inf or NaN');
        end
        if (fu < min(fx))
          % We have a new minimum, keep it
          iter = iter+1;
          if (u < x(1))
            x = [u x(1:2)];
            fx = [fu fx(1:2)];
          elseif (u > x(3))
            x = [x(2:3) u];
            fx = [fx(2:3) fu];
          else
            % This is an interior point, we are done
            if (u < x(2))
              x = [x(1) u x(2)];
              fx = [fx(1) fu fx(2)];
            else
              x = [x(2) u x(3)];
              fx = [fx(2) fu fx(3)];
            end
            return
          end
          continue
        end
      end
    end
    % Get a new point by interval magnification
    if (fx(1)==fx(2) || fx(2)==fx(3))
      x = x([1 3]);
      fx = fx([1 3]);
    elseif (fx(1) < fx(3))
      x = x([1 2]);
      fx = fx([1 2]);
    else
      x = x([2 3]);
      fx = fx([2 3]);
    end
    iter = iter+1;
  end
  error('Failed to bracket minimum. Increase iter_max if necessary.');
end
  
function [x,fx] = get_triple(f,x,fx)
  % x & fx will either be 2-vectors or 3-vectors; if they are 2-vectors, we
  % need a 3rd point
  if (length(x) < 3)
    if (fx(1) < fx(2))
      c = x(1) - 1.6*diff(x);
      x = [c x];
      fx = [f(c) fx];
    else
      c = x(2) + 1.6*diff(x);
      x = [x c];
      fx = [fx f(c)];
    end
  end
  % Check 3-vector for NaNs & Infs
  badFlag = isinf(fx) | isnan(fx);
  while (badFlag(1))
    c = (x(1)+x(2))/2;
    x = [c x(2:3)];
    fx = [f(c) fx(2:3)];
    badFlag = isinf(fx) | isnan(fx);
  end
  while (badFlag(3))
    c = (x(2)+x(3))/2;
    x = [x(1:2) c];
    fx = [fx(1:2) f(c)];
    badFlag = isinf(fx) | isnan(fx);
  end
end
