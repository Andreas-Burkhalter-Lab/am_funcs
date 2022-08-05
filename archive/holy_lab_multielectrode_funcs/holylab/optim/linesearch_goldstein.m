function [a,fa] = linesearch_goldstein(f,f0,fp0,a,c)
% LINESEARCH_GOLDSTEIN: one-dimensional minimization for Newton-type descent
% 
% Syntax:
%   [a,fa] = linesearch_goldstein(f,f0,fp0,a0,c)
% where
%   f is a one-dimensional objective function, f(a)
%   f0 is the value at a = 0
%   fp0 is the derivative at a = 0
%   a0 (default 1) is an initial guess for a
%   c (default 0.25) is a scalar, 0 < c < 0.5
% and
%   a is an acceptable step size
%   fa is f(a).
%
% This performs minimization constrained by the de-linked Goldstein
% conditions, as described in:
%   B. Christianson, "Global Convergence using De-linked Goldstein or
%     Wolfe Linesearch Conditions." AMO - Advanced Modeling and
%     Optimization, Volume 11, Number 1, 2009.

% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 5)
    c = 0.25;
  end
  if (nargin < 4)
    a = 1;
  end
  R = 2;
  b = a;
  fa = f(a);
  fb = fa;
  while (fa > f0 + c*a*fp0)
    b = a;
    fb = fa;
    a = b/R;
    fa = f(a);
  end
  while (fb < f0 + (1-c)*b*fp0)
    a = b;
    fa = fb;
    b = R*a;
    fb = f(b);
  end
