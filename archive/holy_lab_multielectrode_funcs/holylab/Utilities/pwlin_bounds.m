function yi = pwlin_bounds(x,y,xi)
% PWLIN_BOUNDS: implement a piecewise linear monotonic function, respecting bounds
%
% Syntax:
%   yi = pwlin_bounds(x,y,xi)
% where
%   the piecewise linear function is defined by the "knots" in the vector x
%     and the values in the vector y;
%   the function is to be evaluated at the positions xi.
%
% See also: INTERP1.

  yi = interp1(x,y,xi,'linear','extrap');
  mxy = max(y);
  mny = min(y);
  yi(yi > mxy) = mxy;
  yi(yi < mny) = mny;
