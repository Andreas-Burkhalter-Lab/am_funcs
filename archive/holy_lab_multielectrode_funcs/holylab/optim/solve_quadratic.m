function x = solve_quadratic(a,b,c,s)
% solve_quadratic: numerically-stable solution to quadratic equation
%
% Syntax:
%   x = solve_quadratic(a,b,c)
%   x = solve_quadratic(a,b,c,s)
% where
%   a, b, and c are the coefficients of the second, first, and zeroth
%     powers of x, respectively.  These can be row vectors of length N if
%     you want to solve multiple equations.
%   s (optional) specifies which solution you want, positive if you want
%     the solution
%         -b + sqrt(b^2-4*a*c)
%         --------------------
%                 2*a
%     and negative if you want the solution
%         -b - sqrt(b^2-4*a*c)
%         --------------------
%                 2*a
% and
%   x is either a 2-by-N matrix (if s is absent) or a 1-by-N vector. If
%     you do not supply s, the top row is the +sqrt solution, and the
%     bottom row is the -sqrt solution.
%
% The main virtue of this function is a numerically-stable solution,
% i.e., no delicate cancelation between b and the sqrt term.

% Copyright 2010 by Timothy E. Holy
  
  N = length(a);
  signb = b >= 0;
  sqrtterm = sqrt(b.^2 - 4*a.*c);
  if (nargin > 3)
    % We want just one solution to the quadratic
    if (s == 0)
      error('The sign (s) cannot be zero');
    end
    if (s > 0)
      x = sq_pos(a,b,c,signb,sqrtterm);
    else
      x = sq_neg(a,b,c,signb,sqrtterm);
    end
  else
    x = zeros(2,N);
    x(1,:) = sq_pos(a,b,c,signb,sqrtterm);
    x(2,:) = sq_neg(a,b,c,signb,sqrtterm);
  end
  
function x = sq_pos(a,b,c,signb,sqrtterm)
  nsb = ~signb;
  x(nsb) = (-b(nsb)+sqrtterm(nsb))/(2*a(nsb));
  x(signb) = -2*c(signb) ./ (b(signb) + sqrtterm(signb));

function x = sq_neg(a,b,c,signb,sqrtterm)
  nsb = ~signb;
  x(signb) = (-b(signb)-sqrtterm(signb))/(2*a(signb));
  x(nsb) = -2*c(nsb) ./ (b(nsb) - sqrtterm(nsb));
  