function [x,lsqerr] = lsq_monotonic(b,sigma)
% LSQ_MONOTONIC: calculate least square error, constrained by monotonicity
% Given a vector b of target values, and another vector sigma of
% uncertainties, calculate the vector x that minimizes
%     sum( (x-b).^2./sigma.^2 )
% subject to the constraint that the values in x are monotonically
% increasing.
%
% Syntax:
%    [x,lsqerr] = lsq_monotonic(b)
%    [x,lsqerr] = lsq_monotonic(b,sigma)
% The first variant assumes all sigma = 1.
%
% The algorithm is a global one, potentially O(N^2) in execution speed.
% It proceeds by induction. Given a value X for which all x_j
% up to j=m satisfy x_j < X.  The error up to the mth term is piecewise
% quadratic. When adding the m+1 term, combine the quadratics
% and compute the minimum for the quadratic in each region of
% piecewise-quadratic support, starting with the leftmost region. If the
% minimum occurs
%   (a) at a point higher than the right boundary of the region: advance to
%       the next region
%   (b) within the region: define a new boundary at the minimum and clear
%       all higher boundaries.
% When all terms have been included, the least-square error is the value at
% the rightmost boundary. The values of x can be read out from the region
% boundaries.
  
% Copyright 2009 by Timothy E. Holy
% Note: there may be a faster algorithm in terms of cumsum(coefbase) and
% min_loc. Look for places where min_loc decreases.
  
  N = length(b);
  if (nargin < 2)
    sigma = ones(1,N);
  end
  % Make sure they are row vectors
  if (size(b,1) > 1)
    b = b';
  end
  if (size(sigma,1) > 1)
    sigma = sigma';
  end

  if any(sigma == 0)
    error('Some sigmas are zero, can''t optimize');
  end

  % Allocate the maximum storage, but keep track of the position of the current end
  % in the variable n. (This prevents a lot of memory re-allocation, which
  % would slow the algorithm.)
  coefp = zeros(3,N);
  breaks = zeros(1,N);
  breaksIndex = zeros(1,N);
  n = 1;
  % Initialize the boundaries and the quadratic coefficients within each
  % region (quadratic coefficients are [2 1 0] in terms of order).
  coefbase = [ones(1,N); b; b.^2] .* repmat(1./sigma.^2,3,1); % note lack of -2 in cross-term
  coefp(:,1) = coefbase(:,1);
  breaks(1) = b(1);
  breaksIndex(1) = 1;
  
  for i = 2:N
    j = 1;
    while (j <= n)
      coef_old = coefp(:,j);
      comb_coef = coefp(:,j) + coefbase(:,i);
      min_loc = comb_coef(2) / comb_coef(1);
      coefp(:,j) = comb_coef;  % combine the coefficients
      if (min_loc < breaks(j))
        % Min occured within the current region
        breaks(j) = min_loc;
        breaksIndex(j) = i;
        n = j;
      end
      j = j+1;
    end
    if (b(i) > breaks(n))
      % The current point is beyond any known regions. Create a new
      % region
      minval_prev = coef_old(1)*breaks(n).^2 - 2 * coef_old(2)*breaks(n) ...
        + coef_old(3);
      n = n+1;
      breaks(n) = b(i);
      coefp(:,n) = coefbase(:,i);
      coefp(3,n) = coefp(3,n)+minval_prev;
      breaksIndex(n) = i;
    end
  end
  % Now unpack the breaks to yield x
  breaksMap = zeros(1,N+1);
  breaksMap(breaksIndex(1:n)+1) = 1;
  breaksMap = cumsum(breaksMap)+1;
  x = breaks(breaksMap(1:N));
  lsqerr = coefp(1,n)*breaks(n)^2 - 2 * coefp(2,n)*breaks(n) ...
    + coefp(3,n);
end
