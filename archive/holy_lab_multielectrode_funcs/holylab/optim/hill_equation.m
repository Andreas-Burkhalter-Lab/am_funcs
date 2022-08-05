function rpred = hill_equation(s,C)
% hill_equation: calculated predicted firing rates from parameters and concentrations
% Syntax:
%   rpred = hill_equation(s,C)
% where
%   s is a structure of the form returned by fit_hill_equation, in
%     results.full:
%        rmax: the maximum rate
%        K: a n_components-by-1 vector of binding coefficients
%        n: the Hill coefficient
%        r0: the offset
%   C is a n_samples-by-n_components matrix, containing the concentration
%     of each component in each sample
% and
%   rpred is a 1-by-n_samples vector, containing the predicted rate for
%     each sample.
%
% See also: fit_hill_equation.

% Copyright 2011 by Timothy E. Holy

  ceff = sum(bsxfun(@rdivide,C,s.K),2);
  ceffn = ceff.^s.n;
  ceffn1 = ceffn ./ (1+ceffn);
  rpred = s.r0 + s.rmax * ceffn1;
  