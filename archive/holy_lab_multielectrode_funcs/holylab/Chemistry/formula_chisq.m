function chi2 = formula_chisq(f,nI,nIcov)
% FORMULA_CHISQ: calculate goodness-of-fit for candidate formulas to isotopologue data
%
% Isotopologue abundance can yield estimates for the numbers of atoms of
% particular elements. This function takes a set of candidates for the
% molecular formula (often derived from the exact mass) and calculates the
% goodness-of-fit to the elemental multiplicity derived from the isotopologues.
%
% Syntax:
%    chi2 = formula_chisq(f,nI,nIcov)
% where
%    f is an n-by-n_elements matrix of integers, listing the numbers of
%      each type of element (usually as listed in isotopeinfo, see
%      ISOTOPES);
%    nI is a vector of estimated multiplicities (can contain NaNs for
%      elements lacking abundant isotopes)
%    nIcov is the covariance matrix for the fit of nI
% and
%    chisq is an n-by-1 vector containing the goodness of fit,
%       v*inv(nIcov)*v', where v = f(i,:)-nI for each row i of f.
%
% See also: ISOTOPES.

% Copyright 2009 by Timothy E. Holy

  [n,n_elements] = size(f);
  if length(nI) ~= n_elements
    error('Mismatch in size of formula and elemental multiplicity')
  end
  goodIndx = ~isnan(nI);
  chi2 = nan(n,1);
  inv_nIcov = inv(nIcov(goodIndx,goodIndx));
  nIg = nI(goodIndx);
  for i = 1:n
    v = f(i,goodIndx) - nIg;
    chi2(i) = v*inv_nIcov*v';
  end
  