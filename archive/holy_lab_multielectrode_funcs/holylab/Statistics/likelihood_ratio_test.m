function p = likelihood_ratio_test(n,q)
  % likelihood_ratio_test: evaluate multinomial (binning) goodness-of-fit
  %
  % For binning, the well-known chi-squared statistic is only an
  % approximate measure of goodness-of-fit; particularly for cases where
  % the bins are of different probabilities, it tends to be too strict,
  % rating random draws from the multinomial distribution as being far less
  % likely than they actually are. The likelihood ratio test, where the
  % statistic is twice the Kullback-Leibler divergence, can be used to
  % develop a more accurate approximation. This function implements the
  % techniques recommended in
  %   P.J. Smith, D. S. Rae, R. W. Manderscheid, and S. Silberg,
  %   "Approximating the moments and distribution of the likelihood ratio
  %   statistic for multinomial goodness of fit." J. Am. Statistical
  %   Assoc. 76: 737-740 (1981).
  % These approximations tend to be OK once there are 30 or more total
  % points distributed among all the bins.
  %
  % For high accuracy, you are still best off if your bins are not of
  % wildly differing probability. You may want to use adaptive binning
  % (see, e.g., bin_evenly).
  %
  % Syntax:
  %   p = likelihood_ratio_test(n,q)
  % where
  %   n is a 1-by-k vector giving the counts per bin
  %   q is a 1-by-k vector giving the model probability per bin
  %     (sum(q) = 1)
  % and on output
  %   p is the probability of observing this, or lower, degree of
  %     disagreement between n/sum(n) and q.
  %
  % See also: bin_evenly.
  
  % Copyright 2012 by Timothy E. Holy
  
  k = length(n);
  N = sum(n);
  switch k
    case 1
      error('Not defined for a single bin');
    case 2
      % For two bins, Smith et al argue that no approximation is
      % acceptable. Use the exact binomial distribution.
      p = binocdf(max(n),sum(n),max(q));
    case 3
      % For three bins, Smith et al argue that the naive chi-squared
      % actually gives the best results. You tend to need N=30 or larger
      % for this to be "good" even for the equiprobability case.
      chi2 = sum((n - N*q).^2./(N*q));
      p = chi2cdf(chi2,k-1);
    case {4,5,6}
      % For 4-6 bins, Smith et al recommend the "G/q' procedure"
      qp = 1 + (sum(1./q) - 1 + sum(1./q-1./q.^2)/N)/(6*N*(k-1));
      keep = n>0;
      G2 = 2*sum(n(keep).*log(n(keep)./(N*q(keep))));
      p = chi2cdf(G2/qp,k-1);
    otherwise
      % Use the "approximate beta fit" technique
      sqim1 = sum(1./q)-1;
      sqq = sum(1./q-1./q.^2);
      M = -2*N*log(min(q));
      E = k-1+sqim1/(6*N) + sqq/(6*N^2);
      V = 2*(k-1)+sqim1/(3*N)+4*sqq/(3*N^2);
      a = (E/M/V)*(E*(M-E)-V);
      b = (M-E)*a/E;
      keep = n>0;
      G2 = 2*sum(n(keep).*log(n(keep)./(N*q(keep))));
      p = betainc(G2/M,a,b);
  end
end
      
      