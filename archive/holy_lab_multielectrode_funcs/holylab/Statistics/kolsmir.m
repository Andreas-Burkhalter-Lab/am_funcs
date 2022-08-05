function [d,pd] = kolsmir(x1,x2,op)
% KOLSMIR: kolmogorov-smirnov statistics for 2 data sets
% [d,pd] = kolsmir(x1,x2)
% where
%   x1, x2 are vectors of measured values
% and
%   d is the K-S statistic, max(cdf1(x)-cdf2(x))
%   pd is Prob(d > observed)

% Copyright 2006 Timothy E. Holy

if nargin<3
    op = struct;
end
op = default(op,'output_NaN_instead_of_errors',0);

%[sx,sindex] = sort([x1(:);x2(:)]);
if (length(x1) < 2 || length(x2) < 2)
  d = NaN;
  pd = NaN;
  return;
end
ux = unique([x1(:);x2(:)]);
pdf1 = zeros(size(ux));
pdf2 = zeros(size(ux));
indx1 = findainb(x1,ux);
for i = 1:length(indx1)
  pdf1(indx1(i)) = pdf1(indx1(i)) + 1/length(x1);  % OK with repeat values
end
indx2 = findainb(x2,ux);
for i = 1:length(indx2)
  pdf2(indx2(i)) = pdf2(indx2(i)) + 1/length(x2);
end
%pdf2(find(sindex > length(x1))) = 1/length(x2);
cdf1 = cumsum(pdf1);
cdf2 = cumsum(pdf2);
d = max(abs(cdf1-cdf2));
if (nargout > 1)
  % Compute significance
  if (d == 0)
    pd = 1;
    return;
  end
  n = length(x1)*length(x2)/(length(x1)+length(x2));
  lam = (sqrt(n) + .12 + .11/sqrt(n))*d;
  lam2 = -2*lam^2;
  pd = 0;
  fac = 2;
  oldterm = 0;
  for j = 1:100
    term = fac*exp(lam2*j^2);
    pd = pd+term;
    if (abs(term) <= .001*oldterm | abs(term) <= 1e-8*pd)
      return;
    end
    fac = -fac;
    oldterm = abs(term);
  end
  if op.output_NaN_instead_of_errors
      d = NaN;
      p = NaN;
  else
    error('Computation of significance did not converge');
  end
end