function chi2expected = avgchi2(psi,n,lambda,A,G0,ops)
% AVGCHI2: expected value of chi-squared, from the data alone
% This computes the expected value of chi-squared, <chi^2>.  The distribution of
% chi-squared is then given by the conventional result for <chi^2> degrees
% of freedom,
%     P(c) = const * c^(<chi^2>/2-1) exp(-c/2).
% The Laplace transform of this distribution is
%            (1+2*alpha)^(-<chi^2>/2).
%
% Syntax:
%   chi2expected = avgchi2(psi,n,lambda,L,G(0),struct('use_a',0))
% OR
%   chi2expected = avgchi2(psi,n,lambda,G,G(0),struct('use_a',1))
% where
%   psi (M-by-1) is the maximum penalized likelihood solution;
%   n (1-by-M) contains the multiplicity of data points;
%   lambda (scalar) contains the Lagrange multiplier for normalization;
%   G or L contains the kernel matrix or its inverse, respectively;
%   G(0) is the value of the kernel at 0;
%   The last input contains the guide to which form, G or L, is being
%     supplied;
% and
%   chi2expected is <chi^2>.
%
% The psi, lambda, and L/G inputs can be conveniently obtained from the
% third output argument of MPL.
%
% For the sparse d=1 formulation (use_a = 0), this function is O(M).
% Otherwise (use_a = 1), it is O(M^3).
%
% See also: MPL.

% A is L (for use_a == 0) or G (for use_a == 1)
  M = length(n);
  Dvec = n'./psi.^2;
  Dmtrx = spdiags(Dvec,0,M,M);
  if ops.use_a
    GD = A*Dmtrx;
    lGD = lambda*eye(N) + GD;
    chi2a = 2*trace(lGD\GD);  % Very slow...O(M^3)
    u = lGD\(A*psi);
  else
    % Don't do it in the following way:
    %chi2a = 2*trace((lambda*A + Dmtrx)\Dmtrx);
    % because (lambda*L+Dmtrx)\Dmtrx is full. Instead, go back to the
    % expression for \tilde P(\alpha) and compute the
    % derivative with respect to alpha numerically.
    chi2guess = G0*sum(Dvec)/lambda;  % order of magnitude, so we know how big to make the finite difference
    dalpha = 1e-6/chi2guess;
    % The technical issue is a way of computing the determinants
    % without overflowing:
    %Pdalpha = sqrt(det((lambda*L+Dmtrx)/...
    %                   (lambda*L + (1+4*dalpha)*Dmtrx)));
    if iscell(A)
      % This will likely be the slow step; if performance is an issue,
      % one can write a direct tridiagonal Cholesky decomposer as a MEX file.
      A = spdiags([[A{1};0],A{2},[0;A{3}]],[-1 0 1],M,M);
    end
    R1 = chol(lambda*A+Dmtrx);  % chol is faster & will also give use our sqrt
    R2 = chol(lambda*A + (1+4*dalpha)*Dmtrx);
    eigrat = spdiags(R1,0)./spdiags(R2,0); % ratio of diagonals
    Pdalpha = prod(eigrat);
    chi2a = (1-Pdalpha)/dalpha;  % Here's the finite difference
    u = (lambda*A+Dmtrx)\psi;  % get ready for chi2b
  end
  chi2b = 2*(u'*Dmtrx*u)/(psi'*u); % Can replace with 1
  chi2expected = chi2a-chi2b;
