function [u,lambda] = mpl_initial_guess(A,V,n,options)
% MPL_INITIAL_GUESS: provide kernel-density estimate for initial guess
% Syntax:
%   [u,lambda] = mpl_initial_guess(A,V,n,options)
% where the meaning of u, A, and V depend upon whether you're in d = 1
% (and can use the O(N) theory) or in d > 1 (and have to use the O(N^2)
% theory).  See MPL_OPT for an explanation of these parameters.
%
% See also: MPL_OPT.

% Copyright 2006 by Timothy E. Holy
    
  if options.use_a
    P = A*n';
    psi = sqrt(P);
    a = n'./psi;  % this is true at optimum (except normalization)
    nrm = a'*V*a; % ... and here's the normalization
    a = a/sqrt(nrm);
    u = a;
  else
    %P = A\n';
    P = tridiag_inv(A{1},A{2},A{3},n');
    psi = sqrt(P);
    %nrm = psi'*V*psi;
    nrm = psi'*tridiag_mult(V{1},V{2},V{3},psi);
    psi = psi/sqrt(nrm);
    u = psi;
  end
  lambda = sum(n);
  %lambda = sum(n)/(u'*A*u);  % no better in practice?
  