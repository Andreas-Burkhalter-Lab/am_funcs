function [u,lambda,action] = mpl_opt(u,lambda,A,V,n,options)
% MPL_OPT: find Maximum Penalized Likelihood solution
%
% This function solves the discrete equations (equivalently, minimizes
% the discrete action) for Maximum Penalized Likelihood in one or more
% dimensions.  The efficient way to do this depends on the
% dimensionality.  In d = 1, the best variables to work in are psi, and
% the algorithm is O(N).  In d > 1, the only "exact" choice is to work in
% terms of a, and the algorithm is O(N^2).  This function handles either
% formulation, using an option.
%
% Syntax:
%   [ufinal,lambda,action] = mpl_opt(u,lambda,A,V,n,options)
% The meaning of u, A, and V depends upon the value of options.use_a:
% When options.use_a = 0 (choose if d = 1),
%  u = psi
%  A = L
%  V = volume operator in psi
% When options.use_a = 1 (choose if d > 1),
%  u = a
%  A = G (= inv(L))
%  V = volume operator in a.
%
% u is the best starting guess (e.g., through kernel density);
% lambda is the best starting guess for lambda (e.g., N);
% n is the multiplicity of each data point;
% options may have the following fields in addition to use_a:
%   tol (default 1e-8): determines the convergence criterion, applied to
%     the action.
%
% See also: MPL_INITIAL_GUESS, MPL_OPTW.
  
% Copyright 2006 by Timothy E. Holy.
  
  if (nargin < 6)
    error('Must supply all 6 arguments (and set use_a)');
  end
  if ~isfield(options,'use_a')
    error(['I don''t know which formulation to choose: set ' ...
           'options.use_a']);
  end
  if ~isfield(options,'tol')
    tol = 1e-8;
    if options.use_a
      tol = 1e-4;  % This version is not fully 2nd order, so relax our requirements
    end
  else
    tol = options.tol;
  end
  
  safe_slow = 0;
  
  N = sum(n);
  itermax = 20+4*round(sqrt(N));
  iter = 0;
  
  % Current normalization
  if iscell(V)
    psinorm = u'*tridiag_mult(V{1},V{2},V{3},u);
  else
    psinorm = u'*V*u;
  end
  u = u/sqrt(psinorm); % Start out normalized
  psinorm = 1;
  action1 = mpl_action(u,lambda,n,A,options);
  isdone = 0;
  while ~isdone
    % Calculate the new a or psi
    % Here we approximate the hessian as lambda*L, throwing out the
    % n./psi.^2 part. The justification is that any positive-definite
    % matrix times the gradient is a descent direction, and this way we
    % don't have to invert a full matrix.
    % In fact, this is necessary only when using the a-formulation, we
    % can handle the n./psi.^2 part when using the psi-formulation with
    % tridiagonal L.
    action0 = action1;
    u0 = u;
    lambda0 = lambda;
    if options.use_a
      psi0 = A*u0;
      psi_inv = n'./psi0;
      a2 = psi_inv/lambda;  % This is the new a if lambda weren't changing
      alpha = (2*(u0'*V*a2) + psinorm - 1)/(2*psinorm);
      u = (1-alpha)*u0+a2;  % This is the new a
      u1 = a2;
    else
      psi0 = u0;
      psi_inv = n'./psi0/lambda0;
      % This next is the new psi if we didn't have to worry about
      % normalization
      %psi = tridiag_inv(A{1},A{2},A{3},psi_inv);
      % Calculate with the full Hessian, since that's simple in this case
      % This will give faster convergence
      Lpsi0 = tridiag_mult(A{1},A{2},A{3},psi0);
      psi1 = tridiag_inv(A{1},A{2}+psi_inv./psi0,A{3},Lpsi0);
      psi2 = tridiag_inv(A{1},A{2}+psi_inv./psi0,A{3},psi_inv);
      pVp1 = psi0'*tridiag_mult(V{1},V{2},V{3},psi1);
      pVp2 = psi0'*tridiag_mult(V{1},V{2},V{3},psi2);
      alpha = (2*pVp2 + psinorm - 1)/(2*pVp1);
      u = psi0 + psi2 - alpha*psi1;
      % The following would have been the new psi if lambda were not also
      % changing. We'll test this one to "insure" global convergence
      u1 = psi0 + psi2 - psi1;
    end
    lambda = alpha*lambda0;
    if safe_slow
      if iscell(A)
        pLp = u1'*tridiag_mult(A{1},A{2},A{3},u1);
        psinorm = u1'*tridiag_mult(V{1},V{2},V{3},u1);
      else
        pLp = u1'*A*u1;
        psinorm = u1'*V*u1;
      end
      u1 = u1/sqrt(psinorm);
      lambda = N*psinorm/pLp;
    end
    % See if this new u1 results in a decreased action
    % Because of the change in lambda, u is not guaranteed to decrease
    % the action for a "good" step, so test only u1.
    action1 = mpl_action(u1,lambda0,n,A,options);
    da = action1-action0;
    dl = lambda-lambda0;
    if (action1 >= action0)
      % No, it didn't decrease the action, so we have to take a smaller
      % step. Let's find the minimum along the line between the current
      % position and the attempted step
      %[step_opt,action1] = fminbnd(...
      %    @(alpha) mpl_line_action(alpha,u0,u,lambda,n,A,options),...
      %    0,1);
      [step_opt,action1] = mpllinesearch(...
        @(alpha) mpl_line_action(alpha,u0,u1,lambda0,lambda0,n,A,options), ...
          0.5,action0);
      if (step_opt)
        % Now apply this step
        u = (1-step_opt)*u0 + step_opt*u1;
        %lambda = (1-step_opt)*lambda0 + step_opt*lambda;
        lambda = lambda0;   % Keep the same lambda
      else
        % There wasn't a good step. We must be within roundoff error of the
        % minimum. Keep the position but update lambda.
        u = u0;
      end
    end
    action1 = mpl_action(u,lambda,n,A,options);
    iter = iter+1;
    %[abs(action0-action1) abs(lambda - lambda0) da dl]
    isdone = ~((abs(da) > tol*N*abs(action0) || ...
               abs(dl) > tol*N) && ...
               iter < itermax);   
    if iscell(V)
      psinorm = u'*tridiag_mult(V{1},V{2},V{3},u);
    else
      psinorm = u'*V*u;
    end
    %u = u/sqrt(psinorm);
    %psinorm = 1;
  end
  if (iter == itermax)
    error('Failed to converge');
  end
  %iter
  
  
function action = mpl_action(u,lambda,n,A,options)
  if options.use_a
    psi = A*u;
    a = u;
  else
    psi = u;
    a = tridiag_mult(A{1},A{2},A{3},psi);
  end
  action = lambda*(psi'*a)-2*(n*log(abs(psi))) - lambda;

function action = mpl_line_action(alpha,u1,u2,lambda1,lambda2,n,A,options)
  u = (1-alpha)*u1 + alpha*u2;
  lambda = (1-alpha)*lambda1 + alpha*lambda2;
  action = mpl_action(u,lambda,n,A,options);

