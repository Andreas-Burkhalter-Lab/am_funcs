function [LLout,gradLLout] = llcv1(w,x,n,options)
% LLCV1: log-likelihood and gradient with respect to 1d projection
% For continuous kernel-density estimator, compute the negative of
% log-likelihood and its first derivative with respect to coordinates of
% projection.
%
% This optionally uses a cross-validated log-likelihoood. With
% cross-validation, you may want to first consolidate points whose
% projections are very close. See CONSOLIDATE_PROJ.
%
% [LLout,gradLLout] = llcv1(w,x,n,options)
% where
%   w is a 1-by-d vector, the direction for projection
%   x is a d-by-M matrix of points
%   n is a 1-by-M vector of multiplicities
%   options is a structure with the following fields:
%     CV (default false): if true, uses cross-validation
% and, when not using cross-validation,
%   LLout = LL = - n*log((L1^-1)(n'/N))   (note the -N*log(w*w')/2 term
%              is not included; add that manually if needed.)
%       (Note this is really the negative log-likelihood!)
%   gradLLout(i) = d_(dw_i)S
% and for cross-validated log-likelihoods
%   LLout = LLCV = - n*log((L1^-1 - 1/2)(n'/N)) (")
%   gradLLout(i) = d_(dw_i)SCV
% Note here that L1 is the nu=1 gradient penalty+normalization, i.e. the
% inverse matrix to G_ij = exp(-|X_i - X_j|)/2.
%
% See also: GRADWLLCV, CONSOLIDATE_PROJ.

% Possible further options?
%   sorted
%   c, s, & cholL1 precalculated?

if (nargin < 4)
  options = struct;
end
if ~isfield(options,'CV')
  options.CV = 0;
end

N = sum(n);
dim = length(w);
M = length(n);

% Project the data points & sort the projections
X = (w*x)';
[X,sort_order] = sort(X);
n = n(sort_order);
dX = diff(X);
npN = n'/N;
permstate = spparms;
spparms('autommd',0);  % This speeds things up later

c = cosh(dX);
s = sinh(dX);
index_inf = find(c == Inf);   % In cases where points are widely separated
% First the value
d0 = c./s;
d0(index_inf) = 1;
B = [[-1./s;0],[1;d0]+[d0;1],[0;-1./s]];
L1 = spdiags(B,[-1 0 1],M,M); % a slow step. Consider MEX solver for P?
cL1 = chol(L1);
P = cL1\(cL1'\npN);
if options.CV
  PCV = P - npN/2;
  LLout = -(n*log(PCV));
else
  LLout = -(n*log(P));
end
if (nargout < 2)
  spparms(permstate);
  return
end
% Now the gradient
o = c./s.^2;
o(index_inf) = 0;
d = -1./s.^2;
dx = diff(x(:,sort_order),1,2);   % In high dimensions, this & the next line are the slow steps
G = dx*dLvec(d,o,P);
if options.CV
  nP = cL1\(cL1'\(n'./PCV));
else
  nP = cL1\(cL1'\(n'./P));
end
gradLLout = G*nP;
spparms(permstate);

function G = dLvec(d,o,vec)
ov0 = o.*vec(1:end-1);
ov1 = o.*vec(2:end);
dv0 = d.*vec(1:end-1);
dv1 = d.*vec(2:end);
B = [ov1+dv0,ov0+dv1];
G = spdiags(B,[0 1],length(d),length(vec));
