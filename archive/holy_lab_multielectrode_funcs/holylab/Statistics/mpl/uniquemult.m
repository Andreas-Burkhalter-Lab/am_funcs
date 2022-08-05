function [Xs,no,uindx] = uniquemult(X,ni,prec,noscale);
% UNIQUEMULT: find the unique points & their multiplicities
% [Xs,no,uindx] = uniquemult(X,ni,prec)
% where
%   X is the matrix of points (N by d)
%   ni is the input multiplicities of the points (default: all 1)
%   prec (optional) is the roundoff precision factor, relative to the
%     overall scale in the coordinate, for judging equality (default = 1e-8)
%   noscale is set to true if you want the precision factor to
%      be absolute, rather than relative to the range in that coordinate
%      (default: false)
%
% and
%   Xs is the sorted set of unique points
%   no is the output multiplicities
%   uindx is the index mapping, Xs = X(uindx,:) up to roundoff

  % One bit of trouble: roundoff error can make 2 "identical"
  % points differ by an insignificant amount. First pre-process
  % points to round off such error
  [N,d] = size(X);
  if (nargin < 2 | isempty(ni))
	  ni = ones(1,N);
  end
  if (nargin < 3 | isempty(prec))
    prec = 1e-8;
  end
  if (nargin < 4 | isempty(noscale))
	  noscale = 0;
  end
  if (prec > 0 & noscale == 0)
    Xmin = repmat(min(X),N,1);
    X2 = X - Xmin;
    Xfac = repmat(prec*max(X2),N,1);
    X2 = round(X2 ./ Xfac);
    Xrnd = X2 .* Xfac + Xmin;
  elseif (prec > 0)
    Xrnd = round(X/prec)*prec;
  else
    Xrnd = X;
  end
  [Xrs,pinv] = sortrows(Xrnd);
  nis = ni(pinv);
  % Now can make the call to unique
  [Xs,uindx,Xj] = unique(Xrs,'rows');
  N = size(X,1);
  M = size(Xs,1);
  no = zeros(1,M);
  for i = 1:N
    no(Xj(i)) = no(Xj(i))+nis(i);
  end
  uindx = pinv(uindx);
