function [dist,w,tout,rout,D]=dtwfe(t,r,options)
%DTWFE: Dynamic Time Warping with one free end
%This takes two vectors as input, and calculates a good alignment
%between them.  The time axis on either is allowed to stretch or shift to
%achieve alignment.  Compared to the conventional DTW algorithm,
%this allows the best alignment path to start and/or end in the interior
%of the first vector; the second must be entirely contained within the first.
%
%The alignment chosen is the one with the smallest square error.
%
%Syntax:
%   [dist,w]=dtw1fe(t,r)
%   [dist,w,tout,rout]=dtw1fe(t,r)
%   [dist,w,tout,rout,D]=dtw1fe(t,r)
%   [...]=dtwfe(t,r,options)
%where
%  t and r are the two vectors you're comparing---the second (r) must be
%    contained entirely within the first (t, "template") 
%and
%  dist is the unnormalized square distance between t and r (sum of the
%    square difference between aligned vectors)
%  w is an noverlap-by-2 vector giving the optimum path
%  tout and rout are the two aligned vectors (for the overlap region
%    only)
%  D is the accumulated distance matrix
%
%options is a structure which may have the following fields:
%  nanpenalty: if you want to have some undefined regions of your signal,
%    you can put in one or more NaNs. This options determines the penalty
%    applied when an NaN is aligned with a finite number. (default: Inf)
%  ndpenalty: non-diagonal penalty, a penalty for stretching either vector.
%    Expressed as an extra contribution to the square difference. (default:
%    0)
%
%See also: DTW, DTWFE.

% Copyright 2005 by Timothy E. Holy

if (nargin < 3)
  options = struct;
end
if ~isfield(options,'nanpenalty')
  options.nanpenalty = Inf;
end
if ~isfield(options,'ndpenalty')
  options.ndpenalty = 0;
end

[rows,N]=size(t);
[rows,M]=size(r);
if (rows == 1)
  % Adjust for NaNs---NaNs match other NaNs exactly, and match finite
  % numbers with an optional penalty.
  tindx_nan = isnan(t);
  rindx_nan = isnan(r);
  ttmp = t(~tindx_nan);
  rtmp = r(~rindx_nan);
  n = length(ttmp);
  m = length(rtmp);
  dtmp = (repmat(ttmp',1,m) - repmat(rtmp,n,1)).^2;
  d = zeros(N,M) + options.nanpenalty;
  d(tindx_nan,rindx_nan) = 0;  % When they're both NaNs, no penalty
  d(~tindx_nan,~rindx_nan) = dtmp; % When both finite, use distance
else
  % This doesn't yet have NaN padding
  for n=1:N
    for m=1:M
      d(n,m)=sum((t(:,n)-r(:,m)).^2);
    end
  end
end

% Calculate the cumulative distance and pathlength matrices
D=zeros(size(d));
D(:,1) = d(:,1);  % First vector can start anywhere
for m = 2:M
  D(1,m) = d(1,m) + D(1,m-1);  % Second vector must start at beginning
end
for n=2:N
  for m=2:M
    D(n,m)=d(n,m)+min([D(n-1,m) + options.ndpenalty,...
                       D(n-1,m-1),...
                       D(n,m-1) + options.ndpenalty]);
  end
end

% Find the enpoint of the best path
[dist,n] = min(D(:,M));
% Traceback
m = M;
w=[];
w(1,:)=[n,m];
while (m > 1)
  if (n == 1)
    m = m-1;
  else
    [values,number]=min([D(n-1,m),D(n,m-1),D(n-1,m-1)]);
    switch number
     case 1
      n=n-1;
     case 2
      m=m-1;
     case 3
      n=n-1;
      m=m-1;
    end
  end
  w=cat(1,w,[n,m]);
end
w = w(end:-1:1,:);  % To put w in increasing order
if (nargout > 2)
  if (rows > 1)
    tout = t(:,w(:,1));
    rout = r(:,w(:,2));
  else
    tout = t(w(:,1));
    rout = r(w(:,2));
  end
end
