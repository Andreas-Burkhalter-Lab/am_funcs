function [Dist,w,tout,rout,k,D]=dtw(t,r,options)
%DTW: Dynamic Time Warping algorithm
%Syntax:
%   [dist,w]=dtw(t,r)
%   [dist,w,tout,rout]=dtw(t,r)
%   [dist,w,tout,rout,k,D]=dtw(t,r)
%   [...]=dtw(t,r,options)
%where
%  t and r are the two vectors you're comparing
%and
%  dist is the unnormalized square distance between t and r (sum of the
%    square difference between aligned vectors)
%  w is an noverlap-by-2 vector giving the optimum path
%  tout and rout are the two aligned vectors (for the overlap region
%    only)
%  k is the path length
%  D is the accumulated distance matrix
%
% options is a structure which may have the following fields:
%   nanpenalty: the penalty added when only one vector has an NaN
% Note: you can pad the beginning and ending with an NaN, and then this
% will allow penalty-free translation in the time domain. Use this with
% some care, as it can cause alignment by shifting one beyond the end of
% the other.  See DTWFE for a better way.
%
% See also: DTWFE.

% Downloaded from Matlab Central, written originally by Timothy Felty
% Modifications 2005 by Timothy E. Holy:
%   Extra outputs to return the aligned vectors
%   Added options
%   Added NaN-padding

if (nargin < 3)
  options = struct;
end
if ~isfield(options,'odpenalty')
  options.odpenalty = 1;
end
if ~isfield(options,'nanpenalty')
  options.nanpenalty = 0;
end

[rows,N]=size(t);
[rows,M]=size(r);
if (rows == 1)
  % Adjust for NaNs---set the distance between an NaN and any number to be
  % zero, so that NaNs can expand to fill gaps, etc.
  tindx = find(~isnan(t));
  rindx = find(~isnan(r));
  ttmp = t(tindx);
  rtmp = r(rindx);
  n = length(tindx);
  m = length(rindx);
  dtmp = zeros(n,m);
  for i=1:n
    for j=1:m
      dtmp(i,j) = (ttmp(i)-rtmp(j))^2;
    end
  end
  d = zeros(N,M) + options.nanpenalty;
  d(isnan(t),isnan(r)) = 0;  % When they're both NaNs, no penalty
  d(tindx,rindx) = dtmp;
else
  % This doesn't yet have NaN padding
  for n=1:N
    for m=1:M
      d(n,m)=sum((t(:,n)-r(:,m)).^2);
    end
  end
end

D=zeros(size(d));
D(1,1)=d(1,1);

for n=2:N
  D(n,1)=d(n,1)+D(n-1,1);
end
for m=2:M
  D(1,m)=d(1,m)+D(1,m-1);
end
for n=2:N
  for m=2:M
    D(n,m)=d(n,m)+min([D(n-1,m),...
                       D(n-1,m-1),...
                       D(n,m-1)]);
  end
end

Dist=D(N,M);
n=N;
m=M;
k=1;
w=[];
w(1,:)=[N,M];
while ((n+m)~=2)
  if (n-1)==0
    m=m-1;
  elseif (m-1)==0
    n=n-1;
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
  k=k+1;
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

