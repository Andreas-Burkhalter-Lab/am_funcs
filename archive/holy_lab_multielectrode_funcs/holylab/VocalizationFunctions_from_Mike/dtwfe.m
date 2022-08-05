function [dist,w,tout,rout,ellmin,D,ell]=dtwfe(t,r,options)
%DTWFE: Dynamic Time Warping with Free Ends
%This takes two vectors as input, and calculates a good alignment
%between them.  The time axis on either is allowed to stretch or shift to
%achieve alignment.  Compared to the conventional DTW algorithm,
%this allows the best alignment path to start and/or end in the interior
%of one of the vectors, so that one can be a shifted version of the
%other.  One does not have to be entirely contained in the other, unless
%this behavior is specified in the options.
%
%The alignment chosen is the one with the smallest square error per unit
%length, with some adjustment available through options (e.g.,
%minoverlap & favoroverlap).  The path length is defined as the shorter of the two
%vectors in the overlap region.
  
%Note this is not guaranteed to find the globally smallest square error
%per unit length, but in practice it seems to do pretty well.
%Error/length is not convex, so it can't easily be found by dynamic
%programming.  A good counterexample might be if the error in one half is
%very different from the error in the other half---you might prefer to
%have a somewhat longer path in the first half if that's the low-error
%portion, just to balance the second half, even if that longer path in
%the first half is not the minimum error/length.
%
%Syntax:
%   [dist,w]=dtwfe(t,r)
%   [dist,w,tout,rout]=dtwfe(t,r)
%   [dist,w,tout,rout,ellmin,D,ell]=dtwfe(t,r)
%   [...]=dtwfe(t,r,options)
%where
%  t and r are the two vectors you're comparing
%and
%  dist is the unnormalized square distance between t and r (sum of the
%    square difference between aligned vectors)
%  w is an noverlap-by-2 vector giving the optimum path
%  tout and rout are the two aligned vectors (for the overlap region
%    only)
%  ellmin is the path length, the minimum of the numbers of unique
%    elements of each vector included in the overlap region
%  D is the accumulated distance matrix
%  ell is the accumulated path length matrix
%
%options is a structure which may have the following fields:
%  minoverlap: consider the best path of length >= minoverlap
%    (default: 1)
%  ellpower: weight the choice of best path by ell^ellpower,
%    to encourage longer paths to be used. In particular, ellpower = 2
%    corresponds to weighting the best path by (length(t)+length(r))/ell,
%    i.e., extrapolating the error into regions of non-overlap.  (default: 1)
%  nanpenalty: if you want to have some undefined regions of your signal,
%    you can put in one or more NaNs. This options determines the penalty
%    applied when an NaN is aligned with a finite number. (default: Inf)
%
%See also: DTW.

% Copyright 2005 by Timothy E. Holy

if (nargin < 3)
  options = struct;
end
if ~isfield(options,'minoverlap')
  options.minoverlap = 1;
end
if ~isfield(options,'nanpenalty')
  options.nanpenalty = Inf;
end
if ~isfield(options,'ellpower')
  options.ellpower = 1;
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
ellt = zeros(size(d));
D(:,1) = d(:,1);
D(1,:) = d(1,:);
ellt(:,1) = 1;
ellt(1,:) = 1;
ellr = ellt;
tinc = [1 1 0];
rinc = [1 0 1];
for n=2:N
  for m=2:M
    Dnext = [D(n-1,m-1) D(n-1,m) D(n,m-1)];
    ellchoices = [ellt(n-1,m-1) ellt(n-1,m) ellt(n,m-1);...
                  ellr(n-1,m-1) ellr(n-1,m) ellr(n,m-1)];
    ellnext = min(ellchoices) + 1;
    [D_over_l_min,indx_min] = min(Dnext./ellnext);
    D(n,m)=d(n,m)+Dnext(indx_min);
    ellt(n,m) = ellchoices(1,indx_min) + tinc(indx_min);
    ellr(n,m) = ellchoices(2,indx_min) + rinc(indx_min);
  end
end

% Find the minimum D/ell at the ending edges, subject to the minimum
% overlap criterion
% ell is defined as min(ellt,ellr)
ell = min(ellt,ellr);
% First, find the best along the last column
indx_ellc = find(ell(:,end) >= options.minoverlap);
[distc,ind_distc] = min(D(indx_ellc,end)./ell(indx_ellc,end).^options.ellpower);
ind_distc = indx_ellc(ind_distc);
% Now the best along the last row
indx_ellr = find(ell(end,:) >= options.minoverlap);
[distr,ind_distr] = min(D(end,indx_ellr)./ell(end,indx_ellr).^options.ellpower);
ind_distr = indx_ellr(ind_distr);
% Find the best overall
if (isempty(distr) | (distc < distr))
  n = ind_distc;
  m = M;
else
  n = N;
  m = ind_distr;
end
dist = D(n,m);
ellmin = ell(n,m);
% Now trace back
w=[];
w(1,:)=[n,m];
if (isempty(n) | isempty(m))
  tout = NaN;
  rout = NaN;
  ellmin = NaN;
  return;
end
while ((n > 1) && (m > 1))
  Dback = [D(n-1,m),D(n,m-1),D(n-1,m-1)];
  ellback = [ell(n-1,m),ell(n,m-1),ell(n-1,m-1)];
  [values,number]=min(Dback./ellback);
  switch number
   case 1
    n=n-1;
   case 2
    m=m-1;
   case 3
    n=n-1;
    m=m-1;
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

% if (nargin > 2)
%   tmp1 = t(w(end:-1:1,1));
%   tmp2 = r(w(end:-1:1,2));
% 
%   if options.plot
%     mplot = [tmp1(:),tmp2(:)];
%     plot(mplot)
%   end
% 
%   if options.normalize
%     dist = dist/sqrt(sum(tmp1.^2)*sum(tmp2.^2));
%   end
% 
%   if options.angle
%     dist = sum(tmp1.*tmp2)/sqrt(sum(tmp1.^2)*sum(tmp2.^2));
%   end
% end
