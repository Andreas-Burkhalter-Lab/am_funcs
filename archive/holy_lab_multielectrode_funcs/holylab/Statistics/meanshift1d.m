function xout = meanshift1d(xin,thresh)
% MEANSHIFT1D: cluster points in one dimension by mean-shift algorithm
% Given an input set of points xin, each point is iteratively moved to the
% center of mass of all points within a given distance of that point.
% When points converge to a final position, the algorithm terminates.
%
% Syntax:
%   xout = meanshift1d(xin,thresh)
% where
%   xin is a vector containing the initial set of locations;
%   thresh is the distance threshold to define the "local neighborhood";
% and
%   xout is a vector containing the final locations; the point starting
%     at xin(j) moves to xout(j).
%
% This uses the full point-movement algorithm; it's O(N^2) if your
% threshold is large enough that much of the sample is encompased, but
% smaller thresholds result in O(N) behavior.
%
% See also: RMEANS.
  
% Copyright 2005 by Timothy E. Holy
  npts = length(xin);
  [sx,sxindx] = sort(xin);
  sx = sx(:)';
  sxo = zeros(size(sx));
  doneflag = 0;
  % Iterate until points have stopped moving
  while ~doneflag
    sxo = meanshift1d1(sx,thresh);
    %doneflag = all(sx == sxo);
    doneflag = max(abs(sx-sxo)) < eps*sx(end)*100;
    %fprintf('%g\n',max(abs(sx-sxo)))
    sx = sxo;
  end
  % Deal with roundoff problems
  ibreak = [1;find(diff(sx(:)) > thresh)+1;length(sx)+1];
  for i = 1:length(ibreak)-1
    sx(ibreak(i):ibreak(i+1)-1) = sx(ibreak(i));
  end
  % Return points in the original order
  [ssxi,invindx] = sort(sxindx);
  xout = sx(invindx);
  xout = reshape(xout,size(xin));
  