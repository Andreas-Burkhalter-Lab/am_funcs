function ny = crossmoveR(x,y,R,usey)
% CROSSMOVER: move a set of representative points to center of mass
% Given two sets of points, x (contains the raw data, N points) and y
% (q points, say q = min(N,100)), move the points y by the following
% algorithm: for a given y_i, find all points x within a radius R, and
% then move y_i to the center of mass of these x.  (Note: y_i is also
% inserted as one of these points, although this can be controlled by an
% option.)
% This algorithm runs in time dNq, so for q << N it avoids the N^2
% behavior of the general point movement algorithm.
% Syntax:
%   newy = crossmoveR(x,y,R)
%   newy = crossmoveR(x,y,R,usey)
% Note R may be either a scalar, or a vector if you want to assign a
% different radius to each point in y.
% 
% The "usey" option controls whether you throw y_i into the collection of
% points within a radius R. Default: true. (If there are no x's within R,
% then y_i will be unchanged in any event.)
%
% See also: MEANSHIFT, RMEANS.
  
% Copyright 2004 by Timothy E. Holy
  
  if (nargin < 4)
    usey = 1;
  end
  [d,q] = size(y);
  [d,N] = size(x);
  ny = zeros(size(y));
  if (length(R) == 1)
    R = repmat(R,1,q);
  end
  for i = 1:q
    %dist = sum((x - repmat(y(:,i),1,N)).^2,1);
    dist = mindist(x,y(:,i));
    indx = find(dist < R(i)^2);
    if ~isempty(indx)
      % Mean of all close x values & y_i
      tmp = sum(x(:,indx),2);
      if usey
        %ny(:,i) = (tmp+y(:,i))/(length(indx)+1);
        nx = length(indx);
        dist = mindist(y,y(:,i));
        indx = find(dist < R(i)^2);
        ny(:,i) = (tmp+sum(y(:,indx),2))/(nx + length(indx));
      else
        ny(:,i) = tmp/length(indx);
      end
    else
      ny(:,i) = y(:,i);
    end
  end