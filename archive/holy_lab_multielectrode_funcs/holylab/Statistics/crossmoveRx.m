function ny = crossmoveRx(x,y,R)
% CROSSMOVERX: move landmarks to weighted center of mass
% Given two sets of points, x (contains the raw data in columns, N
% points) and y (q points), move the points y up the density gradient
% defined by x.  Unlike CROSSMOVER, each data point is  by the following
% algorithm: for a given y_i, find all points x within a radius R, and
% then move y_i to the center of mass of these x.  (Note: y_i is also
% inserted as one of these points, although this can be controlled by an
% option.)
% This algorithm runs in time dNq, so for q << N it avoids the N^2
% behavior of the general point movement algorithm.
% Syntax:
%   newy = crossmoveRx(x,y,R)
% Note R may be either a scalar, or a vector if you want to assign a
% different radius to each point in y.
% 
% See also: MEANSHIFT, RMEANS.
  
% Copyright 2004 by Timothy E. Holy
  
  [d,q] = size(y);
  [d,N] = size(x);
  if (length(R) < N)
    error('crossmoveRx requires one radius for each x');
  end
  R2 = R.^2;
  ny = zeros(size(y));
  for i = 1:q
    %dist = sum((x - repmat(y(:,i),1,N)).^2,1);
    dist = mindist(x,y(:,i));
    indx = find(dist < R2);
    if ~isempty(indx)
      % Mean of all close x values, weighted by the radius
      % In high dimensions this could be delicate, so regularize using
      % the smallest contributing R.
      Rmin = min(R(indx));
      Rw = 1./(R(indx)/Rmin).^(d+2);
      ny(:,i) = sum(x(:,indx).*repmat(Rw,d,1),2)/sum(Rw);
    else
      ny(:,i) = y(:,i);
    end
  end