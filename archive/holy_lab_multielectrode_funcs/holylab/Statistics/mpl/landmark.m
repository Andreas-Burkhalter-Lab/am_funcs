function [xlandmark,n_landmark,landmark_assignment] = landmark(x,num_landmarks)
% LANDMARK: pick landmarks to reduce the size of a data set
% From the data, a random selection of points is drawn. Then each point in
% the original data set is assigned to the closest landmark. The landmarks,
% together with their (new) multiplicities and the assignments, are
% returned.
%
% Syntax:
%   [xlandmark,n_landmark,landmark_assignment] = landmark(x,num_landmarks)
% where
%   x is a d-by-N matrix, each column contains one data point;
%   num_landmarks is the desired number of landmarks;
% and
%   xlandmark is a d-by-K matrix (K = min(num_landmarks,N)) containing the
%     positions of the landmarks;
%   n_landmark is a 1-by-K vector containing the number of points assigned
%     to each landmark;
%   landmark_assignments is a 1-by-N vector indicating the landmark # that
%     each original point is assigned to.

% Copyright 2006 by Timothy E. Holy

  npts_orig = size(x,2);
  num_landmarks_real = min(num_landmarks,npts_orig);
  if (num_landmarks_real < npts_orig)
    rp = randperm(npts_orig);
    landmark_index = rp(1:num_landmarks_real);
    xlandmark = x(:,landmark_index);
    dland = sqrdist(xlandmark,x);
    [sdland,sort_order] = sort(dland);
    landmark_assignment = sort_order(1,:);
    n_landmark = zeros(1,num_landmarks_real);
    for i = 1:npts_orig
      n_landmark(landmark_assignment(i)) = ...
        n_landmark(landmark_assignment(i)) + 1;
    end
  else
    xlandmark = x;
    n_landmark = ones(1,npts_orig);
    landmark_assignment = 1:npts_orig;
  end