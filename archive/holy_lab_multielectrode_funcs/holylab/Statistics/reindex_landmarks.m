function out = reindex_landmarks(in,index)
% REINDEX_LANDMARKS: update landmark structure when points are resampled
  %
% This function is often used for bootstrap resampling of the original
% data points.  The original data points might be described by a d-by-N
% matrix x, with each column of x a different observation.  You might
% want to generate a bootstrap dataset y = x(:,index), where index is a
% 1-by-M vector (perhaps M < N, perhaps M = N) listing the points
% contributing to the "new" dataset.
%
% Syntax:
%   lm_new = reindex_landmarks(lm_old,index)
% where
%   lm_old is the original landmark structure
%   index is a 1-by-M vector listing the data points that are going to be
%     the "new" data points.  (There can be duplicates, it can be shorter
%     than original list, etc.)
% and
%   lm_new is a landmark structure appropriate for the points
%     y = x(:,index).
%
% See also: CHOOSE_LANDMARKS.
  
% Copyright 2008 by Timothy E. Holy
  
  out = in;
  out.landmarkAssignment = in.landmarkAssignment(index);
  out.dist_to_closest_landmark = in.dist_to_closest_landmark(index);
  out.landmarkList = agglabel(out.landmarkAssignment);
  if length(out.landmarkList) < length(in.landmarkList)
    out.landmarkList{length(in.landmarkList)} = [];
  end
  
  