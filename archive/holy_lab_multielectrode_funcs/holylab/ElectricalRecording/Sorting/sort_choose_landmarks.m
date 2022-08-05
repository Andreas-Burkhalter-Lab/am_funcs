function [landmarkIndex,landmarkPosition] = sort_choose_landmarks(x, ...
						  n_landmarks)
% SORT_CHOOSE_LANDMARKS: a spike-sorting wrapper for choose_landmarks
% Syntax:
%   [landmarkIndex,landmarkPosition] = sort_choose_landmarks(x,n_landmarks)
% where
%   x is a ndims-by-npoints matrix of observations, each observation a
%     column of x;
%   n_landmarks is the number of desired landmarks;
% and
%   landmarkIndex is a 1-by-npoints vector, landmarkIndex(i) being the
%     index of the landmark closest to the ith data point x(:,i);
%   landmarkPosition is a ndims-by-nlandmarks matrix,
%     landmarkPosition(:,i) being the location of the ith landmark.
%
% See also: CHOOSE_LANDMARKS.

% Copyright 2006 by Timothy E. Holy

  lminfo = choose_landmarks(x,n_landmarks);
  landmarkIndex = lminfo.landmarkAssignment;
  landmarkPosition = lminfo.landmarks;
  