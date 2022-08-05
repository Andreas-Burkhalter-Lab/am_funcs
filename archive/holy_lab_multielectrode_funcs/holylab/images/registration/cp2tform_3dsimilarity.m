function T = cp2tform_3dsimilarity(input_points,base_points)
% CP2TFORM_3DSIMILARITY: compute a three-dimensional similarity transform from control points
%
% Given a "standard" set of points (e.g., atlas coordinates) and "observed"
% positions of those points in another coordinate system (e.g., measured
% coordinates), compute the transform T that converts (using tformfwd)
% standard coordinates to observed coordinates. The transform is
% constrained to be a similarity transformation. This is like
% cp2tform(input,base,'similarity') except for 3-dimensional inputs.
%
% Syntax:
%   T = cp2tform_3dsimilarity(input_points,base_points)
% where
%   input_points is a 3-by-n_points matrix, each column containing one
%     "observed" point
%   base_points is a 3-by-n_points matrix, each column containing one
%     "standard" point
% and
%   T is a tform (see MAKETFORM) which converts standard points to observed
%     points when using TFORMFWD.
%
% See also: TFORMFWD, MAKETFORM, CP2TFORM, cp2tform_3d_affine.

% Copyright 2009 by Illya Tolokh and Timothy E. Holy

if (size(input_points,1) ~= 3 || size(base_points,1) ~= 3)
  error('Inputs must be matrices with each 3-dimensional point on a column')
end
n_points = size(base_points,2);
mu_base = mean(base_points,2);
mu_input = mean(input_points,2);

sigma_base2=sum(var(base_points,1,2));
%sigma_input2=sum(var(input_points,1,2));

sigma_base_input = zeros(3,3);
for i=1:n_points
  sigma_base_input = sigma_base_input+(input_points(:,i)-mu_input)*(base_points(:,i)-mu_base)';
end
sigma_base_input = sigma_base_input/n_points;

if det(sigma_base_input)>=0
  S = eye(3,3);
else
  S = diag([1,1,-1], 0);
end;
  
[U, D, V] = svd(sigma_base_input);

Rcalc = U*S*V';   %calculated rotation
ccalc = 1/sigma_base2*trace(D*S);  % scaling factor
t = mu_input-ccalc*Rcalc*mu_base; % translation

% Combine these into a transform
A = [ccalc*Rcalc'; t'];
T = maketform('affine',A);
%recalcY = ccalc*Rcalc*base_points + t*ones(1,n_points)