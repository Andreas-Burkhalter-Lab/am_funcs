function P = perpproj(v)
% PERPPROJ: calculate the matrix to project perpendicular to a vector
% Syntax:
%   P = perpproj(v)
% where v is the perpendicular
%
% Example:
% To calculate the component of x which is perpendicular to v, you can say
%   xP = P*x
% where
%   P = perpproj(v)
% Note this is equivalent to saying
%   xP = x - (x dot v)/(v dot v) * v
  
  vnorm = v/sqrt(sum(v.*v));  % Normalize v
  vnorm = vnorm(:);
  P = eye(length(vnorm)) - vnorm * vnorm';
  