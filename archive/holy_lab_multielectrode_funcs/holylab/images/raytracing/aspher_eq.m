function [y,normal] = aspher_eq(x,params)
% ASPHER_EQ: standard equation for an aspheric lens surface
% Syntax:
%   [y,normal] = aspher_eq(x,params)
% where
%   (x,y) is a point on the surface
%   params is a structure describing the surface shape with the following
%     fields:
%       curv: lens curvature (= 1/R for spherical lens)
%       k: ellipticity (k = 0, A-D = 0 is a spherical lens)
%       A: x^4 coefficient
%       B: x^6 coefficient
%       C: x^8 coefficient
%       D: x^10 coefficient
%
% The equation for the surface is
%    y = (1 - sqrt(1-(k+1)*curv^2*x^2))/(curv*(k+1)) +
%        A x^4 + B x^6 + C x^8 + D x^10

% Equation from the Thorlabs website. Note the first term in the equation
% above is sensitive to roundoff error when x is small, so in practice we
% calculate with the stabilized version (mult by
% (1+sqrtterm)/(1+sqrtterm)).
  sqrtterm = sqrt(1-(params.k+1)*params.curv^2*x^2);
  y = params.curv*x^2/(1 + sqrtterm) + ...
      params.A*x^4 + params.B*x^6 + params.C*x^8 + params.D*x^10;
  if (nargout > 1)
    % Compute the normal
    %fp = (2*params.curv*x*(1+sqrtterm) - (params.k+1)*params.curv^3*x^3)/...
    %     ((1+sqrtterm)^2*sqrtterm) + ...
    fp = params.curv*x/sqrtterm + ...
         4*params.A*x^3 + 6*params.B*x^5 + 8*params.C*x^7 + ...
         10*params.D*x^9;
    %fprintf('x %f, fp %f\n',x,fp);
    normal = [-fp 1]/sqrt(1+fp^2);
  end
  