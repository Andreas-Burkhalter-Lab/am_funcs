function [x,sigma] = raywaist(rf,xrange,options)
% RAYWAIST: calculate the width of a collection of rays at various
% positions
% Syntax:
%   [x,sigma] = raywaist(r,xrange)
%   [x,sigma] = raywaist(r,xrange,options)
% where
%   r is a structure array of rays (see RAY)
%   xrange is the range of deviation along the x-axis from the minimum (a
%     2-vector, e.g., [-2 5] will compute waist width between 2 units to
%     the left of the minimum and 5 units to the right of the minimum);
%   options is a structure with the following fields
%     add_diffraction: if true, broadens the waist additively by
%       0.61*lambda/NA (default: true).  See lambda, next field
%     lambda: the wavelength in the same units that all other dimensions
%       are using. This _must_ be supplied if add_diffraction is true.
%       Default: 0.000488.
% and
%   x is a set of x positions
%   sigma is the standard deviation at each position
  
  if (nargin < 3)
      options = struct;
  end
  if ~isfield(options,'add_diffraction')
      options.add_diffraction = 1;
  end
  if ~isfield(options,'lambda')
    options.lambda = .000488;
  end

  coef0 = 0;
  coef1 = 0;
  coef2 = 0;
  
  coef3 = 0; % To calculate y mean
  coef4 = 0; % To calculate y mean
  
  for i = 1:length(rf)
    erat = rf(i).e(2)/rf(i).e(1);
    coef0 = coef0 + (rf(i).x0(2)-erat*rf(i).x0(1))^2 * rf(i).I;
    coef1 = coef1 + 2 * (rf(i).x0(2)-erat*rf(i).x0(1)) * erat * rf(i).I;
    coef2 = coef2 + erat^2 * rf(i).I;
    
    coef3 = coef3 + rf(i).x0(2)*rf(i).I - rf(i).x0(1)*erat*rf(i).I;
    coef4 = coef4 + erat*rf(i).I;
    
  end
  sI = sum([rf.I]);
  coef0 = coef0/sI;
  coef1 = coef1/sI;
  coef2 = coef2/sI;
  coef3 = coef3/sI;
  coef4 = coef4/sI;
  
  xmin = -(coef1-2*coef3*coef4)/(2*(coef2-coef4.^2));
  x = linspace(xmin+xrange(1),xmin+xrange(2),101);
  
  ybar = coef3+coef4*x;
  
  sigma = (coef0 + coef1*x + coef2*x.^2) - ybar.^2;
  
  sigma = sqrt(sigma);
 
  % Add diffraction correction
  if (options.add_diffraction)
    e_all = [rf.e];
    sin_all = e_all(2:2:end);
    NA = (max(sin_all)-min(sin_all))/2;
    sigma = sigma + 0.61*options.lambda/NA;
  end
   