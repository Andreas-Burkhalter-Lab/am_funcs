function z = gauss2Dfunc(x, y, q)
%GAUSS2DFUNC Used by GAUSS2DERR.
%  GAUSS2DFUNC assumes a function of the form
%
%     z = q(1) + q(2) * exp((-0.5*((x - q(3))/ q(4)).^2) + (-0.5*((y - q(5))/ q(6)).^2));
%	thus q(1) is base rate, q(2) is amplitude, q(3) is xcenter, q(4) is xwidth, q(5) is ycenter, and q(6) is ywidth

z = q(1) + q(2) * exp((-0.5*((x - q(3))/ q(4)).^2) + (-0.5*((y - q(5))/ q(6)).^2));


return;
