function z = linearfunc(x, q)
%LINEARFUNC Used by LINEARFIT.
%  LINEARFUNC assumes a function of the form
%
%	  y = q(1) + q(2) * x
%
%	thus q(1) intercept and q(2) is the slope

z = q(1) + q(2) * x;

return;
