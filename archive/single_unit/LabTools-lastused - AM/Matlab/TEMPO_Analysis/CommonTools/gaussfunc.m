function z = gaussfunc(x, q)
%GAUSSFUNC Used by GAUSSFIT.
%  GAUSSFUNC assumes a function of the form
%
%	  y = q(1) + q(2) * exp(-0.5*((x - q(3))/q(4)).^2 )
%
%	thus q(1) is base rate, q(2) is amplitude, q(3) is center, and q(4) is size

z = q(1) + q(2) * exp(-0.5*((x - q(3))/ q(4)).^2);

return;
