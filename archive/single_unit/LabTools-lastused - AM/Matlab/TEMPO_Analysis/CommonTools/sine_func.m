function z = sine_func(x, q)
%SINE_FUNC Used by SINE_FIT.
%  SINE_FUNC assumes a function of the form
%
%	  y = q(1) + q(2)*cos(2*pi*q(3)*x + q(4) )
%
%	thus q(1) is base rate, q(2) is amplitude, q(3) is the frequency, and q(4) is the phase

z = q(1) + q(2)*cos(2*pi*q(3)*x + q(4) );

return;
