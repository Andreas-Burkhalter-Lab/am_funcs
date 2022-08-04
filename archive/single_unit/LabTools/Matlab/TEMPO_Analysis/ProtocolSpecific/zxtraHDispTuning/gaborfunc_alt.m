function z = gaborfunc_alt(x, q)
%GABORFUNC_ALT:  SLightly different parameterization of GABORFUNC
%  assumes a function of the form
%
%	  y = q(1)*(1 + q(2) exp(-(2*(x - q(3))/q(4))^2 )*cos(2*pi*q(5)*(x - q(3)) + q(6)))
%
%	thus ***q(1) is overall scaling factor***, ***q(2) is (amplitude/base rate)***, 
%   q(3) is center, and q(4) is size (sigma)
%	q(5) is the frequency, and q(6) is the phase

z = q(1) * ( 1 + q(2)*exp(-0.5*((x - q(3))/ q(4)).^2).*cos(2*pi*q(5)*(x - q(3)) + q(6)) );

return;
