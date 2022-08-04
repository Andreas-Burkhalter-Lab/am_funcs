function z = gaborgauss_sumfunc(x, q)
%GABORGAUSS_SUMFUNC Used by GABORGAUSSSUM_FIT.
%  GABORGAUSS_SUMFUNC assumes a function of the form
%
%	  y = q(1) + q(2) exp(-(2*(x - q(3))/q(4))^2 )*cos(q(5)*(x - q(3))  ) + q(7) * exp( -((x - q(3) )/q(4)).^2 )
%	thus q(1) is base rate, q(2) is amplitude of the gabor, q(3) is shared center of the gaussian and gabor, 
%   and q(4) is size (sigma) of the gabor
%	q(5) is the frequency
%   q(6) is the amplitude of the gaussian and 
%   q(7) is the exponent of the gabor
% q(8) is the sigma of the gaussian

%z = q(1) + q(2)*exp(-0.5*(abs(x - q(3))/ q(4)).^q(7)).*cos(2*pi*q(5)*(x - q(3)) ) + q(6) * exp(-0.5*(abs(x - q(3))/ q(8)).^2);
%abs needed in event of non-even exponent q(8)

%q(3) is in ms
%q(5) is in cycles/sec
z = q(1) + q(2)*exp( -(  abs(x - q(3)) / q(4)).^q(7)).*cos(2*pi*q(5)/1000 * (x - q(3)) ) + q(6) * exp( -(abs(x - q(3) )/ q(8) ).^2);

%try to constrain exponent between 0.5 and 5
return;
