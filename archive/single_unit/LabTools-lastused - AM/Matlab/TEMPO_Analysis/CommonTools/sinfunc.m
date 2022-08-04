function z = sinfunc(x, q)
%SINFUNC used by SINFIT
%  SINFUNC assumes a function of the form
%
%	  y = q(1) * sin(q(2)*x + q(3)) + q(4)
%
%	thus q(1) is amplitude, q(2) is frequency, q(3) is phase, and q(4) is base line, q(5) is 

z = q(1) * sin(q(2)*x + q(3)) + q(4);

return;
