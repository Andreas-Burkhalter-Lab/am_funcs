function z = gammafunc(x, q)
%GAMMAFUNC Used by GAMMAFIT.
%	GAMMAFUNC(q) returns the error between the data and the
%	values computed by the Gamma function.
%  GAMMAFUNC assumes a function of the form
%
%	  y = q1 + (  (q(2)*t)^q(3) - exp( q(4)*t ) dt)
%
%	thus q(1) is the vertical ofset, q(2) is amplitude, q(3) is alpha, and q(4) is tau, q(5) is exponent
%	The data is in columns such that Data(:,1) is abscissa 
%	Data(:,2) is ordinate
%	The value of err is the sum squared error between data and function
%	given the parameters q.

z = q(1) + q(2)*(( q(3)*(x - q(4)) ).^q(5)).*(exp(-(q(3)*(x - q(4) )))/ ((q(5)^q(5))*exp(-q(5)))   ) ;

return;

