function curve_data = inv_weibull_curve(x,q)
%	returns a vector of values that give the Weibull function
% 	sampled at a range of values given by x but inverted
%	and assumes a function of the form - BJP 1/5/01
%
%	  y = 1.5 - (1 - .5 * exp( -(x/( q(1) ) )^q(2) ))
%
%	thus q(1) is alpha, and q(2) is beta.

curve_data = 1.5 - (1 - .5 * exp( -(x/(q(1) )).^q(2) ));
