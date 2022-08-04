function curve_data = weibull_bs_curve(x,q)
%	returns a vector of values that give the Weibull function
% 	sampled at a range of values given by x
%	assumes a function of the form
%
%	  y = 1 - (1-q(3)) * exp( -(x/q(1))^q(2) ) 
%
%	thus q(1) is alpha, and q(2) is beta.

curve_data = 1 - (1-q(3)) * exp( -(x/q(1)).^q(2) );

