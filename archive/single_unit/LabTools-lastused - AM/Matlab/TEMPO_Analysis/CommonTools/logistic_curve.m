function curve_data = logistic_curve(x,q)
%	returns a vector of values that give the logistic function
%	sampled at a range of values given by x
%	assumes a function of the form
%
%	  y = 1 / (1 + exp( -(x-q(2))/q(1) ))
%
%	thus q(1) is alpha, and q(2) is beta.

curve_data = 1 ./ (1 + exp( -(x-q(2))./q(1) ));
