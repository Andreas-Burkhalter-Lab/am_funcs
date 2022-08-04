function err = weibull_func(q)
%WEIBULL_FUNC Used by WEIBULL_FIT.
%	returns the error between the data and the
%	values computed by the current function of weibul params.
%	assumes a function of the form
%
%	  y = 1 - .5 * exp( -(x/q(1))^q(2) )
%
%	thus q(1) is alpha, and q(2) is beta.
%	The data is in columns such that Data(:,1) is abscissa 
%	Data(:,2) is observed percent correct (0..1)
%	Data(:,3) is number of observations.
%	The value of err is the -log likelihood of obtaining Data
%	given the parameters q.
%
global Data

TINY = 1e-10;
x = Data(:,1);
y = Data(:,2);
n = Data(:,3);
z = 1 - .5 * exp( -(x/q(1)).^q(2) );
z = z - (z > .999999)*TINY + (z < .0000001)*TINY; 
llik = n .* y .* log(z) +  n .* (1-y) .* log(1-z);
% for moment, just return the square error
% err = norm(z-y);
err = -sum(llik);

