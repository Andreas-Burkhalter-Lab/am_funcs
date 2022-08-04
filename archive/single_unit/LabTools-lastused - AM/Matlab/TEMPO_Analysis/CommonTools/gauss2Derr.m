function err = gauss2Derr(q)
%GAUSS2DERR Used by GAUSS2DFIT.
%	GAUSSERR(q) returns the error between the data and the
%	values computed by the GAUSSFUNC function.

global Data RawData

x = RawData(:,1);
y = RawData(:,2);
spikes = RawData(:,3);

z = gauss2Dfunc(x, y, q);

%threshold the fitted values (don't allow less than zero)
z(z < 0) = 0;

% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.  GCD, 1/31/01
err = norm(sqrt(z)-sqrt(spikes));

return;